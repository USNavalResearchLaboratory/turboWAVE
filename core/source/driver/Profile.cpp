module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"
#ifndef USE_STD_MODULE
#include <algorithm>
#endif

export module driver:profile;
import :engine;

import input;
import functions;
import logger;

using namespace tw::bc;

export namespace tw
{
	namespace profile
	{
		enum class quantity { density,energy,px,py,pz };
		enum class timing { triggered,maintained };
		enum class loading { statistical,deterministic };
		enum class shape { triangle,sin2,quartic,quintic,sech };
	}
}

export struct Profile : Engine
{
	tw::profile::shape segmentShape;
	tw::grid::geometry symmetry;
	tw::vec3 centerPt;
	tw::vec3 modeNumber;
	tw::Float modeAmplitude;
	tw::basis orientation;
	tw::Float gammaBoost;

	// items needed for particle/fluid loading
private:
	tw::vec3 driftMomentum;
public:
	tw::profile::quantity whichQuantity;
	tw::vec3 thermalMomentum;
	tw::Float temperature;
	bool neutralize,variableCharge;
	tw::profile::loading loadingMethod;
	tw::profile::timing timingMethod;
	tw::Float t0,t1;
	bool wasTriggered;

	Profile(const std::string& name,MetricSpace *m,Task *tsk) : Engine(name,m,tsk) {
		symmetry = tw::grid::cartesian;
		segmentShape = tw::profile::shape::triangle;
		neutralize = true;
		timingMethod = tw::profile::timing::triggered;
		variableCharge = false;
		loadingMethod = tw::profile::loading::deterministic;
		whichQuantity = tw::profile::quantity::density;
		modeAmplitude = 0.0;
		modeNumber = 0.0;
		temperature = 0.0;
		wasTriggered = false;
		t0 = 0.0;
		t1 = tw::big_pos;
		orientation.u = tw::vec3(1,0,0);
		orientation.v = tw::vec3(0,1,0);
		orientation.w = tw::vec3(0,0,1);
		gammaBoost = 1.0;

		// Lots of enumerated-type hash tables involved in Profile directives.
		// Do them first, then setup directives.

		std::map<std::string,tw::profile::quantity> qty = {
			{ "density",tw::profile::quantity::density },
			{ "energy",tw::profile::quantity::energy },
			{ "px",tw::profile::quantity::px },
			{ "py",tw::profile::quantity::py },
			{ "pz",tw::profile::quantity::pz }};

		std::map<std::string,tw::profile::shape> shp = {
			{ "triangle",tw::profile::shape::triangle },
			{ "quartic",tw::profile::shape::quartic },
			{ "quintic",tw::profile::shape::quintic }};

		std::map<std::string,tw::profile::timing> tm = {
			{ "triggered",tw::profile::timing::triggered },
			{ "maintained",tw::profile::timing::maintained }};

		std::map<std::string,tw::profile::loading> ld = {
			{ "deterministic",tw::profile::loading::deterministic },
			{ "statistical",tw::profile::loading::statistical }};

		std::map<std::string,tw::grid::geometry> geo = {{"cylindrical",tw::grid::cylindrical},{"spherical",tw::grid::spherical}};

		directives.Add("neutralize",new tw::input::Bool(&neutralize),false); // not used
		directives.Add("position",new tw::input::Vec3(&centerPt),false);
		directives.Add("type",new tw::input::Enums<tw::profile::quantity>(qty,&whichQuantity),false);
		directives.Add("drift momentum",new tw::input::Vec3(&driftMomentum),false);
		directives.Add("thermal momentum",new tw::input::Vec3(&thermalMomentum),false);
		directives.Add("temperature",new tw::input::Float(&temperature),false);
		directives.Add("shape",new tw::input::Enums<tw::profile::shape>(shp,&segmentShape),false);
		directives.Add("timing",new tw::input::Enums<tw::profile::timing>(tm,&timingMethod),false);
		directives.Add("t0",new tw::input::Float(&t0),false);
		directives.Add("t1",new tw::input::Float(&t1),false);
		directives.Add("loading",new tw::input::Enums<tw::profile::loading>(ld,&loadingMethod),false);
		directives.Add("symmetry",new tw::input::Enums<tw::grid::geometry>(geo,&symmetry),false);
		directives.Add("mode amplitude",new tw::input::Float(&modeAmplitude),false);
		directives.Add("mode number",new tw::input::Vec3(&modeNumber),false);
		directives.Add("boosted frame gamma",new tw::input::Float(&gammaBoost),false);
		directives.Add("particle weight",new tw::input::Custom,false);
		directives.Add("euler angles",new tw::input::Custom,false);
	}
	virtual void Initialize() {;}
	tw::vec3 DriftMomentum(const tw::Float& mass) {
		tw::Float p2 = driftMomentum ^ driftMomentum;
		tw::Float p0 = std::sqrt(mass*mass + p2);
		tw::vec4 p4(p0,driftMomentum);
		p4.zBoost(gammaBoost,-1.0);
		return p4.spatial();
		// tw::vec4 v4(0.0,driftMomentum/mass);
		// tw::Float gb2 = v4 ^ v4;
		// v4[0] = std::sqrt(1.0 + gb2);
		// v4.zBoost(gammaBoost,-1.0);
		// return mass*v4.spatial();
	}
	tw::Float Temperature(const tw::Float& mass) {
		if (temperature!=0.0)
			return temperature;
		else
			return sqr(thermalMomentum.x)/mass; // appropriate for exp(-v^2/(2*vth^2)) convention
	}
	tw::vec3 Boost(const tw::vec3& pos)
	{
		// Here the boost only works if the profile is constant in time.
		// The function's caller is giving us boosted frame coordinates.
		// The user is giving us lab frame coordinates.
		// Therefore first transform arguments to lab frame, then proceed as usual.
		tw::vec4 x4(0.0,pos);
		x4.zBoost(gammaBoost,1.0);
		return x4.spatial();
	}
	tw::vec3 Translate_Rotate(const tw::vec3& pos)
	{
		// The argument should already be boosted.
		// N.b. the boost applies to both region and profile, while translate-rotate applies only to the profile.
		tw::vec3 p = pos - centerPt;
		orientation.ExpressInBasis(&p);
		return p;
	}
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds) {
		return theRgn->Inside(Boost(pos),ds) ? 1.0 : 0.0;
	}
	bool TimeGate(tw::Float t,tw::Float *add) {
		bool gateOpen = false;
		switch (timingMethod)
		{
			case tw::profile::timing::triggered:
				*add = 1.0;
				gateOpen = t>=t0 && !wasTriggered;
				break;
			case tw::profile::timing::maintained:
				*add = 0.0;
				gateOpen = t>=t0 && t<=t1;
				break;
		}
		wasTriggered = wasTriggered || gateOpen;
		return gateOpen;
	}
	/// interpolate a value from a 1D non-uniform mesh
	tw::Float Interpolate(tw::Float s0,std::vector<tw::Float>& s,std::vector<tw::Float>& f) {
		auto cmp = [](tw::Float x,tw::Float y) {
			return x < y;
		};
		// first element satisfying s[i] >= s0
		auto it = std::lower_bound(s.begin(),s.end(),s0,cmp);
		auto i = it - s.begin();
		if (i<=0 || i >= f.size()) {
			return 0.0;
		}
		i--;
		auto w = (s[i+1]-s0)/(s[i+1]-s[i]);
		if (segmentShape==tw::profile::shape::quartic)
			w = QuarticPulse(0.5*w);
		if (segmentShape==tw::profile::shape::quintic)
			w = QuinticRise(w);
		return f[i]*w + f[i+1]*(1.0-w);
	};
	virtual bool ReadInputFileDirective(const TSTreeCursor *curs0,const std::string& src) {
		// assume whoever calls this is passing over comments and anonymous
		auto curs = tw::input::Cursor(curs0);
		std::string outer = tw::input::node_kind(curs.get());
		if (outer != "assignment")
			throw tw::FatalError("Expected assignment, got " + outer + ", at " + tw::input::loc_str(curs0));
		ts_tree_cursor_goto_first_child(curs.get());
		std::string com = tw::input::node_text(curs.get(),src);
		if (com=="particle weight") {
			std::string weighting;
			tw::input::String directive(&weighting,true);
			directive.Read(curs.get(),src,"particle weight",native);
			if (weighting!="variable" && weighting!="fixed")
				throw tw::FatalError("Invalid type <"+weighting+"> while processing key <particle weight>.");
			variableCharge = (weighting=="variable" ? true : false);
			return true;
		}
		if (com=="euler angles") // eg, euler angles = ( 45[deg] 90[deg] 30[deg] )
		{
			tw::vec3 angles;
			tw::input::Vec3 directive(&angles);
			directive.Read(curs.get(),src,"euler angles",native);
			orientation.SetWithEulerAngles( angles.x , angles.y , angles.z );
			return true;
		}
		return false;
	}
	virtual void ReadCheckpoint(std::ifstream& inFile)	{
		ComputeTool::ReadCheckpoint(inFile);
		inFile.read((char *)&centerPt,sizeof(tw::vec3));
		inFile.read((char *)&orientation,sizeof(orientation));
		inFile.read((char *)&wasTriggered,sizeof(bool));
		inFile.read((char *)&t0,sizeof(tw::Float));
		inFile.read((char *)&t1,sizeof(tw::Float));
	}

	virtual void WriteCheckpoint(std::ofstream& outFile) {
		ComputeTool::WriteCheckpoint(outFile);
		outFile.write((char *)&centerPt,sizeof(tw::vec3));
		outFile.write((char *)&orientation,sizeof(orientation));
		outFile.write((char *)&wasTriggered,sizeof(bool));
		outFile.write((char *)&t0,sizeof(tw::Float));
		outFile.write((char *)&t1,sizeof(tw::Float));
	}
};

export struct UniformProfile:Profile
{
	tw::Float density;

	UniformProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
};

export struct GaussianProfile:Profile
{
	tw::Float density;
	tw::vec3 beamSize;

	GaussianProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
};

export struct ChannelProfile:Profile
{
	std::vector<tw::Float> z,fz;
	tw::Float coeff[4];

	ChannelProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

export struct ColumnProfile:Profile
{
	std::vector<tw::Float> z,fz;
	tw::vec3 beamSize;

	ColumnProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

export struct PiecewiseProfile:Profile
{
	std::vector<tw::Float> x,fx,y,fy,z,fz;

	PiecewiseProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

export struct CorrugatedProfile:Profile
{
	tw::Float a0,gamma0,w0,wp0,km,delta,rchannel;

	CorrugatedProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
};

UniformProfile::UniformProfile(const std::string& name,MetricSpace *m,Task *tsk):Profile(name,m,tsk)
{
	directives.Add("density",new tw::input::Float(&density));
}

tw::Float UniformProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	return theRgn->Inside(Boost(pos),ds) ? gammaBoost*density : 0.0;
}

GaussianProfile::GaussianProfile(const std::string& name,MetricSpace *m,Task *tsk):Profile(name,m,tsk)
{
	directives.Add("density",new tw::input::Float(&density));
	directives.Add("size",new tw::input::Vec3(&beamSize));
}

tw::Float GaussianProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	tw::Float dens = density;
	tw::vec3 b = Boost(pos);
	tw::vec3 p = Translate_Rotate(b);
	dens *= std::exp(-sqr(p.x/beamSize.x));
	dens *= std::exp(-sqr(p.y/beamSize.y));
	dens *= std::exp(-sqr(p.z/beamSize.z));
	return theRgn->Inside(b,ds) ? gammaBoost*dens : 0.0;
}

ChannelProfile::ChannelProfile(const std::string& name,MetricSpace *m,Task *tsk):Profile(name,m,tsk)
{
	directives.Add("coefficients",new tw::input::Numbers<tw::Float>(&coeff[0],4));
	directives.Add("zpoints",new tw::input::NumberList<std::vector<tw::Float>>(&z));
	directives.Add("zdensity",new tw::input::NumberList<std::vector<tw::Float>>(&fz));
}

void ChannelProfile::ReadCheckpoint(std::ifstream& inFile)
{
	Profile::ReadCheckpoint(inFile);
	inFile.read((char *)&z[0],sizeof(tw::Float)*z.size());
	inFile.read((char *)&fz[0],sizeof(tw::Float)*fz.size());
}

void ChannelProfile::WriteCheckpoint(std::ofstream& outFile)
{
	Profile::WriteCheckpoint(outFile);
	outFile.write((char *)&z[0],sizeof(tw::Float)*z.size());
	outFile.write((char *)&fz[0],sizeof(tw::Float)*fz.size());
}

tw::Float ChannelProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	tw::Int i;
	tw::Float r2,w,dens = 0.0;
	tw::vec3 b = Boost(pos);
	tw::vec3 p = Translate_Rotate(b);
	dens = Interpolate(p.z,z,fz);
	r2 = sqr(p.x) + sqr(p.y);
	dens *= coeff[0] + coeff[1]*r2 + coeff[2]*r2*r2 + coeff[3]*r2*r2*r2;
	return theRgn->Inside(b,ds) ? gammaBoost*dens : 0.0;
}

ColumnProfile::ColumnProfile(const std::string& name,MetricSpace *m,Task *tsk):Profile(name,m,tsk)
{
	directives.Add("size",new tw::input::Vec3(&beamSize));
	directives.Add("zpoints",new tw::input::NumberList<std::vector<tw::Float>>(&z));
	directives.Add("zdensity",new tw::input::NumberList<std::vector<tw::Float>>(&fz));
}

void ColumnProfile::ReadCheckpoint(std::ifstream& inFile)
{
	Profile::ReadCheckpoint(inFile);
	inFile.read((char *)&z[0],sizeof(tw::Float)*z.size());
	inFile.read((char *)&fz[0],sizeof(tw::Float)*fz.size());
}

void ColumnProfile::WriteCheckpoint(std::ofstream& outFile)
{
	Profile::WriteCheckpoint(outFile);
	outFile.write((char *)&z[0],sizeof(tw::Float)*z.size());
	outFile.write((char *)&fz[0],sizeof(tw::Float)*fz.size());
}

tw::Float ColumnProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	tw::Int i;
	tw::Float w,dens = 0.0;
	tw::vec3 b = Boost(pos);
	tw::vec3 p = Translate_Rotate(b);
	dens = Interpolate(p.z,z,fz);
	dens *= std::exp(-sqr(p.x/beamSize.x));
	dens *= std::exp(-sqr(p.y/beamSize.y));
	return theRgn->Inside(b,ds) ? gammaBoost*dens : 0.0;
}

PiecewiseProfile::PiecewiseProfile(const std::string& name,MetricSpace *m,Task *tsk) : Profile(name,m,tsk)
{
	directives.Add("xpoints",new tw::input::NumberList<std::vector<tw::Float>>(&x),false);
	directives.Add("ypoints",new tw::input::NumberList<std::vector<tw::Float>>(&y),false);
	directives.Add("zpoints",new tw::input::NumberList<std::vector<tw::Float>>(&z),false);
	directives.Add("xdensity",new tw::input::NumberList<std::vector<tw::Float>>(&fx),false);
	directives.Add("ydensity",new tw::input::NumberList<std::vector<tw::Float>>(&fy),false);
	directives.Add("zdensity",new tw::input::NumberList<std::vector<tw::Float>>(&fz),false);
}

void PiecewiseProfile::Initialize()
{
	Profile::Initialize();

	if (x.size()<2)
	{
		x.resize(2);
		fx.resize(2);
		theRgn->GetBoxLim(&x[0],&x[1],1);
		fx.assign({1,1});
	}
	if (y.size()<2)
	{
		y.resize(2);
		fy.resize(2);
		theRgn->GetBoxLim(&y[0],&y[1],2);
		fy.assign({1,1});
	}
	if (z.size()<2)
	{
		z.resize(2);
		fz.resize(2);
		theRgn->GetBoxLim(&z[0],&z[1],3);
		fz.assign({1,1});
	}
}

void PiecewiseProfile::ReadCheckpoint(std::ifstream& inFile)
{
	Profile::ReadCheckpoint(inFile);
	inFile.read((char *)&x[0],sizeof(tw::Float)*x.size());
	inFile.read((char *)&fx[0],sizeof(tw::Float)*fx.size());
	inFile.read((char *)&y[0],sizeof(tw::Float)*y.size());
	inFile.read((char *)&fy[0],sizeof(tw::Float)*fy.size());
	inFile.read((char *)&z[0],sizeof(tw::Float)*z.size());
	inFile.read((char *)&fz[0],sizeof(tw::Float)*fz.size());
}

void PiecewiseProfile::WriteCheckpoint(std::ofstream& outFile)
{
	Profile::WriteCheckpoint(outFile);
	outFile.write((char *)&x[0],sizeof(tw::Float)*x.size());
	outFile.write((char *)&fx[0],sizeof(tw::Float)*fx.size());
	outFile.write((char *)&y[0],sizeof(tw::Float)*y.size());
	outFile.write((char *)&fy[0],sizeof(tw::Float)*fy.size());
	outFile.write((char *)&z[0],sizeof(tw::Float)*z.size());
	outFile.write((char *)&fz[0],sizeof(tw::Float)*fz.size());
}

tw::Float PiecewiseProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	tw::Float ansX,ansY,ansZ;
	tw::Int i;
	tw::Float r,w;

	tw::vec3 b = Boost(pos);
	tw::vec3 p = Translate_Rotate(b);

	const tw::Float x0 = p.x;
	const tw::Float y0 = p.y;
	const tw::Float z0 = p.z;

	tw::Int xDim = x.size();
	tw::Int yDim = y.size();
	tw::Int zDim = z.size();

	ansX = ansY = ansZ = 1.0;

	switch (symmetry)
	{
		case tw::grid::cartesian:
			ansX = Interpolate(x0,x,fx);
			ansY = Interpolate(y0,y,fy);
			ansZ = Interpolate(z0,z,fz);
			break;
		case tw::grid::cylindrical:
			r = std::sqrt(sqr(x0 - x[0]) + sqr(y0 - 0.5*(y[0] + y[yDim-1])));
			r += x[0];
			ansX = Interpolate(r,x,fx);
			ansY = 1.0;
			ansZ = Interpolate(z0,z,fz);
			break;
		case tw::grid::spherical:
			r = std::sqrt(sqr(x0 - x[0]) + sqr(y0 - 0.5*(y[0] + y[yDim-1])) + sqr(z0 - 0.5*(z[0] + z[zDim-1])));
			r += x[0];
			ansX = Interpolate(r,x,fx);
			ansY = 1.0;
			ansZ = 1.0;
			break;
	}

	tw::Float dens = ansX*ansY*ansZ*sqr(std::cos(0.5*modeNumber.x*p.x)*std::cos(0.5*modeNumber.y*p.y)*std::cos(0.5*modeNumber.z*p.z));
	return theRgn->Inside(b,ds) ? gammaBoost*dens : 0.0;
}

CorrugatedProfile::CorrugatedProfile(const std::string& name,MetricSpace *m,Task *tsk) : Profile(name,m,tsk)
{
	directives.Add("a0",new tw::input::Float(&a0));
	directives.Add("gamma0",new tw::input::Float(&gamma0));
	directives.Add("w0",new tw::input::Float(&w0));
	directives.Add("wp0",new tw::input::Float(&wp0));
	directives.Add("km",new tw::input::Float(&km));
	directives.Add("delta",new tw::input::Float(&delta));
	directives.Add("rchannel",new tw::input::Float(&rchannel));
}

tw::Float CorrugatedProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	tw::Float dens,wp1s,r2,psi,a0Hat,kHat;
	tw::vec3 b = Boost(pos);
	tw::vec3 p = Translate_Rotate(b);
	const tw::Float x = p.x;
	const tw::Float y = p.y;
	const tw::Float z = p.z; // z = 0 is initial injection point
	r2 = x*x + y*y;
	psi = delta*wp0*wp0/(2.0*w0*km);
	a0Hat = 4.0*tw::cyl_bessel_j(1,psi)*a0/(w0*rchannel);
	kHat = w0 + km - 0.5*w0*(sqr(wp0/w0) + 8.0/sqr(w0*rchannel));
	wp1s = 2.0*w0*kHat - 2.0*std::sqrt(sqr(gamma0+a0Hat*w0*z)/(sqr(gamma0+a0Hat*w0*z)-1.0))*w0*w0;
	dens = wp0*wp0*(1.0 + delta*std::sin(km*z)) + wp1s + 4.0*r2/std::pow(rchannel,tw::Float(4.0));
	return theRgn->Inside(b,ds) ? gammaBoost*dens : 0.0;
}

