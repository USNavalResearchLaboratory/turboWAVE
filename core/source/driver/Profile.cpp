module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"
#ifndef USE_STD_MODULE
#include <algorithm>
#endif

export module driver:profile;
import :engine;
import :region;

import input;
import functions;
import logger;

using namespace tw::bc;

export namespace tw
{
	namespace profile
	{
		enum class quantity { density,energy,power,px,py,pz };
		enum class timing { triggered,maintained };
		enum class loading { statistical,deterministic };
		enum class shape { triangle,sin2,quartic,quintic,sech };
	}
}

export struct Profile : Engine
{
	tw::profile::shape segmentShape;
	tw::grid::geometry symmetry;
	tw::vec4 modeNumber;
	tw::Float modeAmplitude;
	/// This is the user's request to boost the profile's parameters
	tw::Float gammaBoost;

	/// displace the body before the rotation
	tw::vec3 origin;
	/// rotate the body
	tw::basis orientation;
	/// displace the body after the rotation
	tw::vec3 translation;
	/// in the active view this gets applied as new_vec = Rz(alpha)*Rx(beta)*Rz(gamma)*old_vec,
	/// i.e. the "last angle" gets applied first.
	tw::vec3 euler;

	// items needed for particle/fluid loading
private:
	tw::vec3 driftMomentum;
protected:
	std::shared_ptr<Region> theRgn;
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
		t0 = tw::big_neg;
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
			{ "power",tw::profile::quantity::power },
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
		directives.Add("origin",new tw::input::Vec3(&origin),false);
		directives.Add("euler angles",new tw::input::Vec3(&euler),false);
		directives.Add("translation",new tw::input::Vec3(&translation),false);
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
		directives.Add("mode number",new tw::input::Vec4(&modeNumber),false);
		directives.Add("boosted frame gamma",new tw::input::Float(&gammaBoost),false);
		directives.Add("particle weight",new tw::input::Custom,false);
	}
	virtual void Initialize() {
		Engine::Initialize();
		orientation.SetWithEulerAngles(euler.x, euler.y, euler.z);
		theRgn = std::dynamic_pointer_cast<Region>(region);
		if (!theRgn) {
			logger::DEBUG(std::format("<{}> using default entire region",name));
			auto temp = std::make_shared<SimpleRegion>("default_entire",space,task,
				std::make_unique<EntireRegion>("entire",space,task));
			theRgn = temp;
		}
	}
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
	/// Boost, translate, and rotate from the simulation frame to the profile's frame
	/// The user specifies a boosted frame for each profile.
	/// If gammaBoost=1 then the profile is in the simulation's frame, which may be considered boosted, or not.
	/// If gammaBoost>1 then the profile is in a different frame, perhaps a lab frame, in this case the
	/// profile will be boosted into the simulation's frame.
	/// The argument is in the simulation's frame, the return value is in the profile's frame.
	tw::vec4 TransformPoint(const tw::vec4& pos) const {
		auto ans(pos);
		ans.zBoost(gammaBoost,1.0);
		ans.Sub3(translation);
		orientation.ExpressInBasis(&ans);
		ans.Sub3(origin);
		return ans;
	}
	virtual tw::Float GetValue(const tw::vec4& pos,const MetricSpace& ds) {
		return theRgn->Inside(TransformPoint(pos),0) ? 1.0 : 0.0;
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
				*add = whichQuantity == tw::profile::quantity::power ? 1.0: 0.0;
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
		return false;
	}
	virtual void ReadCheckpoint(std::ifstream& inFile)	{
		ComputeTool::ReadCheckpoint(inFile);
		inFile.read((char *)&origin,sizeof(tw::vec3));
		inFile.read((char *)&orientation,sizeof(orientation));
		inFile.read((char *)&translation,sizeof(tw::vec3));
		inFile.read((char *)&wasTriggered,sizeof(bool));
		inFile.read((char *)&t0,sizeof(tw::Float));
		inFile.read((char *)&t1,sizeof(tw::Float));
	}

	virtual void WriteCheckpoint(std::ofstream& outFile) {
		ComputeTool::WriteCheckpoint(outFile);
		outFile.write((char *)&origin,sizeof(tw::vec3));
		outFile.write((char *)&orientation,sizeof(orientation));
		outFile.write((char *)&translation,sizeof(tw::vec3));
		outFile.write((char *)&wasTriggered,sizeof(bool));
		outFile.write((char *)&t0,sizeof(tw::Float));
		outFile.write((char *)&t1,sizeof(tw::Float));
	}
};

export struct UniformProfile:Profile
{
	tw::Float density;

	UniformProfile(const std::string& name,MetricSpace *m,Task *tsk):Profile(name,m,tsk) {
		density = 0.0; // hydro is allowed to create this profile automatically
		directives.Add("density",new tw::input::Float(&density));
	}
	virtual tw::Float GetValue(const tw::vec4& pos,const MetricSpace& ds) {
		return theRgn->Inside(TransformPoint(pos),0) ? gammaBoost*density : 0.0;
	}
};

export struct GaussianProfile:Profile
{
	tw::Float density;
	tw::vec4 beamSize;
	GaussianProfile(const std::string& name,MetricSpace *m,Task *tsk):Profile(name,m,tsk) {
		directives.Add("density",new tw::input::Float(&density));
		directives.Add("size",new tw::input::Vec4(&beamSize));
	}
	virtual tw::Float GetValue(const tw::vec4& pos,const MetricSpace& ds) {
		tw::Float dens = gammaBoost*density;
		auto p = TransformPoint(pos);
		for (auto i=0; i<4; i++) {
			dens *= std::exp(-sqr(p[i]/beamSize[i]));
		}
		return theRgn->Inside(p,0) ? dens : 0.0;
	}
};

export struct ChannelProfile:Profile
{
	std::vector<tw::Float> z,fz;
	tw::Float coeff[4];
	ChannelProfile(const std::string& name,MetricSpace *m,Task *tsk):Profile(name,m,tsk) {
		directives.Add("coefficients",new tw::input::Numbers<tw::Float>(&coeff[0],4));
		directives.Add("zpoints",new tw::input::NumberList<std::vector<tw::Float>>(&z));
		directives.Add("zdensity",new tw::input::NumberList<std::vector<tw::Float>>(&fz));
	}

	virtual tw::Float GetValue(const tw::vec4& pos,const MetricSpace& ds) {
		tw::Int i;
		tw::Float r2,w,dens = 0.0;
		auto p = TransformPoint(pos);
		dens = Interpolate(p[3],z,fz);
		r2 = sqr(p[1]) + sqr(p[2]);
		dens *= coeff[0] + coeff[1]*r2 + coeff[2]*r2*r2 + coeff[3]*r2*r2*r2;
		return theRgn->Inside(p,0) ? gammaBoost*dens : 0.0;
	}
};

export struct ColumnProfile:Profile
{
	std::vector<tw::Float> z,fz;
	tw::vec4 beamSize;
	ColumnProfile(const std::string& name,MetricSpace *m,Task *tsk):Profile(name,m,tsk) {
		directives.Add("size",new tw::input::Vec4(&beamSize));
		directives.Add("zpoints",new tw::input::NumberList<std::vector<tw::Float>>(&z));
		directives.Add("zdensity",new tw::input::NumberList<std::vector<tw::Float>>(&fz));
	}
	virtual tw::Float GetValue(const tw::vec4& pos,const MetricSpace& ds) {
		tw::Int i;
		tw::Float w,dens = 0.0;
		auto p = TransformPoint(pos);
		dens = Interpolate(p[3],z,fz);
		dens *= std::exp(-sqr(p[0]/beamSize[0]));
		dens *= std::exp(-sqr(p[1]/beamSize[1]));
		dens *= std::exp(-sqr(p[2]/beamSize[2]));
		return theRgn->Inside(p,0) ? gammaBoost*dens : 0.0;
	}
};

export struct PiecewiseProfile:Profile
{
	std::vector<tw::Float> t,ft,x,fx,y,fy,z,fz;
	PiecewiseProfile(const std::string& name,MetricSpace *m,Task *tsk) : Profile(name,m,tsk) {
		directives.Add("tpoints",new tw::input::NumberList<std::vector<tw::Float>>(&t),false);
		directives.Add("xpoints",new tw::input::NumberList<std::vector<tw::Float>>(&x),false);
		directives.Add("ypoints",new tw::input::NumberList<std::vector<tw::Float>>(&y),false);
		directives.Add("zpoints",new tw::input::NumberList<std::vector<tw::Float>>(&z),false);
		directives.Add("tdensity",new tw::input::NumberList<std::vector<tw::Float>>(&ft),false);
		directives.Add("xdensity",new tw::input::NumberList<std::vector<tw::Float>>(&fx),false);
		directives.Add("ydensity",new tw::input::NumberList<std::vector<tw::Float>>(&fy),false);
		directives.Add("zdensity",new tw::input::NumberList<std::vector<tw::Float>>(&fz),false);
	}
	void Initialize() {
		Profile::Initialize();
		auto bounds = theRgn->Bounds(0);

		if (t.size()<2)
		{
			t.resize(2);
			ft.resize(2);
			t.assign({tw::big_neg,tw::big_pos});
			ft.assign({1,1});
		}
		if (x.size()<2)
		{
			x.resize(2);
			fx.resize(2);
			x.assign({bounds[0],bounds[1]});
			fx.assign({1,1});
		}
		if (y.size()<2)
		{
			y.resize(2);
			fy.resize(2);
			y.assign({bounds[2],bounds[3]});
			fy.assign({1,1});
		}
		if (z.size()<2)
		{
			z.resize(2);
			fz.resize(2);
			z.assign({bounds[4],bounds[5]});
			fz.assign({1,1});
		}
	}
	virtual tw::Float GetValue(const tw::vec4& pos,const MetricSpace& ds) {
		tw::Float r;
		tw::vec4 ans(1,1,1,1);

		auto p = TransformPoint(pos);

		const tw::Int yDim = y.size();
		const tw::Int zDim = z.size();

		ans[0] = Interpolate(p[0],t,ft);
		switch (symmetry)
		{
			case tw::grid::cartesian:
				ans[1] = Interpolate(p[1],x,fx);
				ans[2] = Interpolate(p[2],y,fy);
				ans[3] = Interpolate(p[3],z,fz);
				break;
			case tw::grid::cylindrical:
				r = std::sqrt(sqr(p[1] - x[0]) + sqr(p[2] - 0.5*(y[0] + y[yDim-1])));
				r += x[0];
				ans[1] = Interpolate(r,x,fx);
				ans[3] = Interpolate(p[3],z,fz);
				break;
			case tw::grid::spherical:
				r = std::sqrt(sqr(p[1] - x[0]) + sqr(p[2] - 0.5*(y[0] + y[yDim-1])) + sqr(p[3] - 0.5*(z[0] + z[zDim-1])));
				r += x[0];
				ans[1] = Interpolate(r,x,fx);
				break;
		}

		tw::Float dens = gammaBoost * ans[0]*ans[1]*ans[2]*ans[3];
		for (auto i=0; i<4 ; i++) {
			dens *= sqr(std::cos(0.5*modeNumber[i]*p[i]));
		}
		return theRgn->Inside(p,0) ? dens : 0.0;
	}
};

export struct CorrugatedProfile:Profile
{
	tw::Float a0,gamma0,w0,wp0,km,delta,rchannel;
	CorrugatedProfile(const std::string& name,MetricSpace *m,Task *tsk) : Profile(name,m,tsk) {
		directives.Add("a0",new tw::input::Float(&a0));
		directives.Add("gamma0",new tw::input::Float(&gamma0));
		directives.Add("w0",new tw::input::Float(&w0));
		directives.Add("wp0",new tw::input::Float(&wp0));
		directives.Add("km",new tw::input::Float(&km));
		directives.Add("delta",new tw::input::Float(&delta));
		directives.Add("rchannel",new tw::input::Float(&rchannel));
	}
	virtual tw::Float GetValue(const tw::vec4& pos,const MetricSpace& ds) {
		tw::Float dens,wp1s,r2,psi,a0Hat,kHat;
		auto p = TransformPoint(pos);
		const tw::Float x = p[1];
		const tw::Float y = p[2];
		const tw::Float z = p[3]; // z = 0 is initial injection point
		r2 = x*x + y*y;
		psi = delta*wp0*wp0/(2.0*w0*km);
		a0Hat = 4.0*tw::cyl_bessel_j(1,psi)*a0/(w0*rchannel);
		kHat = w0 + km - 0.5*w0*(sqr(wp0/w0) + 8.0/sqr(w0*rchannel));
		wp1s = 2.0*w0*kHat - 2.0*std::sqrt(sqr(gamma0+a0Hat*w0*z)/(sqr(gamma0+a0Hat*w0*z)-1.0))*w0*w0;
		dens = wp0*wp0*(1.0 + delta*std::sin(km*z)) + wp1s + 4.0*r2/std::pow(rchannel,tw::Float(4.0));
		return theRgn->Inside(p,0) ? gammaBoost*dens : 0.0;
	}
};





