#include "meta_base.h"
#include "computeTool.h"
#include "injection.h"
using namespace tw::bc;

///////////////////
//               //
// PROFILE CLASS //
//               //
///////////////////

Profile::Profile(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
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

void Profile::Initialize()
{
}

void Profile::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	std::string word;
	if (com=="particle weight")
	{
		inputString >> word >> word;
		if (word!="variable" && word!="fixed")
			throw tw::FatalError("Invalid type <"+word+"> while processing key <particle weight>.");
		variableCharge = (word=="variable" ? true : false);
	}
	if (com=="euler angles") // eg, euler angles = ( %45deg %90deg %30deg )
	{
		tw::Float alpha,beta,gamma;
		inputString >> alpha >> beta >> gamma;
		orientation.SetWithEulerAngles(alpha,beta,gamma);
	}
}

void Profile::ReadCheckpoint(std::ifstream& inFile)
{
	ComputeTool::ReadCheckpoint(inFile);
	inFile.read((char *)&centerPt,sizeof(tw::vec3));
	inFile.read((char *)&orientation,sizeof(orientation));
	inFile.read((char *)&wasTriggered,sizeof(bool));
	inFile.read((char *)&t0,sizeof(tw::Float));
	inFile.read((char *)&t1,sizeof(tw::Float));
}

void Profile::WriteCheckpoint(std::ofstream& outFile)
{
	ComputeTool::WriteCheckpoint(outFile);
	outFile.write((char *)&centerPt,sizeof(tw::vec3));
	outFile.write((char *)&orientation,sizeof(orientation));
	outFile.write((char *)&wasTriggered,sizeof(bool));
	outFile.write((char *)&t0,sizeof(tw::Float));
	outFile.write((char *)&t1,sizeof(tw::Float));
}

tw::vec3 Profile::DriftMomentum(const tw::Float& mass)
{
	tw::vec4 v4(0.0,driftMomentum/mass);
	tw::Float gb2 = v4 ^ v4;
	v4[0] = sqrt(1.0 + gb2);
	v4.zBoost(gammaBoost,-1.0);
	return mass*v4.spatial();
}

tw::Float Profile::Temperature(const tw::Float& mass)
{
	if (temperature!=0.0)
		return temperature;
	else
	 	return sqr(thermalMomentum.x)/mass; // appropriate for exp(-v^2/(2*vth^2)) convention
}

tw::vec3 Profile::Boost(const tw::vec3& pos)
{
	// Here the boost only works if the profile is constant in time.
	// The function's caller is giving us boosted frame coordinates.
	// The user is giving us lab frame coordinates.
	// Therefore first transform arguments to lab frame, then proceed as usual.
	tw::vec4 x4(0.0,pos);
	x4.zBoost(gammaBoost,1.0);
	return x4.spatial();
}

tw::vec3 Profile::Translate_Rotate(const tw::vec3& pos)
{
	// The argument should already be boosted.
	// N.b. the boost applies to both region and profile, while translate-rotate applies only to the profile.
	tw::vec3 p = pos - centerPt;
	orientation.ExpressInBasis(&p);
	return p;
}

tw::Float Profile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	return theRgn->Inside(Boost(pos),ds) ? 1.0 : 0.0;
}

bool Profile::TimeGate(tw::Float t,tw::Float *add)
{
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
	dens *= exp(-sqr(p.x/beamSize.x));
	dens *= exp(-sqr(p.y/beamSize.y));
	dens *= exp(-sqr(p.z/beamSize.z));
	return theRgn->Inside(b,ds) ? gammaBoost*dens : 0.0;
}

ChannelProfile::ChannelProfile(const std::string& name,MetricSpace *m,Task *tsk):Profile(name,m,tsk)
{
	directives.Add("coefficients",new tw::input::Numbers<tw::Float>(&coeff[0],4));
	directives.Add("zpoints",new tw::input::List<std::valarray<tw::Float>>(&z));
	directives.Add("zdensity",new tw::input::List<std::valarray<tw::Float>>(&fz));
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
	for (i=0;i<z.size()-1;i++)
	{
		if (p.z>=z[i] && p.z<=z[i+1])
		{
			w = (z[i+1]-p.z)/(z[i+1]-z[i]);
			if (segmentShape==tw::profile::shape::quartic)
				w = QuarticPulse(0.5*w);
			if (segmentShape==tw::profile::shape::quintic)
				w = QuinticRise(w);
			dens = fz[i]*w + fz[i+1]*(1.0-w);
		}
	}
	r2 = sqr(p.x) + sqr(p.y);
	dens *= coeff[0] + coeff[1]*r2 + coeff[2]*r2*r2 + coeff[3]*r2*r2*r2;
	return theRgn->Inside(b,ds) ? gammaBoost*dens : 0.0;
}

ColumnProfile::ColumnProfile(const std::string& name,MetricSpace *m,Task *tsk):Profile(name,m,tsk)
{
	directives.Add("size",new tw::input::Vec3(&beamSize));
	directives.Add("zpoints",new tw::input::List<std::valarray<tw::Float>>(&z));
	directives.Add("zdensity",new tw::input::List<std::valarray<tw::Float>>(&fz));
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
	for (i=0;i<z.size()-1;i++)
	{
		if (p.z>=z[i] && p.z<=z[i+1])
		{
			w = (z[i+1]-p.z)/(z[i+1]-z[i]);
			if (segmentShape==tw::profile::shape::quartic)
				w = QuarticPulse(0.5*w);
			if (segmentShape==tw::profile::shape::quintic)
				w = QuinticRise(w);
			dens = fz[i]*w + fz[i+1]*(1.0-w);
		}
	}
	dens *= exp(-sqr(p.x/beamSize.x));
	dens *= exp(-sqr(p.y/beamSize.y));
	return theRgn->Inside(b,ds) ? gammaBoost*dens : 0.0;
}

PiecewiseProfile::PiecewiseProfile(const std::string& name,MetricSpace *m,Task *tsk) : Profile(name,m,tsk)
{
	directives.Add("xpoints",new tw::input::List<std::valarray<tw::Float>>(&x),false);
	directives.Add("ypoints",new tw::input::List<std::valarray<tw::Float>>(&y),false);
	directives.Add("zpoints",new tw::input::List<std::valarray<tw::Float>>(&z),false);
	directives.Add("xdensity",new tw::input::List<std::valarray<tw::Float>>(&fx),false);
	directives.Add("ydensity",new tw::input::List<std::valarray<tw::Float>>(&fy),false);
	directives.Add("zdensity",new tw::input::List<std::valarray<tw::Float>>(&fz),false);
}

void PiecewiseProfile::Initialize()
{
	Profile::Initialize();
	if (x.size()<2)
	{
		x.resize(2);
		fx.resize(2);
		theRgn->GetBoxLim(&x[0],&x[1],1);
		fx = 1.0;
	}
	if (y.size()<2)
	{
		y.resize(2);
		fy.resize(2);
		theRgn->GetBoxLim(&y[0],&y[1],2);
		fy = 1.0;
	}
	if (z.size()<2)
	{
		z.resize(2);
		fz.resize(2);
		theRgn->GetBoxLim(&z[0],&z[1],3);
		fz = 1.0;
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
			ansX = 0.0;
			for (i=0;i<xDim-1;i++)
			{
				if (x0>=x[i] && x0<=x[i+1])
				{
					w = (x[i+1]-x0)/(x[i+1]-x[i]);
					if (segmentShape==tw::profile::shape::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==tw::profile::shape::quintic)
						w = QuinticRise(w);
					ansX = fx[i]*w + fx[i+1]*(1.0-w);
				}
			}

			ansY = 0.0;
			for (i=0;i<yDim-1;i++)
			{
				if (y0>=y[i] && y0<=y[i+1])
				{
					w = (y[i+1]-y0)/(y[i+1]-y[i]);
					if (segmentShape==tw::profile::shape::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==tw::profile::shape::quintic)
						w = QuinticRise(w);
					ansY = fy[i]*w + fy[i+1]*(1.0-w);
				}
			}

			ansZ = 0.0;
			for (i=0;i<zDim-1;i++)
			{
				if (z0>=z[i] && z0<=z[i+1])
				{
					w = (z[i+1]-z0)/(z[i+1]-z[i]);
					if (segmentShape==tw::profile::shape::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==tw::profile::shape::quintic)
						w = QuinticRise(w);
					ansZ = fz[i]*w + fz[i+1]*(1.0-w);
				}
			}
			break;
		case tw::grid::cylindrical:
			r = sqrt(sqr(x0 - x[0]) + sqr(y0 - 0.5*(y[0] + y[yDim-1])));
			r += x[0];
			ansX = 0.0;
			for (i=0;i<xDim-1;i++)
			{
				if (r>=x[i] && r<=x[i+1])
				{
					w = (x[i+1]-r)/(x[i+1]-x[i]);
					if (segmentShape==tw::profile::shape::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==tw::profile::shape::quintic)
						w = QuinticRise(w);
					ansX = fx[i]*w + fx[i+1]*(1.0-w);
				}
			}

			ansZ = 0.0;
			for (i=0;i<zDim-1;i++)
			{
				if (z0>=z[i] && z0<=z[i+1])
				{
					w = (z[i+1]-z0)/(z[i+1]-z[i]);
					if (segmentShape==tw::profile::shape::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==tw::profile::shape::quintic)
						w = QuinticRise(w);
					ansZ = fz[i]*w + fz[i+1]*(1.0-w);
				}
			}

			ansY = 1.0;
			break;
		case tw::grid::spherical:
			r = sqrt(sqr(x0 - x[0]) + sqr(y0 - 0.5*(y[0] + y[yDim-1])) + sqr(z0 - 0.5*(z[0] + z[zDim-1])));
			r += x[0];
			ansX = 0.0;
			for (i=0;i<xDim-1;i++)
			{
				if (r>=x[i] && r<=x[i+1])
				{
					w = (x[i+1]-r)/(x[i+1]-x[i]);
					if (segmentShape==tw::profile::shape::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==tw::profile::shape::quintic)
						w = QuinticRise(w);
					ansX = fx[i]*w + fx[i+1]*(1.0-w);
				}
			}

			ansY = 1.0;
			ansZ = 1.0;
			break;
	}

	tw::Float dens = ansX*ansY*ansZ*sqr(cos(0.5*modeNumber.x*p.x)*cos(0.5*modeNumber.y*p.y)*cos(0.5*modeNumber.z*p.z));
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
	wp1s = 2.0*w0*kHat - 2.0*sqrt(sqr(gamma0+a0Hat*w0*z)/(sqr(gamma0+a0Hat*w0*z)-1.0))*w0*w0;
	dens = wp0*wp0*(1.0 + delta*sin(km*z)) + wp1s + 4.0*r2/pow(rchannel,tw::Float(4.0));
	return theRgn->Inside(b,ds) ? gammaBoost*dens : 0.0;
}



///////////////////////////
//                       //
// PULSE and BEAM SHAPES //
//                       //
///////////////////////////

PulseShape::PulseShape()
{
	whichProfile = tw::profile::shape::quintic;
	delay = 0.0;
	risetime = 1.0;
	holdtime = 0.0;
	falltime = 1.0;
	t1 = 0.0;
	t2 = 1.0;
	t3 = 0.0;
	t4 = 2.0;
}

void PulseShape::Initialize(const tw::Float time_origin)
{
	t1 = delay - time_origin;
	t2 = delay + risetime - time_origin;
	t3 = delay + risetime + holdtime - time_origin;
	t4 = delay + risetime + holdtime + falltime - time_origin;
}

tw::Float PulseShape::PulseShapeFactor(const tw::Float t) const
{
	const tw::Float hold = tw::Float(t > t2 && t <= t3);
	const tw::Float tau_rise = tw::Float(t > t1 && t <= t2) * (t-t1) / (t2-t1);
	const tw::Float tau_fall = tw::Float(t > t3 && t <= t4) * (1.0 - (t-t3) / (t4-t3));
	const tw::Float tau = tau_rise + hold + tau_fall;

	switch (whichProfile)
	{
		case tw::profile::shape::triangle:
			return tau;
		case tw::profile::shape::sech:
			return SechRise(tau);
		case tw::profile::shape::quartic:
			return QuarticRise(tau);
		case tw::profile::shape::quintic:
			return QuinticRise(tau);
		case tw::profile::shape::sin2:
			return Sin2Rise(tau);
	}
	return 0.0;
}

tw::Float PulseShape::D1Amplitude(const tw::Float t) const
{
	const tw::Float hold = tw::Float(t > t2 && t <= t3);
	const tw::Float tau_rise = tw::Float(t > t1 && t <= t2) * (t-t1) / (t2-t1);
	const tw::Float tau_fall = tw::Float(t > t3 && t <= t4) * (1.0 - (t-t3) / (t4-t3));
	const tw::Float tau = tau_rise + hold + tau_fall;
	tw::Float w = tw::Float(t > t1 && t <= t2) / (t2 - t1);
	w += tw::Float(t > t3 && t <= t4) / (t4 - t3);

	switch (whichProfile)
	{
		case tw::profile::shape::triangle:
			return 0.0;
		case tw::profile::shape::sech:
			return w*D1SechRise(tau);
		case tw::profile::shape::quartic:
			return w*D1QuarticRise(tau);
		case tw::profile::shape::quintic:
			return w*D1QuinticRise(tau);
		case tw::profile::shape::sin2:
      return w*D1Sin2Rise(tau);
	}
	return 0.0;
}

tw::Float PulseShape::D1Intensity(const tw::Float t) const
{
    return 2.0*PulseShapeFactor(t)*D1Amplitude(t);
}


/////////////////////////
//                     //
// Radiation Injection //
//                     //
/////////////////////////


Wave::Wave(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	direction = tw::vec3(0.0,0.0,1.0);
	focusPosition = tw::vec3(0.0,0.0,0.0);
	a = tw::vec3(0.1,0.0,0.0);
	nrefr = 1.0;
	w = 1.0;
	chirp = 0.0;
	phase = 0.0;
	randomPhase = 0.0;
	gammaBoost = 1.0;
	zones = 1;
	modeData.order[0] = 0;
	modeData.order[1] = 0;
	modeData.scale[0] = 1.0;
	modeData.scale[1] = 1.0;
	modeData.exponent[0] = 2;
	modeData.exponent[1] = 2;
	directives.Add("direction",new tw::input::Vec3(&direction),false);
	directives.Add("focus position",new tw::input::Vec3(&focusPosition),false);
	directives.Add("a",new tw::input::Vec3(&a));
	directives.Add("w",new tw::input::Float(&w));
	directives.Add("refractive index",new tw::input::Float(&nrefr),false);
	directives.Add("chirp",new tw::input::Float(&chirp),false);
	directives.Add("phase",new tw::input::Float(&phase),false);
	directives.Add("random phase",new tw::input::Float(&randomPhase),false);
	directives.Add("delay",new tw::input::Float(&pulseShape.delay));
	directives.Add("risetime",new tw::input::Float(&pulseShape.risetime));
	directives.Add("holdtime",new tw::input::Float(&pulseShape.holdtime));
	directives.Add("falltime",new tw::input::Float(&pulseShape.falltime));
	directives.Add("boosted frame gamma",new tw::input::Float(&gammaBoost),false);
	directives.Add("exponent",new tw::input::Numbers<tw::Int>(&modeData.exponent[0],2),false);
	directives.Add("mode",new tw::input::Numbers<tw::Int>(&modeData.order[0],2),false);
	std::map<std::string,tw::profile::shape> shape = {{"quintic",tw::profile::shape::quintic},{"sech",tw::profile::shape::sech},{"sin2",tw::profile::shape::sin2}};
	directives.Add("shape",new tw::input::Enums<tw::profile::shape>(shape,&pulseShape.whichProfile),false);
	directives.Add("zones",new tw::input::Int(&zones),false);
}

void Wave::Initialize()
{
	if (pulseShape.risetime<=0.0)
		throw tw::FatalError("Pulse rise time must be positive and non-zero.");
	if (pulseShape.holdtime<0.0)
		throw tw::FatalError("Pulse hold time must be positive.");
	if (pulseShape.falltime<=0.0)
		throw tw::FatalError("Pulse fall time must be positive and non-zero.");
	// set up a transformation so we can work in a coordinate system
	// where the polarization direction is x and the propagation direction is z
	vg = nrefr>1.0 ? 1.0/nrefr : nrefr;
	laserFrame.w = direction;
	laserFrame.v = direction | a;
	laserFrame.u = laserFrame.v | laserFrame.w;
	Normalize(laserFrame);
	a0 = Magnitude(a);
	pulseShape.Initialize(pulseShape.delay + pulseShape.risetime);
}

PlaneWave::PlaneWave(const std::string& name,MetricSpace *m,Task *tsk) : Wave(name,m,tsk)
{
}

tw::Complex PlaneWave::PrimitivePhasor(const tw::vec4& x) const
{
	// t,r are expected to be in the Cartesian laser basis for all primitive functions
	const tw::Float tau = x[0] - nrefr*x[3];
	const tw::Float psi = phase - w*tau - chirp*tau*tau;
	return a0 * pulseShape.PulseShapeFactor(tau) * std::exp(ii*psi);
}

tw::vec3 PlaneWave::PrimitiveVector(const tw::vec4& x4) const
{
	return tw::vec3(std::real(PrimitivePhasor(x4)),0.0,0.0);
}

BesselBeam::BesselBeam(const std::string& name,MetricSpace *m,Task *tsk) : Wave(name,m,tsk)
{
	directives.Add("r0",new tw::input::Numbers<tw::Float>(&modeData.scale[0],2));
}

tw::Complex BesselBeam::PrimitivePhasor(const tw::vec4& x) const
{
	const tw::Float rbar = sqrt(sqr(x[1]) + sqr(x[2]))/modeData.scale[0];
	const tw::Float tau = x[0] - nrefr*x[3];
	const tw::Float psi = phase - w*tau - chirp*tau*tau;
	return a0 * tw::cyl_bessel_j(0,rbar) * pulseShape.PulseShapeFactor(tau) * std::exp(ii*psi);
}

AiryDisc::AiryDisc(const std::string& name,MetricSpace *m,Task *tsk) : Wave(name,m,tsk)
{
	directives.Add("r0",new tw::input::Numbers<tw::Float>(&modeData.scale[0],2));
}

tw::Complex AiryDisc::PrimitivePhasor(const tw::vec4& x) const
{
	// User input is the radius of the first zero.
	// This is convenient because then the energy corresponds to a Gaussian with the same radius.
	const tw::Float rbar = sqrt(sqr(x[1]) + sqr(x[2]))*3.83171/modeData.scale[0];
	const tw::Float tau = x[0] - nrefr*x[3];
	const tw::Float psi = phase - w*tau - chirp*tau*tau;
	return 2*a0 * (tw::cyl_bessel_j(1,rbar)/rbar) * pulseShape.PulseShapeFactor(tau) * std::exp(ii*psi);
}

LaguerreGauss::LaguerreGauss(const std::string& name,MetricSpace *m,Task *tsk) : Wave(name,m,tsk)
{
	directives.Add("r0",new tw::input::Numbers<tw::Float>(&modeData.scale[0],2));
}

void LaguerreGauss::Initialize()
{
	Wave::Initialize();
	if (w==0.0)
		throw tw::FatalError("Zero frequency not allowed for Laguerre-Gauss mode.");
}

tw::Complex LaguerreGauss::PrimitivePhasor(const tw::vec4& x) const
{
	const tw::Float r0 = modeData.scale[0];
	const tw::Int nexp = modeData.exponent[0];
	const tw::Int mu = modeData.order[1];
	const tw::Int mv = modeData.order[1];

	tw::Float Ax=a0, tau=x[0]-nrefr*x[3], guoy_shift=0.0;
	tw::Float zR,rm,rho,phi;

	// setup for Laguerre-Gaussian mode
	rho = sqrt(x[1]*x[1] + x[2]*x[2]);
	phi = atan2(x[2],x[1]);
	zR = 0.5*nrefr*w*r0*r0;
	rm = r0*sqrt(one + x[3]*x[3]/(zR*zR));
	guoy_shift = -(one + two*mu + mv)*atan(x[3]/zR);

	// compute factors for r-phi profile
	Ax *= r0/rm;
	Ax *= tw::assoc_laguerre(mu,mv,two*sqr(rho/rm));
	Ax *= nexp%2==0 ? exp(-pow(rho/rm,nexp)) : pow(cos(0.5*pi*rho/rm),nexp+1)*tw::Float(rho<rm);
	Ax *= pow(sqrt(two)*rho/rm,tw::Float(mv)); // put appropriate hole on axis
	Ax *= sqrt(Factorial(mu)/Factorial(mu+mv)); // Laguerre normalization
	tau += guoy_shift/w - 0.5*nrefr*rho*rho*x[3]/(x[3]*x[3] + zR*zR);
	tw::Float psi = phase + 2.0*guoy_shift - mv*phi - w*tau - chirp*tau*tau;

	return Ax * pulseShape.PulseShapeFactor(tau) * std::exp(ii*psi);
}

HermiteGauss::HermiteGauss(const std::string& name,MetricSpace *m,Task *tsk) : Wave(name,m,tsk)
{
	directives.Add("r0",new tw::input::Numbers<tw::Float>(&modeData.scale[0],2));
}

void HermiteGauss::Initialize()
{
	Wave::Initialize();
	if (w==0.0)
		throw tw::FatalError("Zero frequency not allowed for Hermite-Gauss mode.");
}

tw::Complex HermiteGauss::PrimitivePhasor(const tw::vec4& x) const
{
	const tw::Float r0[2] = { modeData.scale[0] , modeData.scale[1] };
	const tw::Int nexp[2] = { modeData.exponent[0] , modeData.exponent[1] };
	const tw::Int order[2] = { modeData.order[0] , modeData.order[1] };

	tw::Float Ax=a0, tau=x[0]-nrefr*x[3], guoy_shift=0.0;

	// lambda function to compute Hermite amplitude and phase per axis
	auto herm = [&] (tw::Float q,tw::Int ax)
	{
		const tw::Float m = order[ax-1];
		const tw::Float r00 = r0[ax-1];
		const tw::Float zR = 0.5*nrefr*w*r00*r00;
		const tw::Float rm = r00*sqrt(1.0 + x[3]*x[3]/(zR*zR));
		guoy_shift -= (0.5 + m)*atan(x[3]/zR);
		Ax *= sqrt(r00/rm);
		Ax *= tw::hermite(m,1.41421356*q/rm);
		Ax *= nexp[ax-1]%2==0 ? exp(-pow(q/rm,nexp[ax-1])) : pow(cos(0.5*pi*q/rm),nexp[ax-1]+1)*tw::Float(q*q<rm*rm);
		Ax /= sqrt(pow(two,tw::Float(m))*Factorial(m)); // Hermite normalization
		tau += guoy_shift/w - 0.5*nrefr*q*q*x[3]/(x[3]*x[3] + zR*zR);
	};

	herm(x[1],1);
	herm(x[2],2);
	tw::Float psi = phase + 2.0*guoy_shift - w*tau - chirp*tau*tau;

	return Ax * pulseShape.PulseShapeFactor(tau) * std::exp(ii*psi);
}

Multipole::Multipole(const std::string& name,MetricSpace *m,Task *tsk) : Wave(name,m,tsk)
{
}

void Multipole::Initialize()
{
	Wave::Initialize();
	if (w==0.0)
		throw tw::FatalError("Zero frequency not allowed for multipole mode.");
}

tw::Complex Multipole::PrimitivePhasor(const tw::vec4& x) const
{
	// For now hard code in magnetic dipole radiation
	const tw::Float j1max = 0.436182;
	const tw::Float rho = sqrt(x[1]*x[1] + x[2]*x[2]);
	const tw::Float R = sqrt(rho*rho + x[3]*x[3]);
	const tw::Float stheta = rho/R;
	const tw::Float wR = w*R;
	//auto spherical_bessel = [&] (tw::Float x) { return sin(x)/sqr(x+tw::small_pos) - cos(x)/(x+tw::small_pos); };

	tw::Complex Aphi;
	tw::Float tau_in,tau_out;
	tau_out = x[0] - nrefr*R;
	tau_in = x[0] + nrefr*R;
	Aphi = (-one/wR + one/(ii*wR*wR)) * pulseShape.PulseShapeFactor(tau_out) * std::exp(-ii*w*tau_out);
	Aphi += (-one/wR - one/(ii*wR*wR)) * pulseShape.PulseShapeFactor(tau_in) * std::exp(-ii*w*tau_in);
	Aphi *= 0.5 * a0 * (stheta/j1max);
	return Aphi;
}

tw::vec3 Multipole::PrimitiveVector(const tw::vec4& x) const
{
	tw::Complex Ax = PrimitivePhasor(x);
	const tw::Float rho = tw::small_pos + sqrt(x[1]*x[1] + x[2]*x[2]);
	return tw::vec3( -std::real(Ax)*x[2]/rho , std::real(Ax)*x[1]/rho ,0.0 );
}


////////////////////
//                //
//   CONDUCTOR    //
//                //
////////////////////


Conductor::Conductor(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	affectsPhi = true;
	affectsA = true;
	currentType = EM::current::electric;
	pulseShape.delay = 0.0;
	pulseShape.risetime = 0.01*timestep(*m);
	pulseShape.holdtime = 1e10*timestep(*m);
	pulseShape.falltime = 0.01*timestep(*m);
	gaussianRadius = tw::big_pos;
	f = tw::big_pos;
	ks = 0.0;
	temperature = 0.0;
	directives.Add("px",new tw::input::List<std::valarray<tw::Float>>(&Px),false);
	directives.Add("py",new tw::input::List<std::valarray<tw::Float>>(&Py),false);
	directives.Add("pz",new tw::input::List<std::valarray<tw::Float>>(&Pz),false);
	directives.Add("potential",new tw::input::List<std::valarray<tw::Float>>(&potential),false);
	directives.Add("w",new tw::input::List<std::valarray<tw::Float>>(&angFreq),false);
	directives.Add("phase",new tw::input::List<std::valarray<tw::Float>>(&phase),false);
	directives.Add("delay",new tw::input::Float(&pulseShape.delay),false);
	directives.Add("risetime",new tw::input::Float(&pulseShape.risetime),false);
	directives.Add("holdtime",new tw::input::Float(&pulseShape.holdtime),false);
	directives.Add("falltime",new tw::input::Float(&pulseShape.falltime),false);
	std::map<std::string,tw::profile::shape> shapeMap = {{"quintic",tw::profile::shape::quintic},{"sech",tw::profile::shape::sech},{"sin2",tw::profile::shape::sin2}};
	directives.Add("shape",new tw::input::Enums<tw::profile::shape>(shapeMap,&pulseShape.whichProfile),false);
	directives.Add("enable electrostatic",new tw::input::Bool(&affectsPhi),false);
	directives.Add("enable electromagnetic",new tw::input::Bool(&affectsA),false);
	std::map<std::string,EM::current> currentMap = {{"none",EM::current::none},{"electric",EM::current::electric},{"magnetic",EM::current::magnetic}};
	directives.Add("current type",new tw::input::Enums<EM::current>(currentMap,&currentType),false);
	directives.Add("gaussian size",new tw::input::Vec3(&gaussianRadius),false);
	directives.Add("f",new tw::input::Float(&f),false);
	directives.Add("ks",new tw::input::Vec3(&ks),false);
	directives.Add("temperature",new tw::input::Float(&temperature),false);
}

void Conductor::Initialize()
{
	if (currentType!=EM::current::electric)
		throw tw::FatalError("Current type must be electric in this version.");
	if (pulseShape.risetime<=0.0)
		throw tw::FatalError("Pulse rise time must be positive and non-zero.");
	if (pulseShape.holdtime<0.0)
		throw tw::FatalError("Pulse hold time must be positive.");
	if (pulseShape.falltime<=0.0)
		throw tw::FatalError("Pulse fall time must be positive and non-zero.");

	tw::Int maxArraySize;
	pulseShape.Initialize(0.0);

	maxArraySize = Px.size();
	maxArraySize = Py.size()>maxArraySize ? Py.size() : maxArraySize;
	maxArraySize = Pz.size()>maxArraySize ? Py.size() : maxArraySize;
	maxArraySize = potential.size()>maxArraySize ? Py.size() : maxArraySize;
	maxArraySize = angFreq.size()>maxArraySize ? Py.size() : maxArraySize;
	maxArraySize = phase.size()>maxArraySize ? Py.size() : maxArraySize;
	maxArraySize = maxArraySize==0 ? 1 : maxArraySize;

	if (Px.size()==0)
	{
		Px.resize(maxArraySize);
		Px = 0.0;
	}
	if (Py.size()==0)
	{
		Py.resize(maxArraySize);
		Py = 0.0;
	}
	if (Pz.size()==0)
	{
		Pz.resize(maxArraySize);
		Pz = 0.0;
	}
	if (potential.size()==0)
	{
		potential.resize(maxArraySize);
		potential = 0.0;
	}
	if (angFreq.size()==0)
	{
		angFreq.resize(maxArraySize);
		angFreq = 0.0;
	}
	if (phase.size()==0)
	{
		phase.resize(maxArraySize);
		phase = 0.0;
	}
}

tw::Float Conductor::Voltage(tw::Float t)
{
	tw::Float ans = 0.0;
	for (tw::Int i=0;i<potential.size();i++)
		ans += pulseShape.PulseShapeFactor(t)*potential[i]*cos(angFreq[i]*t + phase[i]);
	return ans;
}

tw::Float Conductor::VoltageRate(tw::Float t)
{
	tw::Float ans = 0.0;
	for (tw::Int i=0;i<potential.size();i++)
	{
		ans -= pulseShape.PulseShapeFactor(t)*angFreq[i]*potential[i]*sin(angFreq[i]*t + phase[i]);
		ans += pulseShape.D1Amplitude(t)*potential[i]*cos(angFreq[i]*t + phase[i]);
	}
	return ans;
}

tw::vec3 Conductor::PolarizationDensity(const tw::vec3& pos,tw::Float t)
{
	tw::vec3 P0;
	tw::vec3 ans(0.0,0.0,0.0);
	tw::vec3 rc = pos - theRgn->center;
	theRgn->orientation.ExpressInBasis(&rc);
	for (tw::Int i=0;i<potential.size();i++)
	{
		P0 = tw::vec3(Px[i],Py[i],Pz[i]);
		ans += P0*sin(angFreq[i]*t + phase[i] + (0.5*angFreq[i]/f)*(rc.x*rc.x + rc.y*rc.y) + (ks^rc));
	}
	ans *= pulseShape.PulseShapeFactor(t + (0.5/f)*(rc.x*rc.x + rc.y*rc.y));
	ans *= exp(-sqr(rc.x/gaussianRadius.x)-sqr(rc.y/gaussianRadius.y)-sqr(rc.z/gaussianRadius.z));
	theRgn->orientation.ExpressInStdBasis(&ans);
	return ans;
}

void Conductor::DepositSources(Field& sources,tw::Float t,tw::Float dt)
{
	tw::Int x0,x1,y0,y1,z0,z1;
	const MetricSpace& m = *space;
	theRgn->GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (tw::Int i=x0;i<=x1;i++)
		for (tw::Int j=y0;j<=y1;j++)
			for (tw::Int k=z0;k<=z1;k++)
			{
				const tw::vec3 pos = m.Pos(i,j,k);
				if (theRgn->Inside(pos,m))
				{
					const tw::Float dV = m.dS(i,j,k,0);
					// To conserve charge we have to center arc lengths in arc direction
					const tw::vec3 dl = 0.5 * tw::vec3(m.dL(i,j,k,1),m.dL(i,j,k,2),m.dL(i,j,k,3));
					const tw::vec3 Q0 = 0.5 * dV * PolarizationDensity(pos,t-dt) / dl;
					const tw::vec3 Q1 = 0.5 * dV * PolarizationDensity(pos,t) / dl;
					const tw::vec3 I3 = (Q1 - Q0)/dt;

					// Deposit cell charge using small displacement model
					sources(i-1,j,k,0) -= Q1.x;
					sources(i,j-1,k,0) -= Q1.y;
					sources(i,j,k-1,0) -= Q1.z;
					sources(i,j,k+1,0) += Q1.z;
					sources(i,j+1,k,0) += Q1.y;
					sources(i+1,j,k,0) += Q1.x;

					// Deposit wall current using small displacement model
					sources(i,j,k,1) += I3.x;
					sources(i+1,j,k,1) += I3.x;
					sources(i,j,k,2) += I3.y;
					sources(i,j+1,k,2) += I3.y;
					sources(i,j,k,3) += I3.z;
					sources(i,j,k+1,3) += I3.z;
				}
			}
}

//////////////////////
//                  //
// LINDMAN BOUNDARY //
//                  //
//////////////////////


void LindmanBoundary::Initialize(Task *task,MetricSpace *ms,std::vector<Wave*> *wave,const tw::grid::axis& axis,const tw::grid::side& side,tw::Int c0,tw::Int c1)
{
	this->ms = ms;
	this->wave = wave;
	this->axis = axis;
	this->side = side;
	// Set() assumes that c0 is an ordered 4-vector index
	this->c0 = c0;
	this->c1 = c1;

	DiscreteSpace bm_layout;

	if (axis==tw::grid::x && ms->Dim(1)==1)
		throw tw::FatalError("Lindman boundary geometry is not consistent");

	if (axis==tw::grid::y && ms->Dim(2)==1)
		throw tw::FatalError("Lindman boundary geometry is not consistent");

	if (axis==tw::grid::z && ms->Dim(3)==1)
		throw tw::FatalError("Lindman boundary geometry is not consistent");

	switch (axis)
	{
		case tw::grid::x:
			bm_layout.Resize(1,ms->Dim(2),ms->Dim(3),Corner(*ms),PhysicalSize(*ms));
			break;
		case tw::grid::y:
			bm_layout.Resize(ms->Dim(1),1,ms->Dim(3),Corner(*ms),PhysicalSize(*ms));
			break;
		case tw::grid::z:
			bm_layout.Resize(ms->Dim(1),ms->Dim(2),1,Corner(*ms),PhysicalSize(*ms));
			break;
		default:
			throw tw::FatalError("Lindman boundary axis must be one of x,y,z");
	}

	boundaryMemory.Initialize( 9 , bm_layout , task );
	boundaryMemory.SetBoundaryConditions(tw::grid::x,fld::none,fld::none);
	boundaryMemory.SetBoundaryConditions(tw::grid::y,fld::none,fld::none);
	boundaryMemory.SetBoundaryConditions(tw::grid::z,fld::none,fld::none);
}

void LindmanBoundary::UpdateBoundaryMemory(Field& A,tw::Float dt)
{
	tw::Float alpha[3] = {0.3264,0.1272,0.0309};
	tw::Float beta[3] = {0.7375,0.98384,0.9996472};

	tw::Int ax = tw::grid::naxis(axis);
	tw::Float source,dtt;
	tw::Int i,j,k,s,s0,s1;

	s0 = side==tw::grid::low ? 0 : A.Dim(ax);
	s1 = side==tw::grid::low ? 1 : A.Dim(ax)+1;

	for (auto strip : StripRange(A,ax,strongbool::no))
	{
		source = A.d2(strip,s1,c0,1) + A.d2(strip,s1,c0,2) + A.d2(strip,s1,c0,3) - A.d2(strip,s1,c0,ax);
		source -= A.d2(strip,s0,c0,1) + A.d2(strip,s0,c0,2) + A.d2(strip,s0,c0,3) - A.d2(strip,s0,c0,ax);
		source = A.d2(strip,s1,c1,1) + A.d2(strip,s1,c1,2) + A.d2(strip,s1,c1,3) - A.d2(strip,s1,c1,ax);
		source -= A.d2(strip,s0,c1,1) + A.d2(strip,s0,c1,2) + A.d2(strip,s0,c1,3) - A.d2(strip,s0,c1,ax);
		strip.Decode(0,&i,&j,&k);
		for (s=0;s<3;s++)
		{
			dtt =  dxi(A)*dxi(A)*(boundaryMemory(i-1,j,k,3+s) - 2.0*boundaryMemory(i,j,k,3+s) + boundaryMemory(i+1,j,k,3+s));
			dtt += dyi(A)*dyi(A)*(boundaryMemory(i,j-1,k,3+s) - 2.0*boundaryMemory(i,j,k,3+s) + boundaryMemory(i,j+1,k,3+s));
			dtt += dzi(A)*dzi(A)*(boundaryMemory(i,j,k-1,3+s) - 2.0*boundaryMemory(i,j,k,3+s) + boundaryMemory(i,j,k+1,3+s));
			dtt =  dt*dt*(beta[s]*dtt + alpha[s]*source);
			boundaryMemory(i,j,k,6+s) = 2.0*boundaryMemory(i,j,k,3+s) - boundaryMemory(i,j,k,s) + dtt;
		}
	}

	//boundaryMemory.ApplyBoundaryCondition(Element(6,8));
	boundaryMemory.CopyFromNeighbors(Element(6,8));

	boundaryMemory.Swap(Element(0,2),Element(3,5));
	boundaryMemory.Swap(Element(3,5),Element(6,8));
}

void LindmanBoundary::Set(Field& A,tw::Float t0,tw::Float dt)
{
	tw::Int ax = tw::grid::naxis(axis);
	tw::Float sgn,correction,ds[4] = {dt,dx(A),dy(A),dz(A)};
	tw::Int i,j,k,s,s0,s1;
	tw::vec3 pos,Ain;

	tw::Float propagator = (1.0 - ds[0]/ds[ax])/(1.0 + ds[0]/ds[ax]);
	tw::Float injector = 1.0/(1.0 + ds[0]/ds[ax]);

	tw::vec3 direction(tw::Int(ax==1),tw::Int(ax==2),tw::Int(ax==3));
	tw::vec3 offset = 0.5 * direction * ds[ax];
	sgn = side==tw::grid::high ? -1.0 : 1.0;

	s0 = side==tw::grid::low ? 0 : A.Dim(ax)+1;
	s1 = side==tw::grid::low ? 1 : A.Dim(ax);
	for (auto strip : StripRange(A,ax,strongbool::no))
	{
		strip.Decode(s0,&i,&j,&k);
		correction = boundaryMemory(i,j,k,3) + boundaryMemory(i,j,k,4) + boundaryMemory(i,j,k,5);
		pos = ms->Pos(i,j,k) + sgn*offset;
		Ain = 0.0;
		for (s=0;s<wave->size();s++)
		{
			if ( (((*wave)[s]->direction) ^ (sgn*direction)) > 0.0)
			{
				Ain += 4.0*(*wave)[s]->VectorPotential(t0+dt,pos);
				Ain -= 4.0*(*wave)[s]->VectorPotential(t0,pos);
			}
		}
		A(strip,s0,c0) = A(strip,s1,c1) + propagator*(A(strip,s0,c1)-A(strip,s1,c0));
		A(strip,s0,c0) += injector*(Ain[c0-1] + sgn*dt*correction*ds[ax]);
	}
}

void LindmanBoundary::ReadCheckpoint(std::ifstream& inFile)
{
	inFile.read((char*)&c0,sizeof(tw::Int));
	inFile.read((char*)&c1,sizeof(tw::Int));
	boundaryMemory.ReadCheckpoint(inFile);
}

void LindmanBoundary::WriteCheckpoint(std::ofstream& outFile)
{
	outFile.write((char*)&c0,sizeof(tw::Int));
	outFile.write((char*)&c1,sizeof(tw::Int));
	boundaryMemory.WriteCheckpoint(outFile);
}


//////////////////////////////
//                          //
// Mora & Antonsen BOUNDARY //
//                          //
//////////////////////////////


MABoundary::MABoundary(tw::Int numCells)
{
	tw::Int i;
	for (i=0;i<10;i++)
	{
		g[i].resize(numCells);
		g[i] = tw::Complex(0.0,0.0);
	}

	gamma[0] = .05;
	gamma[1] = .0813625303;
	gamma[2] = .132397227;
	gamma[3] = .215443469;
	gamma[4] = .35058516;
	gamma[5] = .570482359;
	gamma[6] = .928317767;
	gamma[7] = 1.51060565;
	gamma[8] = 2.45813397;
	gamma[9] = 4.0;

	sigma[0] = -.0530263170;
	sigma[1] = .188607311;
	sigma[2] = -.386762908;
	sigma[3] = .479633721;
	sigma[4] = -.39959789;
	sigma[5] = -.409263775;
	sigma[6] = 2.41066015;
	sigma[7] = -8.90422832;
	sigma[8] = 20.6407061;
	sigma[9] = -30.0417443;

	scaleFactor = 1.0;
}

void MABoundary::AdvanceLeft(std::valarray<tw::Complex>& amplitude,tw::Float dt)
{
	const tw::Float dts = scaleFactor*dt;  // scale factor is A&M's "w"
	tw::Int i,j;
	for (i=0;i<10;i++)
		for (j=0;j<g[i].size();j++)
			g[i][j] = (g[i][j]/dts - sigma[i]*amplitude[j])/(gamma[i] + one/dts);
}

void MABoundary::AdvanceRight(std::valarray<tw::Complex>& amplitude,tw::Float dt)
{
	const tw::Float dts = scaleFactor*dt;  // scale factor is A&M's "w"
	tw::Int i,j;
	for (i=0;i<10;i++)
		for (j=0;j<g[i].size();j++)
			g[i][j] = (g[i][j]/dts + sigma[i]*amplitude[j])/(gamma[i] + one/dts);
}

tw::Complex MABoundary::NormalDerivativeLeft(tw::Int j,tw::Complex amplitude,tw::Float carrierFrequency)
{
	tw::Int i;
	tw::Complex ans = 0.0;
	for (i=0;i<10;i++)
		ans += g[i][j] + amplitude*sigma[i]/gamma[i];
	return ans*std::sqrt(-two*ii*carrierFrequency*scaleFactor);
}

tw::Complex MABoundary::NormalDerivativeRight(tw::Int j,tw::Complex amplitude,tw::Float carrierFrequency)
{
	tw::Int i;
	tw::Complex ans = 0.0;
	for (i=0;i<10;i++)
		ans += g[i][j] - amplitude*sigma[i]/gamma[i];
	return ans*std::sqrt(-two*ii*carrierFrequency*scaleFactor);
}

//////////////////
//  GRID WARPS  //
//////////////////

Warp::Warp(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	ax = tw::grid::z;
	increasing = true;
	directives.Add("axis",new tw::input::Enums<tw::grid::axis>(tw::grid::axis_map(),&ax));
	directives.Add("increasing",new tw::input::Bool(&increasing));
	directives.Add("index range",new tw::input::Numbers<tw::Int>(&rng[0],2));
	directives.Add("length",new tw::input::Float(&L));
}

void Warp::Initialize()
{
	const tw::Int N = rng[1] - rng[0] + 1;
	gridSum = 0.0;
	for (tw::Int i=1;i<=N;i++)
		gridSum += QuinticRise(tw::Float(i-1)/tw::Float(N-1));
}

tw::Float Warp::AddedCellWidth(tw::Int globalCell)
{
	const tw::Int N = rng[1] - rng[0] + 1;
	const tw::Float h = space->dx0(tw::grid::naxis(ax));
	const tw::Float A = (1.0/gridSum)*(L/h - N);
	if (globalCell>=rng[0] && globalCell<=rng[1])
	{
		if (increasing)
			return h*A*QuinticRise(tw::Float(globalCell-rng[0])/tw::Float(N-1));
		else
			return h*A*QuinticFall(tw::Float(globalCell-rng[0])/tw::Float(N-1));
	}
	else
		return 0.0;
}

void Warp::ReadCheckpoint(std::ifstream& inFile)
{
	ComputeTool::ReadCheckpoint(inFile);
	inFile.read((char*)&ax,sizeof(ax));
	inFile.read((char*)&increasing,sizeof(increasing));
	inFile.read((char*)&rng[0],sizeof(rng));
	inFile.read((char*)&L,sizeof(L));
	inFile.read((char*)&gridSum,sizeof(gridSum));
}

void Warp::WriteCheckpoint(std::ofstream& outFile)
{
	ComputeTool::WriteCheckpoint(outFile);
	outFile.write((char*)&ax,sizeof(ax));
	outFile.write((char*)&increasing,sizeof(increasing));
	outFile.write((char*)&rng[0],sizeof(rng));
	outFile.write((char*)&L,sizeof(L));
	outFile.write((char*)&gridSum,sizeof(gridSum));
}
