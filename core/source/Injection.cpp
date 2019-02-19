#include "simulation.h"

///////////////////
//               //
// PROFILE CLASS //
//               //
///////////////////


Profile::Profile(std::vector<Region*>& rlist) : rgnList(rlist)
{
	profileSpec = nullProfile;
	symmetry = cartesian;
	segmentShape = loading::triangle;
	neutralize = true;
	timingMethod = triggeredProfile;
	variableCharge = false;
	loadingMethod = deterministic;
	whichQuantity = densityProfile;
	modeAmplitude = 0.0;
	modeNumber = 0.0;
	temperature = 0.0;
	theRgn = rgnList[0];
	wasTriggered = false;
	t0 = 0.0;
	t1 = tw::big_pos;
	orientation.u = tw::vec3(1,0,0);
	orientation.v = tw::vec3(0,1,0);
	orientation.w = tw::vec3(0,0,1);
}

Profile* Profile::CreateObjectFromFile(std::vector<Region*>& rgnList,std::ifstream& inFile)
{
	Profile *ans;
	tw_profile_spec theSpec;
	inFile.read((char *)&theSpec,sizeof(tw_profile_spec));
	switch (theSpec)
	{
		case nullProfile:
			break;
		case uniformProfile:
			ans = new UniformProfile(rgnList);
			break;
		case gaussianProfile:
			ans = new GaussianProfile(rgnList);
			break;
		case channelProfile:
			ans = new ChannelProfile(rgnList);
			break;
		case columnProfile:
			ans = new ColumnProfile(rgnList);
			break;
		case piecewiseProfile:
			ans = new PiecewiseProfile(rgnList);
			break;
		case corrugatedProfile:
			ans = new CorrugatedProfile(rgnList);
			break;
	}
	ans->ReadData(inFile);
	return ans;
}

void Profile::ReadInputFileBlock(std::stringstream& inputString,bool neutralize)
{
	this->neutralize = neutralize;
	std::string com;
	do
	{
		inputString >> com;
		ReadInputFileDirective(inputString,com);
	} while (com!="}");
}

void Profile::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
    std::string word;

	if (com=="use")
		throw tw::FatalError("'use' directives not supported.  See docs for alternative.");

	theRgn = Region::ReadRegion(rgnList,theRgn,inputString,com);

	if (com=="position")
		inputString >> word >> centerPt.x >> centerPt.y >> centerPt.z;
	if (com=="type") // eg, type = energy
	{
		inputString >> word >> word;
		if (word=="density")
			whichQuantity = densityProfile;
		if (word=="energy")
			whichQuantity = energyProfile;
		if (word=="px")
			whichQuantity = pxProfile;
		if (word=="py")
			whichQuantity = pyProfile;
		if (word=="pz")
			whichQuantity = pzProfile;
	}
	if (com=="drift")
		inputString >> word >> word >> driftMomentum.x >> driftMomentum.y >> driftMomentum.z;
	if (com=="thermal")
		inputString >> word >> word >> thermalMomentum.x >> thermalMomentum.y >> thermalMomentum.z;
	if (com=="temperature_eV" || com=="temperature_ev")
		throw tw::FatalError("Outdated temperature input.  Use unit conversion macros, e.g., temperature = %1.0eV or temperature = %11600K");
	if (com=="temperature")
		inputString >> word >> temperature;
	if (com=="shape") // eg, shape = triangle
	{
		inputString >> word >> word;
		if (word=="triangle")
			segmentShape = loading::triangle;
		if (word=="quartic")
			segmentShape = loading::quartic;
		if (word=="quintic")
			segmentShape = loading::quintic;
	}
	if (com=="timing")
	{
		inputString >> word >> word;
		if (word=="trigger" || word=="triggered")
			timingMethod = triggeredProfile;
		else
			timingMethod = maintainedProfile;
	}
    if (com=="t0")
        inputString >> word >> t0;
    if (com=="t1")
        inputString >> word >> t1;
	if (com=="loading")
	{
		inputString >> word >> word;
		if (word=="deterministic")
			loadingMethod = deterministic;
		else
			loadingMethod = statistical;
	}
	if (com=="symmetry")
	{
		inputString >> word >> word;
		if (word=="axisymmetric" || word=="cylindrical")
			symmetry = cylindrical;
		if (word=="spherical")
			symmetry = spherical;
	}
	if (com=="particle")
	{
		inputString >> word >> word >> word;
		if (word=="variable")
			variableCharge = true;
		else
			variableCharge = false;
	}
	if (com=="mode")
	{
		inputString >> word;
		if (word=="amplitude") // eg, mode amplitude = 0.1
			inputString >> word >> modeAmplitude;
		if (word=="number") // eg, mode number = ( 0 1 1 )
			inputString >> word >> modeNumber.x >> modeNumber.y >> modeNumber.z;
	}
	if (com=="euler") // eg, euler angles = ( 45 90 30 )
	{
		tw::Float alpha,beta,gamma;
		inputString >> word >> word >> alpha >> beta >> gamma;
		orientation.SetWithEulerAngles(alpha*pi/180.0,beta*pi/180.0,gamma*pi/180.0);
	}
}

void Profile::ReadData(std::ifstream& inFile)
{
    tw::Int rgnIndex;
    inFile.read((char *)&rgnIndex,sizeof(tw::Int));
    theRgn = rgnList[rgnIndex];
	inFile.read((char *)&centerPt,sizeof(tw::vec3));
	inFile.read((char *)&symmetry,sizeof(tw_geometry));
	inFile.read((char *)&segmentShape,sizeof(segmentShape));
	inFile.read((char *)&timingMethod,sizeof(tw_profile_timing));
	inFile.read((char *)&thermalMomentum,sizeof(tw::vec3));
	inFile.read((char *)&driftMomentum,sizeof(tw::vec3));
	inFile.read((char *)&neutralize,sizeof(bool));
	inFile.read((char *)&variableCharge,sizeof(bool));
	inFile.read((char *)&loadingMethod,sizeof(tw_load_method));
	inFile.read((char *)&whichQuantity,sizeof(tw_profile_quantity));
	inFile.read((char *)&modeNumber,sizeof(tw::vec3));
	inFile.read((char *)&modeAmplitude,sizeof(tw::Float));
	inFile.read((char *)&temperature,sizeof(tw::Float));
	inFile.read((char *)&orientation,sizeof(orientation));
	inFile.read((char *)&wasTriggered,sizeof(bool));
	inFile.read((char *)&t0,sizeof(tw::Float));
	inFile.read((char *)&t1,sizeof(tw::Float));
}

void Profile::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&profileSpec,sizeof(tw_profile_spec)); // must be first
    tw::Int rgnIndex = std::find(rgnList.begin(),rgnList.end(),theRgn) - rgnList.begin();
	outFile.write((char *)&rgnIndex,sizeof(tw::Int));
	outFile.write((char *)&centerPt,sizeof(tw::vec3));
	outFile.write((char *)&symmetry,sizeof(tw_geometry));
	outFile.write((char *)&segmentShape,sizeof(segmentShape));
	outFile.write((char *)&timingMethod,sizeof(tw_profile_timing));
	outFile.write((char *)&thermalMomentum,sizeof(tw::vec3));
	outFile.write((char *)&driftMomentum,sizeof(tw::vec3));
	outFile.write((char *)&neutralize,sizeof(bool));
	outFile.write((char *)&variableCharge,sizeof(bool));
	outFile.write((char *)&loadingMethod,sizeof(tw_load_method));
	outFile.write((char *)&whichQuantity,sizeof(tw_profile_quantity));
	outFile.write((char *)&modeNumber,sizeof(tw::vec3));
	outFile.write((char *)&modeAmplitude,sizeof(tw::Float));
	outFile.write((char *)&temperature,sizeof(tw::Float));
	outFile.write((char *)&orientation,sizeof(orientation));
	outFile.write((char *)&wasTriggered,sizeof(bool));
	outFile.write((char *)&t0,sizeof(tw::Float));
	outFile.write((char *)&t1,sizeof(tw::Float));
}

tw::Float Profile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	return theRgn->Inside(pos,ds) ? 1.0 : 0.0;
}

void UniformProfile::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
    std::string word;
	Profile::ReadInputFileDirective(inputString,com);
	if (com=="density")
		inputString >> word >> density;
}

void UniformProfile::ReadData(std::ifstream& inFile)
{
	Profile::ReadData(inFile);
	inFile.read((char *)&density,sizeof(tw::Float));
}

void UniformProfile::WriteData(std::ofstream& outFile)
{
	Profile::WriteData(outFile);
	outFile.write((char *)&density,sizeof(tw::Float));
}

tw::Float UniformProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	return theRgn->Inside(pos,ds) ? density : 0.0;
}

void GaussianProfile::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
    std::string word;
	Profile::ReadInputFileDirective(inputString,com);
	if (com=="density")
		inputString >> word >> density;
	if (com=="size")
		inputString >> word >> beamSize.x >> beamSize.y >> beamSize.z;
}

void GaussianProfile::ReadData(std::ifstream& inFile)
{
	Profile::ReadData(inFile);
	inFile.read((char *)&density,sizeof(tw::Float));
	inFile.read((char *)&beamSize,sizeof(tw::vec3));
}

void GaussianProfile::WriteData(std::ofstream& outFile)
{
	Profile::WriteData(outFile);
	outFile.write((char *)&density,sizeof(tw::Float));
	outFile.write((char *)&beamSize,sizeof(tw::vec3));
}

tw::Float GaussianProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	tw::Float dens = density;
	tw::vec3 p = pos - centerPt;
	orientation.ExpressInBasis(&p);
	dens *= exp(-sqr(p.x/beamSize.x));
	dens *= exp(-sqr(p.y/beamSize.y));
	dens *= exp(-sqr(p.z/beamSize.z));
	return theRgn->Inside(pos,ds) ? dens : 0.0;
}

void ChannelProfile::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
    std::string word;
	Profile::ReadInputFileDirective(inputString,com);
	if (com=="radius")
		throw tw::FatalError("Channel radius no longer supported.  Use polynomial coefficients.");
	if (com=="coefficients")
		inputString >> word >> n0 >> n2 >> n4 >> n6;
	if (com=="zpoints") // eg, zpoints = { 0 1 2 3 }
		tw::input::ReadArray(z,inputString);
	if (com=="zdensity") // eg, zdensity = { 0 1 1 0 }
		tw::input::ReadArray(fz,inputString);
}

void ChannelProfile::ReadData(std::ifstream& inFile)
{
	Profile::ReadData(inFile);
	inFile.read((char *)&n0,sizeof(tw::Float));
	inFile.read((char *)&n2,sizeof(tw::Float));
	inFile.read((char *)&n4,sizeof(tw::Float));
	inFile.read((char *)&n6,sizeof(tw::Float));

	tw::Int zDim;
	inFile.read((char *)&zDim,sizeof(tw::Int));
	z.resize(zDim);
	fz.resize(zDim);
	inFile.read((char *)&z[0],sizeof(tw::Float)*zDim);
	inFile.read((char *)&fz[0],sizeof(tw::Float)*zDim);
}

void ChannelProfile::WriteData(std::ofstream& outFile)
{
	Profile::WriteData(outFile);
	outFile.write((char *)&n0,sizeof(tw::Float));
	outFile.write((char *)&n2,sizeof(tw::Float));
	outFile.write((char *)&n4,sizeof(tw::Float));
	outFile.write((char *)&n6,sizeof(tw::Float));

	tw::Int zDim = z.size();
	outFile.write((char *)&zDim,sizeof(tw::Int));
	outFile.write((char *)&z[0],sizeof(tw::Float)*zDim);
	outFile.write((char *)&fz[0],sizeof(tw::Float)*zDim);
}

tw::Float ChannelProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	tw::Int i;
	tw::Float r2,w,dens = 0.0;
	tw::vec3 p = pos - centerPt;
	orientation.ExpressInBasis(&p);
	for (i=0;i<z.size()-1;i++)
	{
		if (p.z>=z[i] && p.z<=z[i+1])
		{
			w = (z[i+1]-p.z)/(z[i+1]-z[i]);
			if (segmentShape==loading::quartic)
				w = QuarticPulse(0.5*w);
			if (segmentShape==loading::quintic)
				w = QuinticRise(w);
			dens = fz[i]*w + fz[i+1]*(1.0-w);
		}
	}
	r2 = sqr(p.x) + sqr(p.y);
	dens *= n0 + n2*r2 + n4*r2*r2 + n6*r2*r2*r2;
	return theRgn->Inside(pos,ds) ? dens : 0.0;
}

void ColumnProfile::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
    std::string word;
	Profile::ReadInputFileDirective(inputString,com);
	if (com=="size")
		inputString >> word >> beamSize.x >> beamSize.y >> beamSize.z;
	if (com=="zpoints") // eg, zpoints = { 0 1 2 3 }
		tw::input::ReadArray(z,inputString);
	if (com=="zdensity") // eg, zdensity = { 0 1 1 0 }
		tw::input::ReadArray(fz,inputString);
}

void ColumnProfile::ReadData(std::ifstream& inFile)
{
	Profile::ReadData(inFile);
	inFile.read((char *)&beamSize,sizeof(tw::vec3));

	tw::Int zDim;
	inFile.read((char *)&zDim,sizeof(tw::Int));
	z.resize(zDim);
	fz.resize(zDim);
	inFile.read((char *)&z[0],sizeof(tw::Float)*zDim);
	inFile.read((char *)&fz[0],sizeof(tw::Float)*zDim);
}

void ColumnProfile::WriteData(std::ofstream& outFile)
{
	Profile::WriteData(outFile);
	outFile.write((char *)&beamSize,sizeof(tw::vec3));

	tw::Int zDim = z.size();
	outFile.write((char *)&zDim,sizeof(tw::Int));
	outFile.write((char *)&z[0],sizeof(tw::Float)*zDim);
	outFile.write((char *)&fz[0],sizeof(tw::Float)*zDim);
}

tw::Float ColumnProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	tw::Int i;
	tw::Float w,dens = 0.0;
	tw::vec3 p = pos - centerPt;
	orientation.ExpressInBasis(&p);
	for (i=0;i<z.size()-1;i++)
	{
		if (p.z>=z[i] && p.z<=z[i+1])
		{
			w = (z[i+1]-p.z)/(z[i+1]-z[i]);
			if (segmentShape==loading::quartic)
				w = QuarticPulse(0.5*w);
			if (segmentShape==loading::quintic)
				w = QuinticRise(w);
			dens = fz[i]*w + fz[i+1]*(1.0-w);
		}
	}
	dens *= exp(-sqr(p.x/beamSize.x));
	dens *= exp(-sqr(p.y/beamSize.y));
	return theRgn->Inside(pos,ds) ? dens : 0.0;
}

void PiecewiseProfile::Initialize(MetricSpace *ds)
{
	Profile::Initialize(ds);
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

void PiecewiseProfile::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
	Profile::ReadInputFileDirective(inputString,com);
	if (com=="xpoints") // eg, xpoints = { 0 1 2 3 }
		tw::input::ReadArray(x,inputString);
	if (com=="ypoints") // eg, ypoints = { 0 1 2 3 }
		tw::input::ReadArray(y,inputString);
	if (com=="zpoints") // eg, zpoints = { 0 1 2 3 }
		tw::input::ReadArray(z,inputString);
	if (com=="xdensity") // eg, xdensity = { 0 1 1 0 }
		tw::input::ReadArray(fx,inputString);
	if (com=="ydensity") // eg, ydensity = { 0 1 1 0 }
		tw::input::ReadArray(fy,inputString);
	if (com=="zdensity") // eg, zdensity = { 0 1 1 0 }
		tw::input::ReadArray(fz,inputString);
}

void PiecewiseProfile::ReadData(std::ifstream& inFile)
{
	Profile::ReadData(inFile);

	tw::Int xDim,yDim,zDim;
	inFile.read((char *)&xDim,sizeof(tw::Int));
	inFile.read((char *)&yDim,sizeof(tw::Int));
	inFile.read((char *)&zDim,sizeof(tw::Int));

	x.resize(xDim);
	fx.resize(xDim);
	inFile.read((char *)&x[0],sizeof(tw::Float)*xDim);
	inFile.read((char *)&fx[0],sizeof(tw::Float)*xDim);

	y.resize(yDim);
	fy.resize(yDim);
	inFile.read((char *)&y[0],sizeof(tw::Float)*yDim);
	inFile.read((char *)&fy[0],sizeof(tw::Float)*yDim);

	z.resize(zDim);
	fz.resize(zDim);
	inFile.read((char *)&z[0],sizeof(tw::Float)*zDim);
	inFile.read((char *)&fz[0],sizeof(tw::Float)*zDim);
}

void PiecewiseProfile::WriteData(std::ofstream& outFile)
{
	Profile::WriteData(outFile);

	tw::Int xDim,yDim,zDim;
	xDim = x.size();
	yDim = y.size();
	zDim = z.size();

	outFile.write((char *)&xDim,sizeof(tw::Int));
	outFile.write((char *)&yDim,sizeof(tw::Int));
	outFile.write((char *)&zDim,sizeof(tw::Int));

	outFile.write((char *)&x[0],sizeof(tw::Float)*xDim);
	outFile.write((char *)&fx[0],sizeof(tw::Float)*xDim);

	outFile.write((char *)&y[0],sizeof(tw::Float)*yDim);
	outFile.write((char *)&fy[0],sizeof(tw::Float)*yDim);

	outFile.write((char *)&z[0],sizeof(tw::Float)*zDim);
	outFile.write((char *)&fz[0],sizeof(tw::Float)*zDim);
}

tw::Float PiecewiseProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	tw::Float ansX,ansY,ansZ;
	tw::Int i;
	tw::Float r,w;

	tw::vec3 p = pos - centerPt;
	orientation.ExpressInBasis(&p);

	const tw::Float x0 = p.x;
	const tw::Float y0 = p.y;
	const tw::Float z0 = p.z;

	tw::Int xDim = x.size();
	tw::Int yDim = y.size();
	tw::Int zDim = z.size();

	ansX = ansY = ansZ = 1.0;

	switch (symmetry)
	{
		case cartesian:
			ansX = 0.0;
			for (i=0;i<xDim-1;i++)
			{
				if (x0>=x[i] && x0<=x[i+1])
				{
					w = (x[i+1]-x0)/(x[i+1]-x[i]);
					if (segmentShape==loading::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==loading::quintic)
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
					if (segmentShape==loading::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==loading::quintic)
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
					if (segmentShape==loading::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==loading::quintic)
						w = QuinticRise(w);
					ansZ = fz[i]*w + fz[i+1]*(1.0-w);
				}
			}
			break;
		case cylindrical:
			r = sqrt(sqr(x0 - x[0]) + sqr(y0 - 0.5*(y[0] + y[yDim-1])));
			r += x[0];
			ansX = 0.0;
			for (i=0;i<xDim-1;i++)
			{
				if (r>=x[i] && r<=x[i+1])
				{
					w = (x[i+1]-r)/(x[i+1]-x[i]);
					if (segmentShape==loading::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==loading::quintic)
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
					if (segmentShape==loading::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==loading::quintic)
						w = QuinticRise(w);
					ansZ = fz[i]*w + fz[i+1]*(1.0-w);
				}
			}

			ansY = 1.0;
			break;
		case spherical:
			r = sqrt(sqr(x0 - x[0]) + sqr(y0 - 0.5*(y[0] + y[yDim-1])) + sqr(z0 - 0.5*(z[0] + z[zDim-1])));
			r += x[0];
			ansX = 0.0;
			for (i=0;i<xDim-1;i++)
			{
				if (r>=x[i] && r<=x[i+1])
				{
					w = (x[i+1]-r)/(x[i+1]-x[i]);
					if (segmentShape==loading::quartic)
						w = QuarticPulse(0.5*w);
					if (segmentShape==loading::quintic)
						w = QuinticRise(w);
					ansX = fx[i]*w + fx[i+1]*(1.0-w);
				}
			}

			ansY = 1.0;
			ansZ = 1.0;
			break;
	}

	tw::Float dens = ansX*ansY*ansZ*sqr(cos(0.5*modeNumber.x*p.x)*cos(0.5*modeNumber.y*p.y)*cos(0.5*modeNumber.z*p.z));
	return theRgn->Inside(pos,ds) ? dens : 0.0;
}

void CorrugatedProfile::ReadInputFileDirective(std::stringstream& inputString,const std::string& com)
{
    std::string word;
	Profile::ReadInputFileDirective(inputString,com);
	if (com=="a0")
		inputString >> word >> a0;
	if (com=="gamma0")
		inputString >> word >> gamma0;
	if (com=="w0")
		inputString >> word >> w0;
	if (com=="wp0")
		inputString >> word >> wp0;
	if (com=="km")
		inputString >> word >> km;
	if (com=="delta")
		inputString >> word >> delta;
	if (com=="rchannel")
		inputString >> word >> rchannel;
}

void CorrugatedProfile::ReadData(std::ifstream& inFile)
{
	Profile::ReadData(inFile);
	inFile.read((char *)&a0,sizeof(tw::Float));
	inFile.read((char *)&gamma0,sizeof(tw::Float));
	inFile.read((char *)&w0,sizeof(tw::Float));
	inFile.read((char *)&wp0,sizeof(tw::Float));
	inFile.read((char *)&km,sizeof(tw::Float));
	inFile.read((char *)&delta,sizeof(tw::Float));
	inFile.read((char *)&rchannel,sizeof(tw::Float));
}

void CorrugatedProfile::WriteData(std::ofstream& outFile)
{
	Profile::WriteData(outFile);
	outFile.write((char *)&a0,sizeof(tw::Float));
	outFile.write((char *)&gamma0,sizeof(tw::Float));
	outFile.write((char *)&w0,sizeof(tw::Float));
	outFile.write((char *)&wp0,sizeof(tw::Float));
	outFile.write((char *)&km,sizeof(tw::Float));
	outFile.write((char *)&delta,sizeof(tw::Float));
	outFile.write((char *)&rchannel,sizeof(tw::Float));
}

tw::Float CorrugatedProfile::GetValue(const tw::vec3& pos,const MetricSpace& ds)
{
	tw::Float dens,wp1s,r2,psi,a0Hat,kHat;
	tw::vec3 p = pos - centerPt;
	orientation.ExpressInBasis(&p);
	const tw::Float x = p.x;
	const tw::Float y = p.y;
	const tw::Float z = p.z; // z = 0 is initial injection point
	r2 = x*x + y*y;
	psi = delta*wp0*wp0/(2.0*w0*km);
	a0Hat = 4.0*j1(psi)*a0/(w0*rchannel);
	kHat = w0 + km - 0.5*w0*(sqr(wp0/w0) + 8.0/sqr(w0*rchannel));
	wp1s = 2.0*w0*kHat - 2.0*sqrt(sqr(gamma0+a0Hat*w0*z)/(sqr(gamma0+a0Hat*w0*z)-1.0))*w0*w0;
	dens = wp0*wp0*(1.0 + delta*sin(km*z)) + wp1s + 4.0*r2/pow(rchannel,tw::Float(4.0));
	return theRgn->Inside(pos,ds) ? dens : 0.0;
}



///////////////////////////
//                       //
// PULSE and BEAM SHAPES //
//                       //
///////////////////////////

PulseShape::PulseShape()
{
	whichProfile = loading::quintic;
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
		case loading::triangle:
			return tau;
		case loading::sech:
			return SechRise(tau);
		case loading::quartic:
			return QuarticRise(tau);
		case loading::quintic:
			return QuinticRise(tau);
		case loading::sin2:
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
		case loading::triangle:
			return 0.0;
		case loading::sech:
			return w*D1SechRise(tau);
		case loading::quartic:
			return w*D1QuarticRise(tau);
		case loading::quintic:
			return w*D1QuinticRise(tau);
		case loading::sin2:
      return w*D1Sin2Rise(tau);
	}
	return 0.0;
}

tw::Float PulseShape::D1Intensity(const tw::Float t) const
{
    return 2.0*PulseShapeFactor(t)*D1Amplitude(t);
}


////////////////////
//                //
// EXPLICIT WAVES //
//                //
////////////////////


Wave::Wave(GaussianDeviate *gd)
{
	direction = tw::vec3(0.0,0.0,1.0);
	focusPosition = tw::vec3(0.0,0.0,0.0);
	a = tw::vec3(0.1,0.0,0.0);
	nrefr = 1.0;
	w = 1.0;
	chirp = 0.0;
	phase = 0.0;
	randomPhase = 0.0;
	modeType = EM::hermite;
	modeData[0].order = 0;
	modeData[1].order = 0;
	modeData[0].scale = 1.0;
	modeData[1].scale = 1.0;
	modeData[0].exponent = 2;
	modeData[1].exponent = 2;
	deviate = gd;
}

void Wave::Initialize()
{
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

tw::Complex Wave::PlanePrimitive(tw::Float t,const tw::vec3& r) const
{
	// t,r are expected to be in the Cartesian laser basis for all primitive functions
	const tw::Float tau = t - nrefr*r.z;
	const tw::Float psi = phase - w*tau - chirp*tau*tau;
	return a0 * pulseShape.PulseShapeFactor(tau) * std::exp(ii*psi);
}

tw::Complex Wave::BesselPrimitive(tw::Float t,const tw::vec3& r) const
{
	const tw::Float rbar = sqrt(sqr(r.x) + sqr(r.y))/modeData[0].scale;
	const tw::Float tau = t - nrefr*r.z;
	const tw::Float psi = phase - w*tau - chirp*tau*tau;
	return a0 * j0(rbar) * pulseShape.PulseShapeFactor(tau) * std::exp(ii*psi);
}

tw::Complex Wave::AiryDiscPrimitive(tw::Float t,const tw::vec3& r) const
{
	const tw::Float rbar = sqrt(sqr(r.x) + sqr(r.y))/modeData[0].scale;
	const tw::Float tau = t - nrefr*r.z;
	const tw::Float psi = phase - w*tau - chirp*tau*tau;
	return 2*a0 * (j1(rbar)/rbar) * pulseShape.PulseShapeFactor(tau) * std::exp(ii*psi);
}

tw::Complex Wave::LaguerrePrimitive(tw::Float t,const tw::vec3& r) const
{
	const tw::Float r0 = modeData[0].scale;
	const tw::Int nexp = modeData[0].exponent;
	const tw::Int mu = modeData[0].order;
	const tw::Int mv = modeData[1].order;

	tw::Float Ax=a0, tau=t-nrefr*r.z, guoy_shift=0.0;
	tw::Float zR,rm,rho,phi;

	// setup for Laguerre-Gaussian mode
	rho = sqrt(r.x*r.x + r.y*r.y);
	phi = atan2(r.y,r.x);
	zR = 0.5*nrefr*w*r0*r0;
	rm = r0*sqrt(one + r.z*r.z/(zR*zR));
	guoy_shift = -(one + two*mu + mv)*atan(r.z/zR);

	// compute factors for r-phi profile
	Ax *= r0/rm;
	Ax *= Laguerre(mu,mv,two*sqr(rho/rm));
	Ax *= nexp%2==0 ? exp(-pow(rho/rm,nexp)) : pow(cos(0.5*pi*rho/rm),nexp+1)*tw::Float(rho<rm);
	Ax *= pow(sqrt(two)*rho/rm,tw::Float(mv)); // put appropriate hole on axis
	Ax *= sqrt(Factorial(mu)/Factorial(mu+mv)); // Laguerre normalization
	tau += guoy_shift/w - 0.5*nrefr*rho*rho*r.z/(r.z*r.z + zR*zR);
	tw::Float psi = phase + 2.0*guoy_shift - mv*phi - w*tau - chirp*tau*tau;

	return Ax * pulseShape.PulseShapeFactor(tau) * std::exp(ii*psi);
}

tw::Complex Wave::HermitePrimitive(tw::Float t,const tw::vec3& r) const
{
	const tw::Float r0[2] = { modeData[0].scale , modeData[1].scale };
	const tw::Int nexp[2] = { modeData[0].exponent , modeData[1].exponent };
	const tw::Int order[2] = { modeData[0].order , modeData[1].order };

	tw::Float Ax=a0, tau=t-nrefr*r.z, guoy_shift=0.0;

	// lambda function to compute Hermite amplitude and phase per axis
	auto herm = [&] (tw::Float x,tw::Int ax)
	{
		const tw::Float m = order[ax-1];
		const tw::Float r00 = r0[ax-1];
		const tw::Float zR = 0.5*nrefr*w*r00*r00;
		const tw::Float rm = r00*sqrt(1.0 + r.z*r.z/(zR*zR));
		guoy_shift -= (0.5 + m)*atan(r.z/zR);
		Ax *= sqrt(r00/rm);
		Ax *= Hermite(m,1.414*x/rm);
		Ax *= nexp[ax-1]%2==0 ? exp(-pow(x/rm,nexp[ax-1])) : pow(cos(0.5*pi*x/rm),nexp[ax-1]+1)*tw::Float(x*x<rm*rm);
		Ax /= sqrt(pow(two,tw::Float(m))*Factorial(m)); // Hermite normalization
		tau += guoy_shift/w - 0.5*nrefr*x*x*r.z/(r.z*r.z + zR*zR);
	};

	herm(r.x,1);
	herm(r.y,2);
	tw::Float psi = phase + 2.0*guoy_shift - w*tau - chirp*tau*tau;

	return Ax * pulseShape.PulseShapeFactor(tau) * std::exp(ii*psi);
}

tw::Complex Wave::MultipolePrimitive(tw::Float t, const tw::vec3& r) const
{
	// For now hard code in magnetic dipole radiation
	const tw::Float j1max = 0.436182;
	const tw::Float rho = sqrt(r.x*r.x + r.y*r.y);
	const tw::Float R = sqrt(rho*rho + r.z*r.z);
	const tw::Float stheta = rho/R;
	const tw::Float wR = w*R;
	//auto spherical_bessel = [&] (tw::Float x) { return sin(x)/sqr(x+tw::small_pos) - cos(x)/(x+tw::small_pos); };

	tw::Complex Aphi;
	tw::Float tau_in,tau_out;
	tau_out = t - nrefr*R;
	tau_in = t + nrefr*R;
	Aphi = (-one/wR + one/(ii*wR*wR)) * pulseShape.PulseShapeFactor(tau_out) * std::exp(-ii*w*tau_out);
	Aphi += (-one/wR - one/(ii*wR*wR)) * pulseShape.PulseShapeFactor(tau_in) * std::exp(-ii*w*tau_in);
	Aphi *= 0.5 * a0 * (stheta/j1max);
	return Aphi;
}

void Wave::ReadInputFile(std::stringstream& inputString,std::string& command)
{
	std::string word;

	do
	{
		inputString >> word;

		if (word=="use")
			throw tw::FatalError("'use' directives not supported.  See docs for alternative.");

		if (word=="direction") // eg, direction = (1.0,1.0,0.0)
		{
			inputString >> word;
			inputString >> direction.x >> direction.y >> direction.z;
		}
		if (word=="focus") // eg, focus position = (0.0,0.0,0.0)
		{
			inputString >> word >> word;
			inputString >> focusPosition.x >> focusPosition.y >> focusPosition.z;
		}
		if (word=="a") // eg, a = 0.0 0.0 0.1
		{
			inputString >> word;
			inputString >> a.x >> a.y >> a.z;
		}
		if (word=="w") // eg, w = 2.0
		{
			inputString >> word;
			inputString >> w;
		}
		if (word=="refractiveindex") // eg, refractiveindex = 1.5
		{
			inputString >> word;
			inputString >> nrefr;
		}
		if (word=="chirp") // eg, chirp = 0.01
		{
			inputString >> word;
			inputString >> chirp;
		}
		if (word=="phase") // eg, phase = 90
		{
			inputString >> word;
			inputString >> phase;
			phase *= pi/180.0;
		}
		if (word=="random") // eg, random phase = 10
		{
			inputString >> word >> word;
			inputString >> randomPhase;
			randomPhase *= pi/180.0;
		}
		if (word=="delay") // eg, delay = 10
		{
			inputString >> word;
			inputString >> pulseShape.delay;
		}
		if (word=="risetime") // eg, risetime = 1.0
		{
			inputString >> word;
			inputString >> pulseShape.risetime;
			if (pulseShape.risetime<=0.0)
				throw tw::FatalError("Pulse rise time must be positive and non-zero.");
		}
		if (word=="holdtime") // eg, holdtime = 10
		{
			inputString >> word;
			inputString >> pulseShape.holdtime;
			if (pulseShape.holdtime<0.0)
				throw tw::FatalError("Pulse hold time must be positive.");
		}
		if (word=="falltime") // eg, falltime = 10
		{
			inputString >> word;
			inputString >> pulseShape.falltime;
			if (pulseShape.falltime<=0.0)
				throw tw::FatalError("Pulse fall time must be positive and non-zero.");
		}
		if (word=="r0") // eg, r0 = 1.0 1.0
		{
			inputString >> word >> modeData[0].scale >> modeData[1].scale;
		}
		if (word=="exponent") // eg, exponent = 2 2
		{
			inputString >> word >> modeData[0].exponent >> modeData[1].exponent;
		}
		if (word=="mode")
		{
			inputString >> word >> modeData[0].order >> modeData[1].order;
		}
		if (word=="type")
		{
			inputString >> word >> word;
			if (word=="plane")
				modeType = EM::plane;
			if (word=="multipole")
				modeType = EM::multipole;
			if (word=="hermite")
				modeType = EM::hermite;
			if (word=="laguerre")
				modeType = EM::laguerre;
			if (word=="bessel")
				modeType = EM::bessel;
			if (word=="airy_disc")
				modeType = EM::airy_disc;
		}
		if (word=="shape")
		{
			inputString >> word >> word;
			if (word=="quintic")
				pulseShape.whichProfile = loading::quintic;
			if (word=="sin2")
				pulseShape.whichProfile = loading::sin2;
			if (word=="sech")
				pulseShape.whichProfile = loading::sech;
		}
	} while (word!="}");

	if (w==0.0 && (modeType==EM::hermite || modeType==EM::laguerre || modeType==EM::multipole))
		throw tw::FatalError("Zero frequency requested for an EM mode that does not support it.");
}

void Wave::ReadData(std::ifstream& inFile)
{
	inFile.read((char *)&direction,sizeof(tw::vec3));
	inFile.read((char *)&focusPosition,sizeof(tw::vec3));
	inFile.read((char *)&a,sizeof(tw::vec3));
	inFile.read((char *)&a0,sizeof(tw::Float));
	inFile.read((char *)&w,sizeof(tw::Float));
	inFile.read((char *)&nrefr,sizeof(tw::Float));
	inFile.read((char *)&phase,sizeof(tw::Float));
	inFile.read((char *)&vg,sizeof(tw::Float));
	inFile.read((char *)&chirp,sizeof(tw::Float));
	inFile.read((char *)&randomPhase,sizeof(tw::Float));
	inFile.read((char *)&pulseShape,sizeof(PulseShape));
	inFile.read((char *)&modeType,sizeof(modeType));
	inFile.read((char *)modeData,sizeof(modeData));
	inFile.read((char *)&laserFrame,sizeof(laserFrame));
}

void Wave::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&direction,sizeof(tw::vec3));
	outFile.write((char *)&focusPosition,sizeof(tw::vec3));
	outFile.write((char *)&a,sizeof(tw::vec3));
	outFile.write((char *)&a0,sizeof(tw::Float));
	outFile.write((char *)&w,sizeof(tw::Float));
	outFile.write((char *)&nrefr,sizeof(tw::Float));
	outFile.write((char *)&phase,sizeof(tw::Float));
	outFile.write((char *)&vg,sizeof(tw::Float));
	outFile.write((char *)&chirp,sizeof(tw::Float));
	outFile.write((char *)&randomPhase,sizeof(tw::Float));
	outFile.write((char *)&pulseShape,sizeof(PulseShape));
	outFile.write((char *)&modeType,sizeof(modeType));
	outFile.write((char *)modeData,sizeof(modeData));
	outFile.write((char *)&laserFrame,sizeof(laserFrame));
}


////////////////////
//                //
//   CONDUCTOR    //
//                //
////////////////////


Conductor::Conductor(std::vector<Region*>& rlist) : rgnList(rlist)
{
	affectsPhi = true;
	affectsA = true;
	magneticCurrent = false;
	electricCurrent = false;
	pulseShape.delay = 0.0;
	pulseShape.risetime = tw::small_pos;
	pulseShape.holdtime = tw::big_pos;
	pulseShape.falltime = tw::small_pos;
	gaussianRadius = tw::big_pos;
	f = tw::big_pos;
	ks = 0.0;
	theRgn = rgnList[0];
}

void Conductor::Initialize(const MetricSpace& ds)
{
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
	tw::Int i;
	tw::Float ans = 0.0;
	for (i=0;i<potential.size();i++)
		ans += pulseShape.PulseShapeFactor(t)*potential[i]*cos(angFreq[i]*t + phase[i]);
	return ans;
}

tw::Float Conductor::VoltageRate(tw::Float t)
{
	tw::Int i;
	tw::Float ans = 0.0;
	for (i=0;i<potential.size();i++)
	{
		ans -= pulseShape.PulseShapeFactor(t)*angFreq[i]*potential[i]*sin(angFreq[i]*t + phase[i]);
		ans += pulseShape.D1Amplitude(t)*potential[i]*cos(angFreq[i]*t + phase[i]);
	}
	return ans;
}

tw::vec3 Conductor::PolarizationDensity(const tw::vec3& pos,tw::Float t)
{
	tw::Int i;
	tw::vec3 P0;
	tw::vec3 ans(0.0,0.0,0.0);
	tw::vec3 rc = pos - theRgn->center;
	theRgn->orientation.ExpressInBasis(&rc);
	for (i=0;i<potential.size();i++)
	{
		P0 = tw::vec3(Px[i],Py[i],Pz[i]);
		ans += P0*sin(angFreq[i]*t + phase[i] + (0.5*angFreq[i]/f)*(rc.x*rc.x + rc.y*rc.y) + (ks^rc));
	}
	ans *= pulseShape.PulseShapeFactor(t + (0.5/f)*(rc.x*rc.x + rc.y*rc.y));
	ans *= exp(-sqr(rc.x/gaussianRadius.x)-sqr(rc.y/gaussianRadius.y)-sqr(rc.z/gaussianRadius.z));
	theRgn->orientation.ExpressInStdBasis(&ans);
	return ans;
}

void Conductor::DepositSources(Field& sources,const MetricSpace& m,tw::Float t,tw::Float dt)
{
	tw::Int x0,x1,y0,y1,z0,z1;
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

void Conductor::ReadInputFile(std::stringstream& inputString,std::string& command)
{
	std::string word;

	do
	{
		inputString >> word;

		if (word=="use")
			throw tw::FatalError("'use' directives not supported.  See docs for alternative.");

		theRgn = Region::ReadRegion(rgnList,theRgn,inputString,word);

		if (word=="px") // eg, px = { 1 , 2 , 3 }
			tw::input::ReadArray(Px,inputString);
		if (word=="py") // eg, py = { 1 , 2 , 3 }
			tw::input::ReadArray(Py,inputString);
		if (word=="pz") // eg, pz = { 1 , 2 , 3 }
			tw::input::ReadArray(Pz,inputString);
		if (word=="potential") // eg, potential = { 1 , 2 , 3 }
			tw::input::ReadArray(potential,inputString);
		if (word=="w") // eg, w = { 0 , 1 , 2 }
			tw::input::ReadArray(angFreq,inputString);
		if (word=="phase") // eg, phase = { 0 , 0 , 90 }
		{
			tw::input::ReadArray(phase,inputString);
			phase *= pi/180.0;
		}
		if (word=="delay") // eg, delay = 10
		{
			inputString >> word;
			inputString >> pulseShape.delay;
		}
		if (word=="risetime") // eg, risetime = 1.0
		{
			inputString >> word;
			inputString >> pulseShape.risetime;
			if (pulseShape.risetime<=0.0)
				throw tw::FatalError("Pulse rise time must be positive and non-zero.");
		}
		if (word=="holdtime") // eg, holdtime = 10
		{
			inputString >> word;
			inputString >> pulseShape.holdtime;
			if (pulseShape.holdtime<0.0)
				throw tw::FatalError("Pulse hold time must be positive.");
		}
		if (word=="falltime") // eg, falltime = 10
		{
			inputString >> word;
			inputString >> pulseShape.falltime;
			if (pulseShape.falltime<=0.0)
				throw tw::FatalError("Pulse fall time must be positive and non-zero.");
		}
		if (word=="shape")
		{
			inputString >> word >> word;
			if (word=="quintic")
				pulseShape.whichProfile = loading::quintic;
			if (word=="sin2")
				pulseShape.whichProfile = loading::sin2;
			if (word=="sech")
				pulseShape.whichProfile = loading::sech;
		}
		if (word=="enable")
		{
			inputString >> word;
			if (word=="electrostatic") // eg, enable electrostatic = true
			{
				inputString >> word >> word;
				affectsPhi = (word=="yes" || word=="true" || word=="on");
			}
			if (word=="electromagnetic") // eg, enable electromagnetic = true
			{
				inputString >> word >> word;
				affectsA = (word=="yes" || word=="true" || word=="on");
			}
		}
		if (word=="current") // eg, current type = electric
		{
			inputString >> word >> word >> word;
			if (word=="none")
				electricCurrent = magneticCurrent = false;
			if (word=="magnetic")
				throw tw::FatalError("Magnetic currents not supported this version.");
			if (word=="electric")
				electricCurrent = true;
		}
		if (word=="gaussian") // eg, gaussian size = 1 1 1
		{
			inputString >> word >> word >> gaussianRadius.x >> gaussianRadius.y >> gaussianRadius.z;
		}
		if (word=="f") // eg, f = 100.0
		{
			inputString >> word >> f;
		}
		if (word=="ks") // eg, ks = 1.0 0.0 0.0
		{
			inputString >> word >> ks.x >> ks.y >> ks.z;
		}
	} while (word!="}");
}

void Conductor::ReadData(std::ifstream& inFile)
{
	tw::Int num,rgnIndex;
    inFile.read((char *)&rgnIndex,sizeof(tw::Int));
    theRgn = rgnList[rgnIndex];
	inFile.read((char *)&pulseShape,sizeof(PulseShape));
	inFile.read((char *)&num,sizeof(tw::Int));
	Px.resize(num);
	Py.resize(num);
	Pz.resize(num);
	potential.resize(num);
	angFreq.resize(num);
	phase.resize(num);
	inFile.read((char *)&Px[0],sizeof(tw::Float)*num);
	inFile.read((char *)&Py[0],sizeof(tw::Float)*num);
	inFile.read((char *)&Pz[0],sizeof(tw::Float)*num);
	inFile.read((char *)&potential[0],sizeof(tw::Float)*num);
	inFile.read((char *)&angFreq[0],sizeof(tw::Float)*num);
	inFile.read((char *)&phase[0],sizeof(tw::Float)*num);
	inFile.read((char *)&affectsPhi,sizeof(bool));
	inFile.read((char *)&affectsA,sizeof(bool));
	inFile.read((char *)&electricCurrent,sizeof(bool));
	inFile.read((char *)&magneticCurrent,sizeof(bool));
	inFile.read((char *)&gaussianRadius,sizeof(tw::vec3));
	inFile.read((char *)&f,sizeof(tw::Float));
	inFile.read((char *)&ks,sizeof(tw::vec3));
}

void Conductor::WriteData(std::ofstream& outFile)
{
    tw::Int rgnIndex = std::find(rgnList.begin(),rgnList.end(),theRgn) - rgnList.begin();
    outFile.write((char *)&rgnIndex,sizeof(tw::Int));
	tw::Int num = Px.size();
	outFile.write((char *)&pulseShape,sizeof(PulseShape));
	outFile.write((char *)&num,sizeof(tw::Int));
	outFile.write((char *)&Px[0],sizeof(tw::Float)*num);
	outFile.write((char *)&Py[0],sizeof(tw::Float)*num);
	outFile.write((char *)&Pz[0],sizeof(tw::Float)*num);
	outFile.write((char *)&potential[0],sizeof(tw::Float)*num);
	outFile.write((char *)&angFreq[0],sizeof(tw::Float)*num);
	outFile.write((char *)&phase[0],sizeof(tw::Float)*num);
	outFile.write((char *)&affectsPhi,sizeof(bool));
	outFile.write((char *)&affectsA,sizeof(bool));
	outFile.write((char *)&electricCurrent,sizeof(bool));
	outFile.write((char *)&magneticCurrent,sizeof(bool));
	outFile.write((char *)&gaussianRadius,sizeof(tw::vec3));
	outFile.write((char *)&f,sizeof(tw::Float));
	outFile.write((char *)&ks,sizeof(tw::vec3));
}

//////////////////////
//                  //
// LINDMAN BOUNDARY //
//                  //
//////////////////////


void LindmanBoundary::Initialize(Task *task,MetricSpace *ms,std::vector<Wave*> *wave,const axisSpec& axis,const sideSpec& side,tw::Int c0,tw::Int c1)
{
	this->ms = ms;
	this->wave = wave;
	this->axis = axis;
	this->side = side;
	// Set() assumes that c0 is an ordered 4-vector index
	this->c0 = c0;
	this->c1 = c1;

	DiscreteSpace bm_layout;

	if (axis==xAxis && ms->Dim(1)==1)
		throw tw::FatalError("Lindman boundary geometry is not consistent");

	if (axis==yAxis && ms->Dim(2)==1)
		throw tw::FatalError("Lindman boundary geometry is not consistent");

	if (axis==zAxis && ms->Dim(3)==1)
		throw tw::FatalError("Lindman boundary geometry is not consistent");

	switch (axis)
	{
		case tAxis:
			break;
		case xAxis:
			bm_layout.Resize(1,ms->Dim(2),ms->Dim(3),Corner(*ms),PhysicalSize(*ms));
			break;
		case yAxis:
			bm_layout.Resize(ms->Dim(1),1,ms->Dim(3),Corner(*ms),PhysicalSize(*ms));
			break;
		case zAxis:
			bm_layout.Resize(ms->Dim(1),ms->Dim(2),1,Corner(*ms),PhysicalSize(*ms));
			break;
	}

	boundaryMemory.Initialize( 9 , bm_layout , task );
	boundaryMemory.SetBoundaryConditions(xAxis,none,none);
	boundaryMemory.SetBoundaryConditions(yAxis,none,none);
	boundaryMemory.SetBoundaryConditions(zAxis,none,none);
}

void LindmanBoundary::UpdateBoundaryMemory(Field& A,tw::Float dt)
{
	tw::Float alpha[3] = {0.3264,0.1272,0.0309};
	tw::Float beta[3] = {0.7375,0.98384,0.9996472};

	tw::Int ax = naxis(axis);
	tw::Float source,dtt;
	tw::Int i,j,k,s,s0,s1;

	s0 = side==lowSide ? 0 : A.Dim(ax);
	s1 = side==lowSide ? 1 : A.Dim(ax)+1;

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
	tw::Int ax = naxis(axis);
	tw::Float sgn,correction,ds[4] = {dt,dx(A),dy(A),dz(A)};
	tw::Int i,j,k,s,s0,s1;
	tw::vec3 pos,Ain;

	tw::Float propagator = (1.0 - ds[0]/ds[ax])/(1.0 + ds[0]/ds[ax]);
	tw::Float injector = 1.0/(1.0 + ds[0]/ds[ax]);

	tw::vec3 direction(tw::Int(ax==1),tw::Int(ax==2),tw::Int(ax==3));
	tw::vec3 offset = 0.5 * direction * ds[ax];
	sgn = side==highSide ? -1.0 : 1.0;

	s0 = side==lowSide ? 0 : A.Dim(ax)+1;
	s1 = side==lowSide ? 1 : A.Dim(ax);
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

void LindmanBoundary::ReadData(std::ifstream& inFile)
{
	inFile.read((char*)&c0,sizeof(tw::Int));
	inFile.read((char*)&c1,sizeof(tw::Int));
	boundaryMemory.ReadData(inFile);
}

void LindmanBoundary::WriteData(std::ofstream& outFile)
{
	outFile.write((char*)&c0,sizeof(tw::Int));
	outFile.write((char*)&c1,sizeof(tw::Int));
	boundaryMemory.WriteData(outFile);
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
