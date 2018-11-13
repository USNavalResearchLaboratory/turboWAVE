enum tw_profile_quantity { densityProfile, energyProfile, pxProfile, pyProfile, pzProfile };
enum tw_profile_spec { nullProfile, uniformProfile, channelProfile, gaussianProfile, columnProfile, piecewiseProfile, corrugatedProfile};
enum tw_profile_timing { triggeredProfile, maintainedProfile };
enum tw_load_method {statistical,deterministic};

namespace loading
{
	enum ramp_shape { triangle , sin2 , quartic , quintic , sech };
}

namespace EM
{
	enum mode { plane , hermite , laguerre , bessel , multipole };
	struct mode_data
	{
		tw::Int order,exponent;
		tw::Float scale;
	};
}

struct Profile
{
	tw_profile_spec profileSpec;
	loading::ramp_shape segmentShape;
	tw_geometry symmetry;
	tw::vec3 centerPt;
	tw::vec3 modeNumber;
	tw::Float modeAmplitude;
	tw::basis orientation;
	Region *theRgn;
	std::vector<Region*>& rgnList;

	// items needed for particle/fluid loading
	tw_profile_quantity whichQuantity;
	tw::vec3 thermalMomentum,driftMomentum;
	tw::Float temperature;
	bool neutralize,variableCharge;
	tw_load_method loadingMethod;
	tw_profile_timing timingMethod;
	tw::Float t0,t1;
	bool wasTriggered;

	Profile(std::vector<Region*>& rgnList);
	virtual void Initialize(MetricSpace *ds) {;}
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadInputFileBlock(std::stringstream& inputString,bool neutralize);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	static Profile* CreateObjectFromFile(std::vector<Region*>& rgnList,std::ifstream& inFile);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct UniformProfile:Profile
{
	tw::Float density;

	UniformProfile(std::vector<Region*>& rgnList):Profile(rgnList) { profileSpec = uniformProfile; }
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct GaussianProfile:Profile
{
	tw::Float density;
	tw::vec3 beamSize;

	GaussianProfile(std::vector<Region*>& rgnList):Profile(rgnList) { profileSpec = gaussianProfile; }
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct ChannelProfile:Profile
{
	std::valarray<tw::Float> z,fz;
	tw::Float n0,n2,n4,n6;

	ChannelProfile(std::vector<Region*>& rgnList):Profile(rgnList) { profileSpec = channelProfile; }
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct ColumnProfile:Profile
{
	std::valarray<tw::Float> z,fz;
	tw::vec3 beamSize;

	ColumnProfile(std::vector<Region*>& rgnList):Profile(rgnList) { profileSpec = columnProfile; }
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct PiecewiseProfile:Profile
{
	std::valarray<tw::Float> x,fx,y,fy,z,fz;

	PiecewiseProfile(std::vector<Region*>& rgnList):Profile(rgnList) { profileSpec = piecewiseProfile; }
	virtual void Initialize(MetricSpace *ds);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct CorrugatedProfile:Profile
{
	tw::Float a0,gamma0,w0,wp0,km,delta,rchannel;

	CorrugatedProfile(std::vector<Region*>& rgnList):Profile(rgnList) { profileSpec = corrugatedProfile; }
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct PulseShape
{
	tw::Float delay,risetime,holdtime,falltime;
	tw::Float t1,t2,t3,t4;
	loading::ramp_shape whichProfile;

	PulseShape();
	void Initialize(const tw::Float time_origin);
	tw::Float PulseShapeFactor(const tw::Float t) const;
	tw::Float FastPulseShapeFactor(const tw::Float t) const
	{
		const tw::Float hold = tw::Float(t > t2 && t <= t3);
		const tw::Float tau_rise = tw::Float(t > t1 && t <= t2) * (t-t1) / (t2-t1);
		const tw::Float tau_fall = tw::Float(t > t3 && t <= t4) * (1.0 - (t-t3) / (t4-t3));
		return QuinticRise(tau_rise + tau_fall + hold);
	}

	tw::Float D1Intensity(const tw::Float t) const;
	tw::Float D1Amplitude(const tw::Float t) const;
};

struct Wave
{
	// Eventually break this out into polymorphic ComputeTool
	tw::vec3 direction;
	tw::vec3 focusPosition;
	tw::vec3 a;
	tw::Float a0,w,nrefr,w0,k0,phase,vg,chirp,randomPhase;
	PulseShape pulseShape;
	EM::mode modeType;
	EM::mode_data modeData[2];
	tw::basis laserFrame;

	GaussianDeviate *deviate;

	Wave(GaussianDeviate *deviate);
	void Initialize();

	tw::Complex PlanePrimitive(tw::Float t,const tw::vec3& pos) const;
	tw::Complex BesselPrimitive(tw::Float t,const tw::vec3& pos) const;
	tw::Complex MultipolePrimitive(tw::Float t,const tw::vec3& pos) const;
	tw::Complex HermitePrimitive(tw::Float t,const tw::vec3& pos) const;
	tw::Complex LaguerrePrimitive(tw::Float t,const tw::vec3& pos) const;

	// Dispatch function for all mode types
	tw::vec3 VectorPotential(tw::Float time,const tw::vec3& pos) const
	{
		tw::Complex Ax;
		tw::vec3 ans;
		tw::vec3 r = pos - focusPosition;
		laserFrame.ExpressInBasis(&r);

		const tw::Float rho = tw::small_pos + sqrt(r.x*r.x+r.y*r.y);
		const tw::Float t = time - pulseShape.delay - pulseShape.risetime;
		switch (modeType)
		{
			case EM::plane:
				Ax = PlanePrimitive(t,r);
				break;
			case EM::bessel:
				Ax = BesselPrimitive(t,r);
				break;
			case EM::multipole:
				Ax = MultipolePrimitive(t,r);
				break;
			case EM::hermite:
				Ax = HermitePrimitive(t,r);
				break;
			case EM::laguerre:
				Ax = LaguerrePrimitive(t,r);
				break;
		}
		// is this thread safe?
		// Ax *= std::exp(ii*randomPhase*deviate->Next());
		switch (modeType)
		{
			case EM::multipole:
				ans = tw::vec3( -std::real(Ax)*r.y/rho , std::real(Ax)*r.x/rho ,0.0 );
				break;
			default:
				ans = tw::vec3( std::real(Ax) , 0.0 ,0.0 );
				break;
		}
		laserFrame.ExpressInStdBasis(&ans);
		return ans;
	}

	void ReadInputFile(std::stringstream& inputString,std::string& command);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
};

struct Pulse:Wave
{
	Pulse(GaussianDeviate *deviate);
	void Initialize();

	// Dispatch function for all mode types
	tw::Complex VectorPotentialEnvelope(tw::Float time,const tw::vec3& pos) const
	{
		tw::Complex Ax;
		tw::vec3 r = pos - focusPosition;
		laserFrame.ExpressInBasis(&r);

		const tw::Float t = time - pulseShape.delay - pulseShape.risetime;
		switch (modeType)
		{
			case EM::plane:
				Ax = PlanePrimitive(t,r);
				break;
			case EM::bessel:
				Ax = BesselPrimitive(t,r);
				break;
			case EM::multipole:
				Ax = MultipolePrimitive(t,r);
				break;
			case EM::hermite:
				Ax = HermitePrimitive(t,r);
				break;
			case EM::laguerre:
				Ax = LaguerrePrimitive(t,r);
				break;
		}
		// is this thread safe?
		// Ax *= std::exp(ii*randomPhase*deviate->Next());
		Ax *= std::exp(-ii*(k0*r.z - w0*t));
		return Ax;
	}
};

struct Conductor
{
	std::vector<Region*>& rgnList;
	Region *theRgn;
	PulseShape pulseShape;
	std::valarray<tw::Float> Px,Py,Pz,potential,angFreq,phase;
	bool affectsPhi,affectsA,electricCurrent,magneticCurrent;
	tw::vec3 gaussianRadius,ks;
	tw::Float f;

	Conductor(std::vector<Region*>& rgnList);
	void Initialize(const MetricSpace& ds);
	tw::Float Voltage(tw::Float t);
	tw::Float VoltageRate(tw::Float t);
	tw::vec3 PolarizationDensity(const tw::vec3& pos,tw::Float t);

	void DepositSources(Field& j4,const MetricSpace& ds,tw::Float t,tw::Float dt);

	void ReadInputFile(std::stringstream& inputString,std::string& command);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
};

struct LindmanBoundary
{
	MetricSpace *ms;
	axisSpec axis;
	sideSpec side;
	Field boundaryMemory;
	tw::Int c0,c1;
	std::vector<Wave*> *wave;

	void Initialize(Task *task,MetricSpace *ms,std::vector<Wave*> *waves,const axisSpec& axis,const sideSpec& side,tw::Int c0,tw::Int c1);
	void UpdateBoundaryMemory(Field& A,tw::Float dt);
	void Set(Field& A,tw::Float t0,tw::Float dt);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
};

struct MABoundary
{
	std::valarray<tw::Complex> g[10];
	tw::Float gamma[10];
	tw::Float sigma[10];
	tw::Float scaleFactor;

	MABoundary(tw::Int numCells);
	void AdvanceLeft(std::valarray<tw::Complex>& amplitude,tw::Float dt);
	void AdvanceRight(std::valarray<tw::Complex>& amplitude,tw::Float dt);
	tw::Complex NormalDerivativeLeft(tw::Int index,tw::Complex amplitude,tw::Float carrierFrequency);
	tw::Complex NormalDerivativeRight(tw::Int index,tw::Complex amplitude,tw::Float carrierFrequency);
};
