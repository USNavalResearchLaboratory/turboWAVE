namespace tw
{
	namespace profile
	{
		enum class quantity { density,energy,px,py,pz };
		enum class timing { triggered,maintained };
		enum class loading { statistical,deterministic };
		enum class shape { triangle,sin2,quartic,quintic,sech };
	}
}

namespace EM
{
	enum mode { plane , hermite , laguerre , airy_disc , bessel , multipole };
	struct mode_data
	{
		tw::Int order,exponent;
		tw::Float scale;
	};
}

struct Profile : ComputeTool
{
	tw::profile::shape segmentShape;
	tw::dom::geometry symmetry;
	tw::vec3 centerPt;
	tw::vec3 modeNumber;
	tw::Float modeAmplitude;
	tw::basis orientation;
	tw::Float gamma_boost;

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

	Profile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	tw::vec3 DriftMomentum(const tw::Float& mass);
	tw::vec3 Boost(const tw::vec3& pos);
	tw::vec3 Translate_Rotate(const tw::vec3& pos);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct UniformProfile:Profile
{
	tw::Float density;

	UniformProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct GaussianProfile:Profile
{
	tw::Float density;
	tw::vec3 beamSize;

	GaussianProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct ChannelProfile:Profile
{
	std::valarray<tw::Float> z,fz;
	tw::Float coeff[4];

	ChannelProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct ColumnProfile:Profile
{
	std::valarray<tw::Float> z,fz;
	tw::vec3 beamSize;

	ColumnProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct PiecewiseProfile:Profile
{
	std::valarray<tw::Float> x,fx,y,fy,z,fz;

	PiecewiseProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct CorrugatedProfile:Profile
{
	tw::Float a0,gamma0,w0,wp0,km,delta,rchannel;

	CorrugatedProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct PulseShape
{
	tw::Float delay,risetime,holdtime,falltime;
	tw::Float t1,t2,t3,t4;
	tw::profile::shape whichProfile;

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
	tw::Float a0,w,nrefr,phase,vg,chirp,randomPhase,gamma_boost;
	PulseShape pulseShape;
	EM::mode modeType;
	EM::mode_data modeData[2];
	tw::basis laserFrame;

	GaussianDeviate *deviate;

	Wave(GaussianDeviate *deviate);
	void Initialize();

	tw::Complex PlanePrimitive(tw::Float t,const tw::vec3& pos) const;
	tw::Complex BesselPrimitive(tw::Float t,const tw::vec3& pos) const;
	tw::Complex AiryDiscPrimitive(tw::Float t,const tw::vec3& pos) const;
	tw::Complex MultipolePrimitive(tw::Float t,const tw::vec3& pos) const;
	tw::Complex HermitePrimitive(tw::Float t,const tw::vec3& pos) const;
	tw::Complex LaguerrePrimitive(tw::Float t,const tw::vec3& pos) const;

	// Dispatch function for all mode types
	tw::vec3 VectorPotential(tw::Float time,const tw::vec3& pos) const
	{
		tw::vec3 ans;
		tw::Complex Ax;
		tw::vec4 x4(time,pos);
		// N.b. at present boost only works if polarization is orthogonal to z.
		// Otherwise we would have to bring in the scalar potential.
		// The function's caller is giving us boosted frame coordinates.
		// The user is giving us lab frame coordinates.
		// Therefore first transform arguments to lab frame, then proceed as usual.
		x4.zBoost(gamma_boost,1.0);
		tw::vec3 r = x4.spatial() - focusPosition;
		laserFrame.ExpressInBasis(&r);

		const tw::Float rho = tw::small_pos + sqrt(r.x*r.x+r.y*r.y);
		const tw::Float t = x4[0] - pulseShape.delay - pulseShape.risetime;
		switch (modeType)
		{
			case EM::plane:
				Ax = PlanePrimitive(t,r);
				break;
			case EM::bessel:
				Ax = BesselPrimitive(t,r);
				break;
			case EM::airy_disc:
				Ax = AiryDiscPrimitive(t,r);
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
	Pulse(GaussianDeviate *deviate) : Wave(deviate)
	{
	}
	// Dispatch function for all mode types
	tw::Complex VectorPotentialEnvelope(tw::Float time,const tw::vec3& pos,tw::Float w0) const
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
			case EM::airy_disc:
				Ax = AiryDiscPrimitive(t,r);
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
		Ax *= std::exp(-ii*(w0*nrefr*r.z - w0*t));
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
	tw::dom::axis axis;
	tw::dom::side side;
	Field boundaryMemory;
	tw::Int c0,c1;
	std::vector<Wave*> *wave;

	void Initialize(Task *task,MetricSpace *ms,std::vector<Wave*> *waves,const tw::dom::axis& axis,const tw::dom::side& side,tw::Int c0,tw::Int c1);
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
