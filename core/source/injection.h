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
	enum class current { none,electric,magnetic };
	struct mode_data
	{
		tw::Int order[2],exponent[2];
		tw::Float scale[2];
	};
}

struct Profile : ComputeTool
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

	Profile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	tw::vec3 DriftMomentum(const tw::Float& mass);
	tw::vec3 Boost(const tw::vec3& pos);
	tw::vec3 Translate_Rotate(const tw::vec3& pos);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct UniformProfile:Profile
{
	tw::Float density;

	UniformProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct GaussianProfile:Profile
{
	tw::Float density;
	tw::vec3 beamSize;

	GaussianProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct ChannelProfile:Profile
{
	std::valarray<tw::Float> z,fz;
	tw::Float coeff[4];

	ChannelProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct ColumnProfile:Profile
{
	std::valarray<tw::Float> z,fz;
	tw::vec3 beamSize;

	ColumnProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct PiecewiseProfile:Profile
{
	std::valarray<tw::Float> x,fx,y,fy,z,fz;

	PiecewiseProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct CorrugatedProfile:Profile
{
	tw::Float a0,gamma0,w0,wp0,km,delta,rchannel;

	CorrugatedProfile(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Float GetValue(const tw::vec3& pos,const MetricSpace& ds);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct PulseShape
{
	tw::Float delay,risetime,holdtime,falltime;
	tw::Float t1,t2,t3,t4;
	tw::profile::shape whichProfile;

	PulseShape();
	void Initialize(const tw::Float time_origin);
	tw::Float PulseShapeFactor(const tw::Float t) const;
	tw::Float D1Intensity(const tw::Float t) const;
	tw::Float D1Amplitude(const tw::Float t) const;
};

struct Wave : ComputeTool
{
	tw::vec3 direction;
	tw::vec3 focusPosition;
	tw::vec3 a;
	tw::Float a0,w,nrefr,phase,vg,chirp,randomPhase,gammaBoost;
	PulseShape pulseShape;
	EM::mode_data modeData;
	tw::basis laserFrame;

	GaussianDeviate *deviate;

	Wave(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Initialize();
	virtual tw::Complex PrimitivePhasor(const tw::vec4& x4) const = 0;
	virtual tw::vec3 PrimitiveVector(const tw::vec4& x4) const
	{
		// Estimate from kz*Az = i*dAx/dx
		// This is refined later by field solvers
		tw::Float kz = w*nrefr;
		tw::Float dx = 0.01*modeData.scale[0];
		tw::vec4 dxh(0.0,0.5*dx,0.0,0.0);
		tw::Complex Ax1 = PrimitivePhasor(x4-dxh);
		tw::Complex Ax2 = PrimitivePhasor(x4+dxh);
		tw::Complex Az = ii*(Ax2-Ax1)/(dx*kz);
		return tw::vec3(std::real(0.5*(Ax1+Ax2)),0.0,std::real(Az));
	}
	void ToLaserFrame(tw::vec4 *x4) const
	{
		// The function's caller is giving us boosted frame coordinates.
		// The user is describing the laser in the lab frame.
		// To get to laser's frame: boost, then translate, then rotate.
		x4->zBoost(gammaBoost,1.0);
		tw::vec4 displ(pulseShape.delay+pulseShape.risetime,focusPosition);
		*x4 -= displ;
		laserFrame.ExpressInBasis(x4);
	}
	void ToBoostedFrame(tw::vec4 *A4) const
	{
		laserFrame.ExpressInStdBasis(A4);
		A4->zBoost(gammaBoost,-1.0);
	}
	tw::Complex VectorPotentialEnvelope(tw::Float time,const tw::vec3& pos,tw::Float w0) const
	{
		tw::vec4 x4(time,pos);
		ToLaserFrame(&x4);
		return PrimitivePhasor(x4)*std::exp(-ii*(w0*nrefr*x4[3] - w0*x4[0]));
	}
	tw::vec3 VectorPotential(tw::Float time,const tw::vec3& pos) const
	{
		tw::vec4 x4(time,pos);
		ToLaserFrame(&x4);
		tw::vec4 A4(0.0,PrimitiveVector(x4));
		ToBoostedFrame(&A4);
		return A4.spatial(); // For certain boost geometries we will need to keep scalar potential
	}

	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct PlaneWave : Wave
{
	PlaneWave(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Complex PrimitivePhasor(const tw::vec4& x4) const;
	virtual tw::vec3 PrimitiveVector(const tw::vec4& x4) const;
};

struct BesselBeam : Wave
{
	BesselBeam(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Complex PrimitivePhasor(const tw::vec4& x4) const;
};

struct AiryDisc : Wave
{
	AiryDisc(const std::string& name,MetricSpace *m,Task *tsk);
	virtual tw::Complex PrimitivePhasor(const tw::vec4& x4) const;
};

struct Multipole : Wave
{
	Multipole(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Complex PrimitivePhasor(const tw::vec4& x4) const;
	virtual tw::vec3 PrimitiveVector(const tw::vec4& x4) const;
};

struct HermiteGauss : Wave
{
	HermiteGauss(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Complex PrimitivePhasor(const tw::vec4& x4) const;
};

struct LaguerreGauss : Wave
{
	LaguerreGauss(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Complex PrimitivePhasor(const tw::vec4& x4) const;
};

struct Conductor : ComputeTool
{
	PulseShape pulseShape;
	std::valarray<tw::Float> Px,Py,Pz,potential,angFreq,phase;
	bool affectsPhi,affectsA;
	EM::current currentType;
	tw::vec3 gaussianRadius,ks;
	tw::Float f;

	Conductor(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	tw::Float Voltage(tw::Float t);
	tw::Float VoltageRate(tw::Float t);
	tw::vec3 PolarizationDensity(const tw::vec3& pos,tw::Float t);
	void DepositSources(Field& j4,tw::Float t,tw::Float dt);

	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct LindmanBoundary
{
	MetricSpace *ms;
	tw::grid::axis axis;
	tw::grid::side side;
	Field boundaryMemory;
	tw::Int c0,c1;
	std::vector<Wave*> *wave;

	void Initialize(Task *task,MetricSpace *ms,std::vector<Wave*> *waves,const tw::grid::axis& axis,const tw::grid::side& side,tw::Int c0,tw::Int c1);
	void UpdateBoundaryMemory(Field& A,tw::Float dt);
	void Set(Field& A,tw::Float t0,tw::Float dt);
	void ReadCheckpoint(std::ifstream& inFile);
	void WriteCheckpoint(std::ofstream& outFile);
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

struct Warp : ComputeTool
{
	tw::grid::axis ax;
	bool increasing;
	tw::Int rng[2];
	tw::Float L,gridSum;

	Warp(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
	tw::Float AddedCellWidth(tw::Int globalCell);
};
