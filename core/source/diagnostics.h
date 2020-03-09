// Diagnostics are implemented as ComputeTool objects.
// The flow of control starts from Simulation, which calls Module::Report for every Module.
// The Module::Report function recieves a pointer to a given Diagnostic object (one call for every Diagnostic instance).
// The Module::Report function makes a call to a method of Diagnostic for each data set it would like to report.
// The Diagnostic object is responsible for interpreting the data set, extracting relevant subsets, and writing to disk.

void WriteDVHeader(std::ofstream& outFile,tw::Int version,tw::Int xDim,tw::Int yDim,tw::Int zDim,float x0,float x1,float y0,float y1,float z0,float z1);

struct Diagnostic : ComputeTool
{
	std::string filename;
	tw::Int skip[4];
	tw::Float tRef,t0,t1,timePeriod;
	bool moveWithWindow,writeHeader;

	DiagnosticDescriptor(const std::string& name,MetricSpace *ms,Task *tsk);
	bool WriteThisStep(tw::Float elapsedTime,tw::Float dt,tw::Int stepNow);
	virtual void Start(bool writeHeader);
	virtual void Finish();
	virtual void Float(const std::string& label,tw::Float val);
	virtual void Field(const std::string& fieldName,const Field& F,const tw::Int c);
	virtual tw::Float VolumeIntegral(const std::string& fieldName,const Field& F,const tw::Int c);
	virtual void Particle(const Particle& par,tw::Float m0,tw::Float t);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct TextTableBase : Diagnostic
{
	tw::Int numSigFigs;
	std::vector<std::string> labels;
	std::vector<tw::Float> values;
	std::vector<bool> avg;

	TextTableBase(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Start(bool writeHeader);
	virtual void Finish();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct VolumeDiagnostic : TextTableBase
{
	VolumeDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Float(const std::string& label,tw::Float val);
	virtual tw::Float VolumeIntegral(const std::string& fieldName,const Field& F,const tw::Int c);
};

struct PointDiagnostic : TextTableBase
{
	tw::vec3 thePoint;

	PointDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Field(const std::string& fieldName,const Field& F,const tw::Int c);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct BoxDiagnostic : Diagnostic
{
	bool average;

	BoxDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk);
	void GetLocalIndexing(const tw::Int pts[4],const tw::Int glb[6],tw::Int loc[6],const tw::Int coords[4]);
	void GetGlobalIndexing(tw::Int pts[4],tw::Int glb[6]);
	virtual void Field(const std::string& fieldName,const Field& F,const tw::Int c);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct PhaseSpaceDiagnostic : Diagnostic
{
	tw::Float bounds[6];
	tw::Int dims[4];
	tw::grid::axis ax[4];
	ScalarField fxp;

	PhaseSpaceDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Start(bool writeHeader);
	virtual void Finish();
	virtual void Particle(const Particle& par,tw::Float m0,tw::Float t);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct ParticleOrbits : Diagnostic
{
	tw::Float minGamma;
	std::vector<float> parData;

	ParticleOrbits(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Start(bool writeHeader);
	virtual void Finish();
	virtual void Particle(const Particle& par,tw::Float m0,tw::Float t);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

// struct FarFieldDetector : Diagnostic
// {
// 	tw::Float radius,theta0,theta1,phi0,phi1;
// 	tw::Int thetaPts,phiPts,timePts;
// 	std::ofstream AthetaFile,AphiFile;
// 	Vec3Field A; // index as (t,theta,phi)
//
// 	FarFieldDetector(const std::string& name,MetricSpace *ms,Task *tsk);
// 	virtual void ReadCheckpoint(std::ifstream& inFile);
// 	virtual void WriteCheckpoint(std::ofstream& outFile);
// 	void AccumulateField(const tw::Float& elapsedTime,Field& J4);
// };
