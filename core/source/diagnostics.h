// Diagnostics are implemented as ComputeTool objects.
// The flow of control starts from Simulation, which calls Module::Report for every Module.
// The Module::Report function recieves a pointer to a given Diagnostic object (one call for every Diagnostic instance).
// The Module::Report function makes a call to a method of Diagnostic for each data set it would like to report.
// The Diagnostic object is responsible for interpreting the data set, extracting relevant subsets, and writing to disk.

class meta_writer
{
	UnitConverter *units;
public:
	meta_writer(UnitConverter *units);
	std::string s(const std::string& raw);
	void start_entry(const std::string& name,const std::string& diagnostic_name);
	void define_axis(const std::string& name,tw::Int ax,const std::string& label,tw::dimensions units,bool last=false);
	void finish_entry();
};

class npy_writer
{
public:
	std::string form_header(tw::Int shape[4]);
	void write_header(const std::string& name,tw::Int shape[4]);
	void update_shape(const std::string& name,tw::Int shape[4]);
	void add_frame(const std::string& name,const char *gData,tw::Int shape[4]);
};

struct Diagnostic : ComputeTool
{
	std::string filename;
	tw::Int skip[4];
	tw::Float t,tRef,t0,t1,timePeriod,gammaBoost;
	tw::vec3 vGalileo;
	bool headerWritten;

	Diagnostic(const std::string& name,MetricSpace *ms,Task *tsk);
	bool WriteThisStep(tw::Float elapsedTime,tw::Float dt,tw::Int stepNow);
	void StartGridFile(std::ofstream& grid);
	virtual void Start();
	virtual void Finish();
	virtual void Float(const std::string& label,tw::Float val,bool avg);
	virtual void Field(const std::string& fieldName,const struct Field& F,const tw::Int c,const tw::dimensions unit = tw::dimensions::none,const std::string& pretty = "tw::none");
	virtual tw::Float VolumeIntegral(const std::string& fieldName,const struct Field& F,const tw::Int c);
	virtual tw::Float FirstMoment(const std::string& fieldName,const struct Field& F,const tw::Int c,const tw::vec3& r0,const tw::grid::axis axis);
	virtual void Particle(const struct Particle& par,tw::Float m0,tw::Float tp);
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
	virtual void Start();
	virtual void Finish();
	virtual void Float(const std::string& label,tw::Float val,bool avg);
};

struct VolumeDiagnostic : TextTableBase
{
	VolumeDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual tw::Float VolumeIntegral(const std::string& fieldName,const struct Field& F,const tw::Int c);
	virtual tw::Float FirstMoment(const std::string& fieldName,const struct Field& F,const tw::Int c,const tw::vec3& r0,const tw::grid::axis axis);
};

struct PointDiagnostic : TextTableBase
{
	tw::vec3 thePoint;

	PointDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Field(const std::string& fieldName,const struct Field& F,const tw::Int c,const tw::dimensions unit = tw::dimensions::none,const std::string& pretty = "tw::none");
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

struct BoxDiagnostic : Diagnostic
{
	bool average;
	std::vector<std::string> reports;

	BoxDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Finish();
	void GetLocalIndexing(const tw::Int pts[4],const tw::Int glb[6],tw::Int loc[6],const tw::Int coords[4]);
	void GetGlobalIndexing(tw::Int pts[4],tw::Int glb[6]);
	virtual void Field(const std::string& fieldName,const struct Field& F,const tw::Int c,const tw::dimensions unit = tw::dimensions::none,const std::string& pretty = "tw::none");
};

struct PhaseSpaceDiagnostic : Diagnostic
{
	tw::Float bounds[6];
	tw::Int dims[4];
	tw::grid::axis ax[4];
	DiscreteSpace plot_layout;
	ScalarField fxp;

	PhaseSpaceDiagnostic(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Start();
	virtual void Finish();
	virtual void Particle(const struct Particle& par,tw::Float m0,tw::Float tp);
};

struct ParticleOrbits : Diagnostic
{
	tw::Float minGamma;
	std::vector<float> parData;

	ParticleOrbits(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual void Start();
	virtual void Finish();
	virtual void Particle(const struct Particle& par,tw::Float m0,tw::Float tp);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};
