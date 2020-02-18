// ComputeTool objects provide low level functionality that is accessible to all Modules.
// They retain pointers to MetricSpace (grid information) and Task (access to MPI communicators).
// A facility for interfacing with OpenCL kernels is provided.

// Each class of tool has a unique identifier of type tw::tool_type.
// Each instance of a tool has a unique name encoded as a string.
// To avoid duplicate names, there is an automatic name mangling mechanism.
// Automatic mangling relies on always creating a tool via Simulation::CreateTool

// Module and ComputeTool are somewhat similar.
// ComputeTool is intended to be heavy on computations, light on data ownership.
// Module is intended to own and manage data, while delegating heavy computations.

namespace tw
{
	enum class tool_type {		nullTool,
		uniformProfile,channelProfile,gaussianProfile,columnProfile,piecewiseProfile,corrugatedProfile,
		eigenmodePropagator, adiPropagator, isotropicPropagator,
		iterativePoissonSolver, facrPoissonSolver, eigenmodePoissonSolver, ellipticSolver1D,
		generalParabolicPropagator, schroedingerPropagator,
		yeePropagatorPML, lorentzPropagator, eosData, eosIdealGas, eosHotElectrons, // ASHER_MOD
		eosMixture, eosIdealGasMix, eosSimpleMieGruneisen, eosLinearMieGruneisen  };
}

struct ComputeTool
{
	MetricSpace *space;
	Task *task;
	tw::tool_type typeCode;
	std::string name;
	int refCount; // how many modules currently using
	tw::input::DirectiveReader directives;

	// Allow for a region
	std::string region_name;
	Region *theRgn;

	// OpenCL Support
	private:
	std::string programFilename;
	std::string buildLog;
	public:
	#ifdef USE_OPENCL
	cl_program program;
	#endif

	ComputeTool(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual ~ComputeTool();
	virtual void Initialize();
	void InitializeCLProgram(const std::string& filename);
	virtual void WarningMessage(std::ostream *theStream);
	virtual void StatusMessage(std::ostream *theStream) {;}
	virtual void ReadInputFileBlock(std::stringstream& inputString);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
	virtual void SaveToolReference(std::ofstream& outFile);

	static std::map<std::string,tw::tool_type> Map();
	static tw::tool_type CreateTypeFromInput(const std::vector<std::string>& preamble);
	static ComputeTool* CreateObjectFromType(const std::string& name,tw::tool_type theType,MetricSpace *ms,Task *tsk);
	static ComputeTool* CreateObjectFromFile(std::ifstream& inFile,MetricSpace *ms,Task *tsk);
};
