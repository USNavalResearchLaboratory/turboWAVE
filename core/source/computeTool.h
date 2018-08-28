// ComputeTool objects provide low level functionality that is accessible to all Modules.
// They retain pointers to MetricSpace (grid information) and Task (access to MPI communicators).
// A facility for interfacing with OpenCL kernels is provided.

// All ComputeTool types must be named in the enumerated type tw_tool.
// These can be any unique name in the global namespace.

// NEVER create ComputeTool instances directly.
// To create an instance of a specific tool, call Grid::AddSharedTool or Grid::AddPrivateTool
// with a tw_tool type-code as the argument.  The former call creates a tool that is shared
// between modules (subsequent calls return a pointer to the existing tool).  The latter creates
// a new tool each time it is called.

// ComputeTools should typically be instantiated in Module constructors.
// Grid::PrepareSimulation calls ComputeTool::Initialize before Module::Initialize.

// Information from Modules is passed to ComputeTools in arguments of methods of derived classes.
// These methods need not follow any particular protocol.
// A method that advances data in time is conventionally called "Advance."

// ComputeTools should be designed so that data that would be needed for a restart is
// external to the tool (i.e., all state data is directly owned by Modules).

// Adding a new type of tool:
// 1. Give it a name (type-code) and add to the tw_tool enumerated type
// 2. Define the interface in a header file
//   2a. Needs at least constructor and some defining method or methods (e.g., Advance)
//   2b. Overriding Initialize is optional
// 3. Define the implementation in a cpp file that includes the new/modified header
//   3a. Module-owned fields can be passed in as arguments by reference
//   3b. space and task members provide grid and domain decomposition information
// 4. Add a case to static member CreateObjectFromType

enum tw_tool {		nullTool,
					eigenmodePropagator, adiPropagator, isotropicPropagator,
					iterativePoissonSolver, facrPoissonSolver, eigenmodePoissonSolver, ellipticSolver1D,
					generalParabolicPropagator, schroedingerPropagator,
					yeePropagatorPML, lorentzPropagator, eosDataTool, eosIdealGas, eosIdealGasElectrons, // ASHER_MOD
					eosMixture, eosIdealGasMix, eosMieGruneisen, eosMieGruneisen2  };

struct ComputeTool
{
	MetricSpace *space;
	Task *task;
	tw_tool typeCode;
	std::string name;
	bool sharedTool; // if true only one instance per MPI task
	int refCount; // used to manage releasing shared tools from memory

	// OpenCL Support
	private:
	std::string programFilename;
	std::string buildLog;
	public:
	#ifdef USE_OPENCL
	cl_program program;
	#endif

	ComputeTool(MetricSpace *ms,Task *tsk,bool shared);
	virtual ~ComputeTool();
	virtual void Initialize();
	void InitializeCLProgram(const std::string& filename);
	static ComputeTool* CreateObjectFromType(tw_tool theType,MetricSpace *ms,Task *tsk,bool shared);
	virtual void WarningMessage(std::ostream *theStream);
	virtual void StatusMessage(std::ostream *theStream) {;}
};
