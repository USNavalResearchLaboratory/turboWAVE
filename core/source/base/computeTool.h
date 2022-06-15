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
	enum class tool_type
	{
		none,warp,
		// Profiles
		uniformProfile,channelProfile,gaussianProfile,columnProfile,piecewiseProfile,corrugatedProfile,
		// Wave launchers
		conductor, planeWave, besselBeam, airyDisc, hermiteGauss, laguerreGauss, multipole,
		// Laser propagators
		eigenmodePropagator, adiPropagator, isotropicPropagator, schroedingerPropagator,
		// Elliptic solvers
		iterativePoissonSolver, facrPoissonSolver, eigenmodePoissonSolver, ellipticSolver1D,
		// Diffusion
		generalParabolicPropagator,
		// Hyperbolic solvers
		yeePropagatorPML, lorentzPropagator,
		// Equation of state
		eosData, eosIdealGas, eosHotElectrons, eosMixture, eosIdealGasMix, eosSimpleMieGruneisen, eosLinearMieGruneisen,
		eosTillotson,
		// Ionization
		mpi, adk, kyh, ppt_tunneling, ppt, pmpb,
		// Diagnostics
		boxDiagnostic,particleOrbits,phaseSpaceDiagnostic,volumeDiagnostic,pointDiagnostic,
		// Quantum
		randomState,freeState,boundState,tabulatedState,
		// Quantum Electrodynamics
		photonGenerator,pairCreator,
		// Movers
		borisMover,pgcMover,unitaryMover,bohmianMover,photonMover
	};
}

struct ComputeTool
{
	MetricSpace *space;
	Task *task;
	std::string name;
	int refCount; // how many modules currently using
	tw::input::DirectiveReader directives;
	tw::UnitConverter native,natural,atomic,plasma,cgs,mks;

	// Allow for a region
	std::string region_name;
	Region *theRgn;

	// Test info
	std::string testName;

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
	virtual void WarningMessage(std::ostream *theStream);
	virtual void StatusMessage(std::ostream *theStream) {;}
	virtual void ReadInputFileBlock(std::stringstream& inputString);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
	virtual bool Test(tw::Int& id);

	void InitializeCLProgram(const std::string& filename);

	static std::map<std::string,tw::tool_type> Map();
	static tw::tool_type CreateTypeFromInput(const tw::input::Preamble& preamble);
	static ComputeTool* CreateObjectFromType(const std::string& name,tw::tool_type theType,MetricSpace *ms,Task *tsk);
	static bool SetTestGrid(tw::tool_type theType,tw::Int gridId,MetricSpace *ms,Task *tsk);
};

struct BoundedTool : ComputeTool
{
	tw::bc::fld x0,x1,y0,y1,z0,z1;
	tw::bc::fld x0s,x1s,y0s,y1s,z0s,z1s; // saved BC's

	BoundedTool(const std::string& name,MetricSpace *ms,Task *tsk);
	void SetBoundaryConditions(tw::bc::fld x0,tw::bc::fld x1,tw::bc::fld y0,tw::bc::fld y1,tw::bc::fld z0,tw::bc::fld z1);
	void SaveBoundaryConditions();
	void RestoreBoundaryConditions();
	void SetFieldsBoundaryConditions(Field& F,const Element& e);
};
