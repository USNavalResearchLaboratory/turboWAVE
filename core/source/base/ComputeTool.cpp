module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

/// ComputeTool objects provide low level functionality that is accessible to all Modules.
/// They retain pointers to MetricSpace (grid information) and Task (access to MPI communicators).
/// A facility for interfacing with OpenCL kernels is provided.
/// Each class of tool has a unique identifier of type tw::tool_type.
/// Each instance of a tool has a unique name encoded as a string.
/// To avoid duplicate names, there is an automatic name mangling mechanism.
/// Automatic mangling relies on always creating a tool via Simulation::CreateTool
/// Module and ComputeTool are somewhat similar.
/// ComputeTool is intended to be heavy on computations, light on data ownership.
/// Module is intended to own and manage data, while delegating heavy computations.
export module compute_tool;
import base;
import input;
import region;
import static_space;
import metric_space;

export namespace tw
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
		borisMover,hcMover,pgcMover,unitaryMover,bohmianMover,photonMover,
		// Testers
		iteratorTest,metricSpaceTest
	};
}

export struct ComputeTool : Testable
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
	virtual void WarningMessage();
	virtual void StatusMessage(std::ostream *dest) {;}
	virtual void ReadInputFileBlock(TSTreeCursor *curs,const std::string& src);
	virtual bool ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	void InitializeCLProgram(const std::string& filename);

	static std::map<std::string,tw::tool_type> Map();
	static tw::tool_type CreateTypeFromInput(const tw::input::Preamble& preamble);
	/// setup grid, native units, or other environmental factors affecting a test, returns
	/// whether an environment was defined for the given enviro
	static bool SetTestEnvironment(tw::tool_type theType,tw::Int enviro,MetricSpace *ms,Task *tsk);
};

ComputeTool::ComputeTool(const std::string& name,MetricSpace *ms,Task *tsk)
{
	// Do not create ComputeTools directly (unless you are Simulation).
	// Instances are created through Simulation::CreateTool
	this->name = name;
	space = ms;
	task = tsk;
	refCount = 0;
	region_name = "tw::entire";
	theRgn = NULL;
	programFilename = "";
	buildLog = "";
	native = ms->units;
	natural = tw::UnitConverter(tw::units::natural,native);
	atomic = tw::UnitConverter(tw::units::atomic,native);
	plasma = tw::UnitConverter(tw::units::plasma,native);
	cgs = tw::UnitConverter(tw::units::cgs,native);
	mks = tw::UnitConverter(tw::units::mks,native);
	directives.AttachUnits(native);
	directives.Add("clipping region",new tw::input::String(&region_name),false);
}

ComputeTool::~ComputeTool()
{
	#ifdef USE_OPENCL

	if (programFilename!="")
		clReleaseProgram(program);

	#endif
}

void ComputeTool::Initialize()
{
	// Typically constructor should be complete, nothing to do here, unless MPI needed.
}

void ComputeTool::InitializeCLProgram(const std::string& filename)
{
	#ifdef USE_OPENCL
	task->InitializeCLProgram(program,filename,buildLog);
	programFilename = filename;
	#endif
}

/// @brief called if an assignment was not handled normally
/// @param curs on a directive, probably a custom assignment
/// @param src source document
/// @returns whether any directive was handled
bool ComputeTool::ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src)
{
	return false;
}

/// @brief read all directives in the block
/// @param curs can be on block or on first child of block
/// @param src source document
void ComputeTool::ReadInputFileBlock(TSTreeCursor *curs,const std::string& src)
{
	if (tw::input::node_kind(curs)=="block") {
		if (ts_tree_cursor_goto_first_child(curs)) {
			if (tw::input::next_named_node(curs,true)) {
				do
				{
					if (tw::input::node_kind(curs)=="comment") {
						continue;
					}
					if (!directives.ReadNext(curs,src)) {
						if (!ReadInputFileDirective(curs,src)) {
							tw::input::ThrowParsingError(curs,src,"unknown directive");
						}
					}
				} while (tw::input::next_named_node(curs,false));
			}
		}
	}
	directives.ThrowErrorIfMissingKeys(name);
}

void ComputeTool::WarningMessage()
{
	// Save any messages from OpenCL compiler for later output
	if (buildLog.size()>4)
	{
		std::println(std::cout,"{}: Build log for {} is not empty:",term::warning,programFilename);
		std::println(std::cout,"{}",buildLog);
	}
}

void ComputeTool::ReadCheckpoint(std::ifstream& inFile)
{
}

void ComputeTool::WriteCheckpoint(std::ofstream& outFile)
{
	outFile << name << " ";
}

std::map<std::string,tw::tool_type> ComputeTool::Map()
{
	return
	{
		{"warp",tw::tool_type::warp},
		{"conductor",tw::tool_type::conductor},
		{"plane wave",tw::tool_type::planeWave},
		{"hermite gauss pulse",tw::tool_type::hermiteGauss},
		{"laguerre gauss pulse",tw::tool_type::laguerreGauss},
		{"bessel beam",tw::tool_type::besselBeam},
		{"airy disc",tw::tool_type::airyDisc},
		{"multipole",tw::tool_type::multipole},
		{"uniform",tw::tool_type::uniformProfile},
		{"piecewise",tw::tool_type::piecewiseProfile},
		{"channel",tw::tool_type::channelProfile},
		{"column",tw::tool_type::columnProfile},
		{"gaussian",tw::tool_type::gaussianProfile},
		{"corrugated",tw::tool_type::corrugatedProfile},
		{"eigenmode propagator",tw::tool_type::eigenmodePropagator},
		{"adi propagator",tw::tool_type::adiPropagator},
		{"isotropic propagator",tw::tool_type::isotropicPropagator},
		{"parabolic propagator",tw::tool_type::generalParabolicPropagator},
		{"schroedinger propagator",tw::tool_type::schroedingerPropagator},
		{"iterative elliptic solver",tw::tool_type::iterativePoissonSolver},
		{"1d elliptic solver",tw::tool_type::ellipticSolver1D},
		{"facr elliptic solver",tw::tool_type::facrPoissonSolver},
		{"eigenmode elliptic solver",tw::tool_type::eigenmodePoissonSolver},
		{"yee propagator",tw::tool_type::yeePropagatorPML},
		{"lorentz propagator",tw::tool_type::lorentzPropagator},
		{"eos ideal gas tool",tw::tool_type::eosIdealGas},
		{"eos hot",tw::tool_type::eosHotElectrons},
		{"eos mix",tw::tool_type::eosMixture},
		{"eos ideal gas mix",tw::tool_type::eosIdealGasMix},
		{"eos simple mie gruneisen",tw::tool_type::eosSimpleMieGruneisen},
		{"eos linear mie gruneisen",tw::tool_type::eosLinearMieGruneisen},
		{"eos tillotson",tw::tool_type::eosTillotson},
		{"multiphoton ionization",tw::tool_type::mpi},
		{"adk ionization",tw::tool_type::adk},
		{"ppt ionization",tw::tool_type::ppt},
		{"ppt tunneling",tw::tool_type::ppt_tunneling},
		{"kyh ionization",tw::tool_type::kyh},
		{"pmpb ionization",tw::tool_type::pmpb},
		{"box diagnostic",tw::tool_type::boxDiagnostic},
		{"orbit diagnostic",tw::tool_type::particleOrbits},
		{"phase space diagnostic",tw::tool_type::phaseSpaceDiagnostic},
		{"energy diagnostic",tw::tool_type::volumeDiagnostic},
		{"point diagnostic",tw::tool_type::pointDiagnostic},
		{"qstate free",tw::tool_type::freeState},
		{"qstate bound",tw::tool_type::boundState},
		{"qstate random",tw::tool_type::randomState},
		{"qstate tabulated",tw::tool_type::tabulatedState},
		{"photon generator",tw::tool_type::photonGenerator},
		{"pair creator",tw::tool_type::pairCreator},
		{"boris mover",tw::tool_type::borisMover},
		{"hc mover",tw::tool_type::hcMover},
		{"pgc mover",tw::tool_type::pgcMover},
		{"unitary mover",tw::tool_type::unitaryMover},
		{"bohmian mover",tw::tool_type::bohmianMover},
		{"photon mover",tw::tool_type::photonMover},
		{"iterator test",tw::tool_type::iteratorTest},
		{"metric space test",tw::tool_type::metricSpaceTest}
	};
}

/// @brief weak matching of the input file key (legacy compatibility)
/// @param preamble data extracted while parsing object
/// @return the corresponding tool_type enumeration
tw::tool_type ComputeTool::CreateTypeFromInput(const tw::input::Preamble& preamble)
{
	std::map<std::string,tw::tool_type> tool_map = ComputeTool::Map();
	std::stringstream raw_key(preamble.obj_key);
	std::string normalized,word;
	do {
		raw_key >> word;
		if (!normalized.empty()) {
			normalized += " ";
		}
		normalized += word;
	} while (!raw_key.eof());
	if (tool_map.find(normalized)!=tool_map.end()) {
		return tool_map[normalized];
	} else {
		return tw::tool_type::none;
	}
}

bool ComputeTool::SetTestEnvironment(tw::tool_type theType,tw::Int enviro,MetricSpace *ms,Task *tsk)
{
	const tw::node4 cyclic = {0,1,1,0};
	const tw::node4 domains = {1,1,1,2};
	const tw::node5 gdim1d = {1,1,1,4,1};
	const tw::node5 gdim2d = {1,4,1,4,1};
	const tw::node5 gdim3d = {1,4,4,4,1};
	const tw::node4 layers = {0,2,2,2};
	const tw::vec4 gcorner(0,0,0,0);
	tsk->Initialize(domains,cyclic);
	switch (theType)
	{
		case tw::tool_type::ellipticSolver1D:
			if (enviro==1) {
				ms->Resize(tsk,gdim1d,gcorner,tw::vec4(0.1,0.2,0.2,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
		case tw::tool_type::pgcMover:
		case tw::tool_type::hcMover:
		case tw::tool_type::borisMover:
			if (enviro==1) {
				ms->Resize(tsk,gdim1d,gcorner,tw::vec4(0.1,0.2,0.2,0.8),std_packing,layers);
				return true;
			} else if (enviro==2) {
				ms->Resize(tsk,gdim2d,gcorner,tw::vec4(0.1,0.8,0.2,0.8),std_packing,layers);
				return true;
			} else if (enviro==3) {
				ms->Resize(tsk,gdim3d,gcorner,tw::vec4(0.1,0.8,0.8,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
		default:
			if (enviro==1) {
				ms->Resize(tsk,gdim3d,gcorner,tw::vec4(0.1,0.8,0.8,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
	}
}
