module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_test.h"

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
import fields;

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
		borisMover,hcMover,pgcMover,unitaryMover,bohmianMover,photonMover
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
	virtual bool Test(tw::Int& id);

	void InitializeCLProgram(const std::string& filename);

	static std::map<std::string,tw::tool_type> Map();
	static tw::tool_type CreateTypeFromInput(const tw::input::Preamble& preamble);
	static bool SetTestGrid(tw::tool_type theType,tw::Int gridId,MetricSpace *ms,Task *tsk);
};

export struct BoundedTool : ComputeTool
{
	tw::bc::fld x0,x1,y0,y1,z0,z1;
	tw::bc::fld x0s,x1s,y0s,y1s,z0s,z1s; // saved BC's

	BoundedTool(const std::string& name,MetricSpace *ms,Task *tsk);
	void SetBoundaryConditions(tw::bc::fld x0,tw::bc::fld x1,tw::bc::fld y0,tw::bc::fld y1,tw::bc::fld z0,tw::bc::fld z1);
	void SaveBoundaryConditions();
	void RestoreBoundaryConditions();
	void SetFieldsBoundaryConditions(Field& F,const Element& e);
};

//////////////////////////
//                      //
//  COMPUTE TOOL CLASS  //
//                      //
//////////////////////////


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

bool ComputeTool::Test(tw::Int& id)
{
	return false; // did not test, id unchanged means no next test
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
		ts_tree_cursor_goto_first_child(curs);
	}
	do
	{
		if (!directives.ReadNext(curs,src))
			ReadInputFileDirective(curs,src);
	} while (ts_tree_cursor_goto_next_sibling(curs));
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
		{"hermite gauss",tw::tool_type::hermiteGauss},
		{"laguerre gauss",tw::tool_type::laguerreGauss},
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
		{"iterative elliptic",tw::tool_type::iterativePoissonSolver},
		{"1d elliptic",tw::tool_type::ellipticSolver1D},
		{"facr elliptic",tw::tool_type::facrPoissonSolver},
		{"eigenmode elliptic",tw::tool_type::eigenmodePoissonSolver},
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
		{"photon mover",tw::tool_type::photonMover}
	};
}

/// @brief weak matching of the input file key (legacy compatibility)
/// @param preamble data extracted while parsing object
/// @return the corresponding tool_type enumeration
tw::tool_type ComputeTool::CreateTypeFromInput(const tw::input::Preamble& preamble)
{
	// strategy is to lop off trailing words until we get a match
	std::map<std::string,tw::tool_type> tool_map = ComputeTool::Map();
	std::string key_now = preamble.obj_key;
	while (true) {
		if (tool_map.find(key_now)!=tool_map.end()) {
			return tool_map[key_now];
		} else {
			auto p = key_now.rfind(' ');
			if (p != std::string::npos) {
				key_now = key_now.substr(0,p);
			} else {
				break;
			}
		}
	}
	return tw::tool_type::none;
}

bool ComputeTool::SetTestGrid(tw::tool_type theType,tw::Int gridId,MetricSpace *ms,Task *tsk)
{
	switch (theType)
	{
		case tw::tool_type::ellipticSolver1D:
			if (gridId>1)
				return false;
			tsk->Initialize(tw::idx4(1,1,1,2).array,tw::idx4(1,1,1,4).array,tw::idx4(0,1,1,0).array);
			ms->Resize(*tsk,tw::vec4(0,0,0,0),tw::vec4(0.1,0.2,0.2,0.8),2);
			return true;
		case tw::tool_type::pgcMover:
		case tw::tool_type::hcMover:
		case tw::tool_type::borisMover:
			if (gridId==1) {
				tsk->Initialize(tw::idx4(1,1,1,2).array,tw::idx4(1,1,1,4).array,tw::idx4(0,1,1,0).array);
				ms->Resize(*tsk,tw::vec4(0,0,0,0),tw::vec4(0.1,0.2,0.2,0.8),2);
				return true;
			} else if (gridId==2) {
				tsk->Initialize(tw::idx4(1,1,1,2).array,tw::idx4(1,4,1,4).array,tw::idx4(0,1,1,0).array);
				ms->Resize(*tsk,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.2,0.8),2);
				return true;
			} else if (gridId==3) {
				tsk->Initialize(tw::idx4(1,1,1,2).array,tw::idx4(1,4,4,4).array,tw::idx4(0,1,1,0).array);
				ms->Resize(*tsk,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.8,0.8),2);
				return true;
			}
			return false;
		default:
			if (gridId>1)
				return false;
			tsk->Initialize(tw::idx4(1,1,1,2).array,tw::idx4(1,4,4,4).array,tw::idx4(0,1,1,0).array);
			ms->Resize(*tsk,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.8,0.8),2);
			return true;
	}
}

// BoundedTool derived class
// Expect that many tools will want basic access to boundary conditions.
// Therefore provide a derived class at a low level to provide this.

BoundedTool::BoundedTool(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	x0 = x1 = y0 = y1 = z0 = z1 = tw::bc::fld::neumannWall;
	directives.Add("xboundary",new tw::input::Enums<tw::bc::fld>(tw::bc::fld_map(),&x0,&x1),false);
	directives.Add("yboundary",new tw::input::Enums<tw::bc::fld>(tw::bc::fld_map(),&y0,&y1),false);
	directives.Add("zboundary",new tw::input::Enums<tw::bc::fld>(tw::bc::fld_map(),&z0,&z1),false);
}

void BoundedTool::SetBoundaryConditions(tw::bc::fld x0,tw::bc::fld x1,tw::bc::fld y0,tw::bc::fld y1,tw::bc::fld z0,tw::bc::fld z1)
{
	this->x0 = x0;
	this->y0 = y0;
	this->z0 = z0;
	this->x1 = x1;
	this->y1 = y1;
	this->z1 = z1;
}

void BoundedTool::SaveBoundaryConditions()
{
	x0s = x0;
	x1s = x1;
	y0s = y0;
	y1s = y1;
	z0s = z0;
	z1s = z1;
}

void BoundedTool::RestoreBoundaryConditions()
{
	x0 = x0s;
	x1 = x1s;
	y0 = y0s;
	y1 = y1s;
	z0 = z0s;
	z1 = z1s;
}

void BoundedTool::SetFieldsBoundaryConditions(Field& F,const Element& e)
{
	F.SetBoundaryConditions(e,tw::grid::x,x0,x1);
	F.SetBoundaryConditions(e,tw::grid::y,y0,y1);
	F.SetBoundaryConditions(e,tw::grid::z,z0,z1);
}
