module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module driver:tool;

import base;
import input;
import static_space;
import metric_space;

export namespace tw
{
	enum class tool_type
	{
		none,warp,
		// Primitive regions
		entireRegion , rectRegion , prismRegion , circRegion , cylinderRegion , cylindricalShellRegion ,
		roundedCylinderRegion , ellipsoidRegion , boxArrayRegion , torusRegion , coneRegion , tangentOgiveRegion,
		// Regions
		unionRegion, intersectionRegion, differenceRegion,
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
		iteratorTest,metricSpaceTest,fftTest,regionsTest,
		// Drivers
		electrostatic,
		coulombSolver,directSolver,curvilinearDirectSolver,farFieldDiagnostic,
		qsLaser,pgcLaser,
		kinetics,species,fluidFields,equilibriumGroup,chemical,sparcHydroManager,
		boundElectrons,schroedinger,pauli,kleinGordon,dirac,populationDiagnostic
	};
	bool IsDriver(tool_type typ) {
		return typ==tool_type::electrostatic ||
			typ==tool_type::coulombSolver ||
			typ==tool_type::directSolver ||
			typ==tool_type::curvilinearDirectSolver ||
			typ==tool_type::farFieldDiagnostic ||
			typ==tool_type::qsLaser ||
			typ==tool_type::pgcLaser ||
			typ==tool_type::kinetics ||
			typ==tool_type::species ||
			typ==tool_type::fluidFields ||
			typ==tool_type::equilibriumGroup ||
			typ==tool_type::chemical ||
			typ==tool_type::sparcHydroManager ||
			typ==tool_type::boundElectrons ||
			typ==tool_type::schroedinger ||
			typ==tool_type::pauli ||
			typ==tool_type::kleinGordon ||
			typ==tool_type::dirac ||
			typ==tool_type::populationDiagnostic;
	}
}

/// ComputeTool objects are the root of the simulation management hierarchy.
/// They provide basic functions like testing, input file parsing, checkpointing, and
/// offloading calculations to a device.
export struct ComputeTool : Testable
{
	MetricSpace *space;
	Task *task;
	std::string name;
	tw::input::DirectiveReader directives;
	tw::UnitConverter native,natural,atomic,plasma,cgs,mks;

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
	virtual void StatusMessage(std::ostream *dest) {;}
	virtual void ReadInputFileBlock(TSTreeCursor *curs,const std::string& src);
	virtual bool ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	void InitializeCLProgram(const std::string& filename);

	static std::map<std::string,tw::tool_type> Map();
	static std::map<tw::tool_type,std::string> InvMap();
	static tw::tool_type CreateTypeFromInput(const tw::input::Preamble& preamble);
	/// setup grid, native units, or other environmental factors affecting a test, returns
	/// whether an environment was defined for the given enviro
	static bool SetTestEnvironment(tw::tool_type theType,tw::Int enviro,MetricSpace *ms,Task *tsk,tw::Float unitDensityCGS);
};

export using SharedTool = std::shared_ptr<ComputeTool>;

ComputeTool::ComputeTool(const std::string& name,MetricSpace *ms,Task *tsk)
{
	this->name = name;
	space = ms;
	task = tsk;
	programFilename = "";
	buildLog = "";
	native = ms->unitConverter;
	natural = tw::UnitConverter(tw::units::natural,native.UnitDensityCGS());
	atomic = tw::UnitConverter(tw::units::atomic,native.UnitDensityCGS());
	plasma = tw::UnitConverter(tw::units::plasma,native.UnitDensityCGS());
	cgs = tw::UnitConverter(tw::units::cgs,native.UnitDensityCGS());
	mks = tw::UnitConverter(tw::units::mks,native.UnitDensityCGS());
	directives.AttachUnits(native.unit_system,native.UnitDensityCGS());
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

		{"region box_array",tw::tool_type::boxArrayRegion},
		{"region circ",tw::tool_type::circRegion},
		{"region cone",tw::tool_type::coneRegion},
		{"region cylinder",tw::tool_type::cylinderRegion},
		{"region cylindrical_shell",tw::tool_type::cylindricalShellRegion},
		{"region ellipsoid",tw::tool_type::ellipsoidRegion},
		{"region prism",tw::tool_type::prismRegion},
		{"region rect",tw::tool_type::rectRegion},
		{"region rounded_cylinder",tw::tool_type::roundedCylinderRegion},
		{"region tangent_ogive",tw::tool_type::tangentOgiveRegion},
		{"region torus",tw::tool_type::torusRegion},

		{"region union",tw::tool_type::unionRegion},
		{"region intersection",tw::tool_type::intersectionRegion},
		{"region difference",tw::tool_type::differenceRegion},

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
		{"metric space test",tw::tool_type::metricSpaceTest},
		{"fft test",tw::tool_type::fftTest},
		{"regions test",tw::tool_type::regionsTest},

		{"maxwell solver",tw::tool_type::directSolver},
		{"curvilinear maxwell solver",tw::tool_type::curvilinearDirectSolver},
		{"coulomb gauge solver",tw::tool_type::coulombSolver},
		{"far field diagnostic",tw::tool_type::farFieldDiagnostic},
		{"electrostatic solver",tw::tool_type::electrostatic},
		{"quasistatic solver",tw::tool_type::qsLaser},
		{"pgc laser solver",tw::tool_type::pgcLaser},
		{"bound solver",tw::tool_type::boundElectrons},
		{"schroedinger solver",tw::tool_type::schroedinger},
		{"pauli solver",tw::tool_type::pauli},
		{"klein gordon solver",tw::tool_type::kleinGordon},
		{"dirac solver",tw::tool_type::dirac},
		{"fluid",tw::tool_type::fluidFields},
		{"hydrodynamics",tw::tool_type::sparcHydroManager},
		{"group",tw::tool_type::equilibriumGroup},
		{"chemical",tw::tool_type::chemical},
		{"species",tw::tool_type::species},
		{"kinetics",tw::tool_type::kinetics},
		{"population diagnostic",tw::tool_type::populationDiagnostic}
	};
}

std::map<tw::tool_type,std::string> ComputeTool::InvMap() {
	auto ans = std::map<tw::tool_type,std::string>();
	for (auto pair : Map()) {
		ans[pair.second] = pair.first;
	}
	return ans;
}

/// @brief match input file key after normalizing white space
/// @param preamble data extracted while parsing object
/// @return the corresponding module_type enumeration, may return tw::module_type::none
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

bool ComputeTool::SetTestEnvironment(tw::tool_type theType,tw::Int enviro,MetricSpace *ms,Task *tsk,tw::Float unitDensityCGS)
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
		case tw::tool_type::species:
			if (enviro==1) {
				ms->AttachUnits(tw::units::plasma,unitDensityCGS);
				ms->Resize(tsk,gdim2d,gcorner,tw::vec4(0.1,0.8,0.2,0.8),std_packing,layers);
				return true;
			} else if (enviro==2) {
				ms->AttachUnits(tw::units::plasma,unitDensityCGS);
				ms->Resize(tsk,gdim3d,gcorner,tw::vec4(0.1,0.8,0.8,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
		case tw::tool_type::pauli:
			return false;
		case tw::tool_type::schroedinger:
			if (enviro==1) {
				ms->AttachUnits(tw::units::atomic,unitDensityCGS);
				ms->Resize(tsk,gdim3d,gcorner,tw::vec4(0.1,0.8,0.8,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
		case tw::tool_type::dirac:
		case tw::tool_type::kleinGordon:
			if (enviro==1) {
				ms->AttachUnits(tw::units::natural,unitDensityCGS);
				ms->Resize(tsk,gdim3d,gcorner,tw::vec4(0.1,0.8,0.8,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
		default:
			if (enviro==1) {
				ms->AttachUnits(tw::units::plasma,unitDensityCGS);
				ms->Resize(tsk,gdim3d,gcorner,tw::vec4(0.1,0.8,0.8,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
	}
}
