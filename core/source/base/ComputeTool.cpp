#include "meta_base.h"
#include "meta_tools.h"

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

void ComputeTool::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	// Handle any directives that cannot be handled by DirectiveReader
}

void ComputeTool::ReadInputFileBlock(std::stringstream& inputString)
{
	std::string com;
	do
	{
		com = directives.ReadNext(inputString);
		if (com=="tw::EOF")
			throw tw::FatalError("Encountered EOF while processing <"+name+">.");
		ReadInputFileDirective(inputString,com);
	} while (com!="}");
	directives.ThrowErrorIfMissingKeys(name);
}

void ComputeTool::WarningMessage(std::ostream *theStream)
{
	// Save any messages from OpenCL compiler for later output
	if (buildLog.size()>4)
	{
		*theStream << "WARNING : Build log for " << programFilename << " is not empty:" << std::endl;
		*theStream << buildLog << std::endl;
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

tw::tool_type ComputeTool::CreateTypeFromInput(const tw::input::Preamble& preamble)
{
	// Look for a ComputeTool key on a preamble (words between new and opening brace) and return the type of tool.
	const tw::Int max_words = preamble.words.size();
	std::map<std::string,tw::tool_type> tool_map = ComputeTool::Map();
	for (tw::Int i=1;i<=max_words;i++)
	{
		std::string key(tw::input::GetPhrase(preamble.words,i));
		if (tool_map.find(key)!=tool_map.end())
			return tool_map[key];
	}
	return tw::tool_type::none;
}

ComputeTool* ComputeTool::CreateObjectFromType(const std::string& name,tw::tool_type theType,MetricSpace *ms,Task *tsk)
{
	ComputeTool *ans;
	switch (theType)
	{
		case tw::tool_type::none:
			ans = NULL;
			break;
		case tw::tool_type::warp:
			ans = new Warp(name,ms,tsk);
			break;
		case tw::tool_type::conductor:
			ans = new Conductor(name,ms,tsk);
			break;
		case tw::tool_type::planeWave:
			ans = new PlaneWave(name,ms,tsk);
			break;
		case tw::tool_type::hermiteGauss:
			ans = new HermiteGauss(name,ms,tsk);
			break;
		case tw::tool_type::laguerreGauss:
			ans = new LaguerreGauss(name,ms,tsk);
			break;
		case tw::tool_type::besselBeam:
			ans = new BesselBeam(name,ms,tsk);
			break;
		case tw::tool_type::airyDisc:
			ans = new AiryDisc(name,ms,tsk);
			break;
		case tw::tool_type::multipole:
			ans = new Multipole(name,ms,tsk);
			break;
		case tw::tool_type::uniformProfile:
			ans = new UniformProfile(name,ms,tsk);
			break;
		case tw::tool_type::piecewiseProfile:
			ans = new PiecewiseProfile(name,ms,tsk);
			break;
		case tw::tool_type::channelProfile:
			ans = new ChannelProfile(name,ms,tsk);
			break;
		case tw::tool_type::columnProfile:
			ans = new ColumnProfile(name,ms,tsk);
			break;
		case tw::tool_type::gaussianProfile:
			ans = new GaussianProfile(name,ms,tsk);
			break;
		case tw::tool_type::corrugatedProfile:
			ans = new CorrugatedProfile(name,ms,tsk);
			break;
		case tw::tool_type::eigenmodePropagator:
			ans = new EigenmodePropagator(name,ms,tsk);
			break;
		case tw::tool_type::adiPropagator:
			ans = new ADIPropagator(name,ms,tsk);
			break;
		case tw::tool_type::isotropicPropagator:
			ans = new IsotropicPropagator(name,ms,tsk);
			break;
		case tw::tool_type::generalParabolicPropagator:
			ans = new ParabolicSolver(name,ms,tsk);
			break;
		case tw::tool_type::schroedingerPropagator:
			ans = new SchroedingerPropagator(name,ms,tsk);
			break;
		case tw::tool_type::iterativePoissonSolver:
			ans = new IterativePoissonSolver(name,ms,tsk);
			break;
		case tw::tool_type::ellipticSolver1D:
			ans = new EllipticSolver1D(name,ms,tsk);
			break;
		case tw::tool_type::facrPoissonSolver:
			ans = new PoissonSolver(name,ms,tsk);
			break;
		case tw::tool_type::eigenmodePoissonSolver:
			ans = new EigenmodePoissonSolver(name,ms,tsk);
			break;
		case tw::tool_type::yeePropagatorPML:
			ans = new YeePropagatorPML(name,ms,tsk);
			break;
		case tw::tool_type::lorentzPropagator:
			ans = new LorentzPropagator(name,ms,tsk);
			break;
		case tw::tool_type::eosData:
			ans = new EOSComponent(name,ms,tsk);
			break;
		case tw::tool_type::eosIdealGas:
			ans = new EOSIdealGas(name,ms,tsk);
			break;
		case tw::tool_type::eosHotElectrons:
			ans = new EOSHotElectrons(name,ms,tsk);
			break;
		case tw::tool_type::eosMixture:
			ans = new EOSMixture(name,ms,tsk);
			break;
		case tw::tool_type::eosIdealGasMix:
			ans = new EOSIdealGasMix(name,ms,tsk);
			break;
		case tw::tool_type::eosSimpleMieGruneisen:
			ans = new EOSSimpleMieGruneisen(name,ms,tsk);
			break;
		case tw::tool_type::eosLinearMieGruneisen:
			ans = new EOSLinearMieGruneisen(name,ms,tsk);
			break;
		case tw::tool_type::eosTillotson:
			ans = new EOSTillotson(name,ms,tsk);
			break;
		case tw::tool_type::mpi:
			ans = new Multiphoton(name,ms,tsk);
			break;
		case tw::tool_type::adk:
			ans = new ADK(name,ms,tsk);
			break;
		case tw::tool_type::ppt:
			ans = new PPT(name,ms,tsk);
			break;
		case tw::tool_type::ppt_tunneling:
			ans = new PPT_Tunneling(name,ms,tsk);
			break;
		case tw::tool_type::kyh:
			ans = new KYH(name,ms,tsk);
			break;
		case tw::tool_type::pmpb:
			ans = new PMPB(name,ms,tsk);
			break;
		case tw::tool_type::boxDiagnostic:
			ans = new BoxDiagnostic(name,ms,tsk);
			break;
		case tw::tool_type::pointDiagnostic:
			ans = new PointDiagnostic(name,ms,tsk);
			break;
		case tw::tool_type::volumeDiagnostic:
			ans = new VolumeDiagnostic(name,ms,tsk);
			break;
		case tw::tool_type::particleOrbits:
			ans = new ParticleOrbits(name,ms,tsk);
			break;
		case tw::tool_type::phaseSpaceDiagnostic:
			ans = new PhaseSpaceDiagnostic(name,ms,tsk);
			break;
		case tw::tool_type::boundState:
			ans = new BoundState(name,ms,tsk);
			break;
		case tw::tool_type::freeState:
			ans = new FreeState(name,ms,tsk);
			break;
		case tw::tool_type::randomState:
			ans = new RandomState(name,ms,tsk);
			break;
		case tw::tool_type::tabulatedState:
			ans = new TabulatedState(name,ms,tsk);
			break;
		case tw::tool_type::photonGenerator:
			ans = new PhotonGenerator(name,ms,tsk);
			break;
		case tw::tool_type::pairCreator:
			ans = new PairCreator(name,ms,tsk);
			break;
		case tw::tool_type::borisMover:
			ans = new BorisMover(name,ms,tsk);
			break;
		case tw::tool_type::hcMover:
			ans = new HCMover(name,ms,tsk);
			break;
		case tw::tool_type::pgcMover:
			ans = new PGCMover(name,ms,tsk);
			break;
		case tw::tool_type::unitaryMover:
			ans = new UnitaryMover(name,ms,tsk);
			break;
		case tw::tool_type::bohmianMover:
			ans = new BohmianMover(name,ms,tsk);
			break;
		case tw::tool_type::photonMover:
			ans = new PhotonMover(name,ms,tsk);
			break;
	}
	return ans;
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
