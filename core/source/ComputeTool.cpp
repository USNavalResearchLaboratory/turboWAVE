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
	typeCode = tw::tool_type::none;
	refCount = 0;
	region_name = "tw::entire";
	theRgn = NULL;
	programFilename = "";
	buildLog = "";
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
	// Get the region name.  Finding the Region object is deferred to owning Module.
	inFile >> region_name;
	inFile.ignore();
}

void ComputeTool::WriteCheckpoint(std::ofstream& outFile)
{
	outFile.write((char *)&typeCode,sizeof(tw::tool_type));
	outFile << name << " ";
	// Save the region name.  The Region itself is managed (and saved) by Simulation.
	outFile << region_name << " ";
}

void ComputeTool::SaveToolReference(std::ofstream& outFile)
{
	// Parent object should call this to save reference to tool in restart file
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
		{"mpi ionization",tw::tool_type::mpi},
		{"adk ionization",tw::tool_type::adk},
		{"ppt ionization",tw::tool_type::ppt},
		{"box diagnostic",tw::tool_type::boxDiagnostic},
		{"orbit diagnostic",tw::tool_type::particleOrbits},
		{"phase space diagnostic",tw::tool_type::phaseSpaceDiagnostic},
		{"energy diagnostic",tw::tool_type::volumeDiagnostic},
		{"point diagnostic",tw::tool_type::pointDiagnostic},
		{"qstate free",tw::tool_type::freeState},
		{"qstate bound",tw::tool_type::boundState},
		{"qstate random",tw::tool_type::randomState},
		{"qstate tabulated",tw::tool_type::tabulatedState}
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
		case tw::tool_type::mpi:
			ans = new MPI(name,ms,tsk);
			break;
		case tw::tool_type::adk:
			ans = new ADK(name,ms,tsk);
			break;
		case tw::tool_type::ppt:
			ans = new PPT(name,ms,tsk);
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
	}
	return ans;
}

ComputeTool* ComputeTool::CreateObjectFromFile(std::ifstream& inFile,MetricSpace *ms,Task *tsk)
{
	tw::tool_type theType;
	std::string name;
	ComputeTool *ans;
	inFile.read((char*)&theType,sizeof(tw::tool_type));
	inFile >> name;
	inFile.ignore();
	ans = CreateObjectFromType(name,theType,ms,tsk);
	ans->ReadCheckpoint(inFile);
	return ans;
}
