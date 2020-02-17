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
	typeCode = tw::tool_type::nullTool;
	refCount = 0;
	programFilename = "";
	buildLog = "";
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
		ReadInputFileDirective(inputString,com);
	} while (com!="}");
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

void ComputeTool::ReadData(std::ifstream& inFile)
{
}

void ComputeTool::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&typeCode,sizeof(tw::tool_type));
	outFile << name << " ";
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
		{"eos linear mie gruneisen",tw::tool_type::eosLinearMieGruneisen}
	};
}

tw::tool_type ComputeTool::CreateTypeFromInput(const std::vector<std::string>& preamble)
{
	// Look for a ComputeTool key on a preamble (words between new and opening brace) and return the type of tool.
	std::map<std::string,tw::tool_type> tool_map = ComputeTool::Map();
	std::string key2(tw::input::GetPhrase(preamble,2));
	std::string key4(tw::input::GetPhrase(preamble,4));
	if (tool_map.find(key2)!=tool_map.end())
		return tool_map[key2];
	if (tool_map.find(key4)!=tool_map.end())
		return tool_map[key4];
	return tw::tool_type::nullTool;
}

ComputeTool* ComputeTool::CreateObjectFromType(const std::string& name,tw::tool_type theType,MetricSpace *ms,Task *tsk)
{
	ComputeTool *ans;
	switch (theType)
	{
		case tw::tool_type::nullTool:
			ans = NULL;
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
	ans->ReadData(inFile);
	return ans;
}
