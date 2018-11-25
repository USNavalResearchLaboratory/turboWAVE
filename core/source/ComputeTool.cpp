#include "definitions.h"
#include "tasks.h"
#include "ctools.h"
#include "3dmath.h"
#include "metricSpace.h"
#include "3dfields.h"
#include "region.h"
#include "numerics.h"
#include "computeTool.h"
#include "parabolic.h"
#include "elliptic.h"
#include "hyperbolic.h"

#include "physics.h" // ASHER_MOD

//////////////////////////
//                      //
//  COMPUTE TOOL CLASS  //
//                      //
//////////////////////////


ComputeTool::ComputeTool(const std::string& name,MetricSpace *ms,Task *tsk)
{
	// Do not create ComputeTools directly.
	// Instances are created through Grid::CreateTool
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
		// typically constructor should be complete, nothing to do here.
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
}

void ComputeTool::ReadInputFileBlock(std::stringstream& inputString)
{
	std::string com;
	do
	{
		inputString >> com;
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

tw::tool_type ComputeTool::CreateTypeFromInput(const std::vector<std::string>& preamble)
{
	// For creating a named tool at the root level
	std::string type("");
	for (auto word : preamble)
		type = type + word + " ";
	type.pop_back();

	if (type=="eigenmode propagator tool")
		return tw::tool_type::eigenmodePropagator;
	if (type=="adi propagator tool")
		return tw::tool_type::adiPropagator;
	if (type=="isotropic propagator tool")
		return tw::tool_type::isotropicPropagator;
	if (type=="parabolic propagator tool")
		return tw::tool_type::generalParabolicPropagator;
	if (type=="schroedinger propagator tool")
		return tw::tool_type::schroedingerPropagator;
	if (type=="iterative elliptic tool")
		return tw::tool_type::iterativePoissonSolver;
	if (type=="1d elliptic tool")
		return tw::tool_type::ellipticSolver1D;
	if (type=="facr elliptic tool")
		return tw::tool_type::facrPoissonSolver;
	if (type=="eigenmode elliptic tool")
		return tw::tool_type::eigenmodePoissonSolver;
	if (type=="yee propagator tool")
		return tw::tool_type::yeePropagatorPML;
	if (type=="lorentz propagator tool")
		return tw::tool_type::lorentzPropagator;
	if (type=="eos data tool")
		return tw::tool_type::eosData;
	if (type=="eos ideal gas tool")
		return tw::tool_type::eosIdealGas;
	if (type=="eos hot electron tool")
		return tw::tool_type::eosHotElectrons;
	if (type=="eos mix tool")
		return tw::tool_type::eosMixture;
	if (type=="eos ideal gas mix tool")
		return tw::tool_type::eosIdealGasMix;
	if (type=="eos mie gruneisen tool")
		return tw::tool_type::eosMieGruneisen;
	if (type=="eos mie gruneisen2 tool")
		return tw::tool_type::eosMieGruneisen2;
	return tw::tool_type::nullTool;
}

tw::tool_type ComputeTool::CreateTypeFromDirective(std::stringstream& inputString,const std::string& command)
{
	// For creating tools on the fly inside a module block
	std::string word;
	if (command=="elliptic" || command=="elliptical")
	{
		inputString >> word >> word >> word;
		if (word=="1d")
			return tw::tool_type::ellipticSolver1D;
		if (word=="iterative")
			return tw::tool_type::iterativePoissonSolver;
		if (word=="facr")
			return tw::tool_type::facrPoissonSolver;
		if (word=="eigenmode")
			return tw::tool_type::eigenmodePoissonSolver;
	}
	if (command=="propagator")
	{
		inputString >> word >> word;
		if (word=="eigenmode")
			return tw::tool_type::eigenmodePropagator;
		if (word=="adi")
			return tw::tool_type::adiPropagator;
		if (word=="isotropic")
			return tw::tool_type::isotropicPropagator;
		if (word=="parabolic")
			return tw::tool_type::generalParabolicPropagator;
		if (word=="schroedinger")
			return tw::tool_type::schroedingerPropagator;
		if (word=="yee")
			return tw::tool_type::yeePropagatorPML;
		if (word=="lorentz")
			return tw::tool_type::lorentzPropagator;
	}
	if (command=="eos")
	{
		inputString >> word >> word;
		if (word=="ideal-gas")
			return tw::tool_type::eosIdealGas;
		if (word=="hot-electron")
			return tw::tool_type::eosHotElectrons;
		if (word=="mix")
			return tw::tool_type::eosMixture;
		if (word=="ideal-gas-mix")
			return tw::tool_type::eosIdealGasMix;
		if (word=="mie-gruneisen")
			return tw::tool_type::eosMieGruneisen;
		if (word=="mie-gruneisen2")
			return tw::tool_type::eosMieGruneisen2;
	}
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
		case tw::tool_type::eosMieGruneisen:
			ans = new EOSMieGruneisen(name,ms,tsk);
			break;
		case tw::tool_type::eosMieGruneisen2:
			ans = new EOSMieGruneisen2(name,ms,tsk);
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
