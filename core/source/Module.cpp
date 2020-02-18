#include "simulation.h"
#include "fieldSolve.h"
#include "electrostatic.h"
#include "laserSolve.h"
#include "fluid.h"
#include "quantum.h"
#include "particles.h"
#include "solidState.h"


////////////////////////
//                    //
//    MODULE CLASS    //
//                    //
////////////////////////


Module::Module(const std::string& name,Simulation* sim)
{
	this->name = name;
	owner = sim;
	super = NULL;
	programFilename = "";
	buildLog = "";
	updateSequencePriority = 1;
	for (tw::Int i=0;i<4;i++)
		smoothing[i] = compensation[i] = 0;
	typeCode = tw::module_type::nullModule;
	suppressNextUpdate = false;
	DiscreteSpace::operator=(*sim);
	// Any Module will recognize any on-the-fly ComputeTool key
	// for (auto tool_key : ComputeTool::Map())
	// 	directives.Add(tool_key.first,new tw::input::Custom);
	// Any Module recognizes smoothing keys
	directives.Add("smoothing",new tw::input::Numbers<tw::Int>(&smoothing[1],3));
	directives.Add("compensation",new tw::input::Numbers<tw::Int>(&compensation[1],3));
}

Module::~Module()
{
	for (tw::Int i=0;i<profile.size();i++)
		delete profile[i];

	#ifdef USE_OPENCL

	if (programFilename!="")
		clReleaseProgram(program);

	#endif
}

bool Module::ValidSubmodule(Module* sub)
{
	return false;
}

bool Module::AddSubmodule(Module* sub)
{
	if (ValidSubmodule(sub))
	{
		submodule.push_back(sub);
		sub->super = this;
		return true;
	}
	else
		return false;
}

void Module::InitializeCLProgram(const std::string& filename)
{
	#ifdef USE_OPENCL
	owner->InitializeCLProgram(program,filename,buildLog);
	programFilename = filename;
	#endif
}

void Module::PublishResource(void* resource,const std::string& description)
{
	tw::Int i;
	for (i=0;i<owner->module.size();i++)
		owner->module[i]->InspectResource(resource,description);
}

bool Module::InspectResource(void* resource,const std::string& description)
{
	return false;
}

void Module::Initialize()
{
	if (!owner->restarted)
		for (auto p : profile)
			p->Initialize();
}

void Module::ReadData(std::ifstream& inFile)
{
	DiscreteSpace::ReadData(inFile);
	inFile.read((char *)&dt,sizeof(dt));
	inFile.read((char *)&dth,sizeof(dth));
	inFile.read((char *)&smoothing[0],sizeof(smoothing));
	inFile.read((char *)&compensation[0],sizeof(compensation));
	dti = 1.0/dt;

	// Populate the ComputeTool list
	// This relies on all ComputeTools being loaded first.
	tw::Int num;
	inFile.read((char *)&num,sizeof(num));
	for (tw::Int i=0;i<num;i++)
		moduleTool.push_back(owner->GetRestartedTool(inFile));

	// Find supermodule and setup containment
	std::string super_name;
	inFile >> super_name;
	inFile.ignore();
	if (super_name!="NULL")
		owner->GetModule(super_name)->AddSubmodule(this);
}

void Module::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&typeCode,sizeof(typeCode));
	outFile << name << " ";
	DiscreteSpace::WriteData(outFile);
	outFile.write((char *)&dt,sizeof(dt));
	outFile.write((char *)&dth,sizeof(dth));
	outFile.write((char *)&smoothing[0],sizeof(smoothing));
	outFile.write((char *)&compensation[0],sizeof(compensation));

	// Save the ComputeTool list
	tw::Int num = moduleTool.size();
	outFile.write((char *)&num,sizeof(num));
	for (tw::Int i=0;i<num;i++)
		moduleTool[i]->SaveToolReference(outFile);

	// Write the name of the supermodule
	// This can be used to reconstruct the whole hierarchy, assuming modules are sorted correctly
	if (super==NULL)
		outFile << "NULL";
	else
		outFile << super->name;
	outFile << " ";
}

void Module::VerifyInput()
{
	// Carry out checks that are appropriate after the whole input file block has been read in.
	// This includes searching the list of ComputTool pointers for what we need.

	// Profiles are an automatically managed ComputeTool.  There are two aspects of this:
	// 1. The profile list has to be constructed.
	// 2. The named region has to be found.  This is mostly for Profile, but we do it for all ComputeTools while we're at it.

	for (auto tool : moduleTool)
	{
		Profile *theProfile;
		theProfile = dynamic_cast<Profile*>(tool);
		if (theProfile!=NULL)
			profile.push_back(theProfile);

		if (tool->region_name=="tw::entire")
			tool->theRgn = owner->clippingRegion[0];
		else
			tool->theRgn = Region::FindRegion(owner->clippingRegion,tool->region_name);
	}
}

void Module::ReadInputFileBlock(std::stringstream& inputString)
{
	std::string com;
	do
	{
		com = directives.ReadNext(inputString);
		ReadInputFileDirective(inputString,com);
	} while (com!="}");
	VerifyInput();
}

bool Module::ReadQuasitoolBlock(const std::vector<std::string>& preamble,std::stringstream& inputString)
{
	return false;
}

void Module::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	// Handle whatever the DirectiveReader object did not.
	// This includes keywords "get" and "new"

	// Get an existing tool by searching for a name -- ``get tool with name = <name>``
	if (command=="get")
		owner->ToolFromDirective(moduleTool,inputString,command);

	// Read in submodules that are explicitly enclosed in this module's block
	if (command=="new")
		owner->ReadSubmoduleBlock(inputString,this);
}

void Module::StartDiagnostics()
{
}

void Module::WarningMessage(std::ostream *theStream)
{
	if (buildLog.size()>4)
	{
		*theStream << "WARNING : Build log for " << programFilename << " is not empty:" << std::endl;
		*theStream << buildLog << std::endl;
	}
	if (owner->movingWindow)
		for (tw::Int i=0;i<profile.size();i++)
			if (profile[i]->theRgn->moveWithWindow==true && profile[i]->theRgn->name!="entire")
				*theStream << "WARNING : region " << profile[i]->theRgn->name << " in motion in module " << name << std::endl;
}

////// STATIC MEMBERS FOLLOW

bool Module::SingularType(tw::module_type theType)
{
	return theType==tw::module_type::kinetics || theType==tw::module_type::sparcHydroManager;
}

tw::module_type Module::CreateSupermoduleTypeFromSubmoduleKey(const std::string& key)
{
	if (key=="species")
		return tw::module_type::kinetics;
	if (key=="chemical")
		return tw::module_type::equilibriumGroup;
	if (key=="group")
		return tw::module_type::sparcHydroManager;
	return tw::module_type::nullModule;
}

bool Module::QuasitoolNeedsModule(const std::vector<std::string>& preamble)
{
	bool ans = false;
	ans = ans || (preamble[0]=="phase");
	ans = ans || (preamble[0]=="orbit");
	ans = ans || (preamble[0]=="detector");
	ans = ans || (preamble[0]=="reaction");
	ans = ans || (preamble[0]=="collision");
	ans = ans || (preamble[0]=="excitation");
	return ans;
}

tw::module_type Module::CreateTypeFromInput(const std::vector<std::string>& preamble)
{
	if (preamble[0]=="curvilinear" && preamble[1]=="direct")
		return tw::module_type::curvilinearDirectSolver;
	if (preamble[0]=="coulomb")
		return tw::module_type::coulombSolver;
	if (preamble[0]=="electrostatic")
		return tw::module_type::electrostatic;
	if (preamble[0]=="direct")
		return tw::module_type::directSolver;
	if (preamble[0]=="quasistatic")
		return tw::module_type::qsLaser;
	if (preamble[0]=="pgc" || preamble[0]=="PGC")
		return tw::module_type::pgcLaser;
	if (preamble[0]=="bound")
		return tw::module_type::boundElectrons;
	if (preamble[0]=="schroedinger" || preamble[0]=="atomic")
		return tw::module_type::schroedinger;
	if (preamble[0]=="pauli")
		return tw::module_type::pauli;
	if (preamble[0]=="klein")
		return tw::module_type::kleinGordon;
	if (preamble[0]=="dirac")
		return tw::module_type::dirac;
	if (preamble[0]=="fluid")
		return tw::module_type::fluidFields;
	if (preamble[0]=="chemistry" || preamble[0]=="hydro" || preamble[0]=="hydrodynamics")
		return tw::module_type::sparcHydroManager;
	if (preamble[0]=="group")
		return tw::module_type::equilibriumGroup;
	if (preamble[0]=="chemical")
		return tw::module_type::chemical;
	if (preamble[0]=="species")
		return tw::module_type::species;
	if (preamble[0]=="kinetics")
		return tw::module_type::kinetics;
	return tw::module_type::nullModule;
}

Module* Module::CreateObjectFromType(const std::string& name,tw::module_type theType,Simulation* sim)
{
	Module *ans;
	switch (theType)
	{
		case tw::module_type::nullModule:
			ans = NULL;
			break;
		case tw::module_type::curvilinearDirectSolver:
			ans = new CurvilinearDirectSolver(name,sim);
			break;
		case tw::module_type::electrostatic:
			ans = new Electrostatic(name,sim);
			break;
		case tw::module_type::coulombSolver:
			ans = new CoulombSolver(name,sim);
			break;
		case tw::module_type::directSolver:
			ans = new DirectSolver(name,sim);
			break;
		case tw::module_type::qsLaser:
			ans = new QSSolver(name,sim);
			break;
		case tw::module_type::pgcLaser:
			ans = new PGCSolver(name,sim);
			break;
		case tw::module_type::boundElectrons:
			ans = new BoundElectrons(name,sim);
			break;
		case tw::module_type::fluidFields:
			ans = new Fluid(name,sim);
			break;
		case tw::module_type::equilibriumGroup:
			ans = new EquilibriumGroup(name,sim);
			break;
		case tw::module_type::chemical:
			ans = new Chemical(name,sim);
			break;
		case tw::module_type::sparcHydroManager:
			ans = new sparc::HydroManager(name,sim);
			break;
		case tw::module_type::species:
			ans = new Species(name,sim);
			break;
		case tw::module_type::kinetics:
			ans = new Kinetics(name,sim);
			break;
		case tw::module_type::schroedinger:
			ans = new Schroedinger(name,sim);
			break;
		case tw::module_type::pauli:
			ans = new Pauli(name,sim);
			break;
		case tw::module_type::kleinGordon:
			ans = new KleinGordon(name,sim);
			break;
		case tw::module_type::dirac:
			ans = new Dirac(name,sim);
			break;
	}
	return ans;
}

Module* Module::CreateObjectFromFile(std::ifstream& inFile,Simulation* sim)
{
	// This might work, but the containment hierarchy better be sorted.
	tw::module_type theType;
	std::string name;
	Module *ans;
	inFile.read((char*)&theType,sizeof(tw::module_type));
	inFile >> name;
	inFile.ignore();
	ans = CreateObjectFromType(name,theType,sim);
	ans->ReadData(inFile);
	return ans;
}
