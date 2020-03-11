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
	updateSequencePriority = tw::priority::source;
	subSequencePriority = 0;
	for (tw::Int i=0;i<4;i++)
		smoothing[i] = compensation[i] = 0;
	typeCode = tw::module_type::none;
	suppressNextUpdate = false;
	DiscreteSpace::operator=(*sim);
	// Any Module recognizes smoothing keys
	directives.Add("smoothing",new tw::input::Numbers<tw::Int>(&smoothing[1],3));
	directives.Add("compensation",new tw::input::Numbers<tw::Int>(&compensation[1],3));
}

Module::~Module()
{
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

void Module::ReadCheckpoint(std::ifstream& inFile)
{
	DiscreteSpace::ReadCheckpoint(inFile);
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

void Module::WriteCheckpoint(std::ofstream& outFile)
{
	outFile.write((char *)&typeCode,sizeof(typeCode));
	outFile << name << " ";
	DiscreteSpace::WriteCheckpoint(outFile);
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

	// The following captures tools that are handled automatically.
	// N.b. this does not mean they are automatically attached to all modules, only that
	// we automatically populate a strongly typed list.

	for (auto tool : moduleTool)
	{
		Profile *theProfile;
		theProfile = dynamic_cast<Profile*>(tool);
		if (theProfile!=NULL)
			profile.push_back(theProfile);

		Wave *theWave;
		theWave = dynamic_cast<Wave*>(tool);
		if (theWave!=NULL)
			wave.push_back(theWave);

		Conductor *theConductor;
		theConductor = dynamic_cast<Conductor*>(tool);
		if (theConductor!=NULL)
			conductor.push_back(theConductor);
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
}

bool Module::ReadQuasitoolBlock(const tw::input::Preamble& preamble,std::stringstream& inputString)
{
	return false;
}

void Module::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	// Handle whatever the DirectiveReader object did not.
	// This includes keywords "get" and "new"

	// Get an existing tool by searching for a name -- ``get = <name>``
	if (command=="get")
		owner->ToolFromDirective(moduleTool,inputString,command);

	// Take care of nested declarations
	if (command=="new" || command=="generate")
		owner->NestedDeclaration(command,inputString,this);
}

void Module::StartDiagnostics()
{
}

void Module::Report(Diagnostic& diagnostic)
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
	return tw::module_type::none;
}

bool Module::QuasitoolNeedsModule(const tw::input::Preamble& preamble)
{
	bool ans = false;
	ans = ans || (preamble.words[0]=="reaction");
	ans = ans || (preamble.words[0]=="collision");
	ans = ans || (preamble.words[0]=="excitation");
	return ans;
}

std::map<std::string,tw::module_type> Module::Map()
{
	return
	{
		{"direct",tw::module_type::directSolver},
		{"curvilinear direct",tw::module_type::curvilinearDirectSolver},
		{"coulomb",tw::module_type::coulombSolver},
		{"electrostatic",tw::module_type::electrostatic},
		{"quasistatic",tw::module_type::qsLaser},
		{"pgc",tw::module_type::pgcLaser},
		{"bound",tw::module_type::boundElectrons},
		{"schroedinger",tw::module_type::schroedinger},
		{"pauli",tw::module_type::pauli},
		{"klein",tw::module_type::kleinGordon},
		{"dirac",tw::module_type::dirac},
		{"fluid",tw::module_type::fluidFields},
		{"hydrodynamics",tw::module_type::sparcHydroManager},
		{"group",tw::module_type::equilibriumGroup},
		{"chemical",tw::module_type::chemical},
		{"species",tw::module_type::species},
		{"kinetics",tw::module_type::kinetics}
	};
}

tw::module_type Module::CreateTypeFromInput(const tw::input::Preamble& preamble)
{
	// Look for a Module key on a preamble (words between new and opening brace) and return the type of module.
	const tw::Int max_words = preamble.words.size();
	std::map<std::string,tw::module_type> module_map = Module::Map();
	for (tw::Int i=1;i<=max_words;i++)
	{
		std::string key(tw::input::GetPhrase(preamble.words,i));
		if (module_map.find(key)!=module_map.end())
			return module_map[key];
	}
	return tw::module_type::none;
}

Module* Module::CreateObjectFromType(const std::string& name,tw::module_type theType,Simulation* sim)
{
	Module *ans;
	switch (theType)
	{
		case tw::module_type::none:
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
	ans->ReadCheckpoint(inFile);
	return ans;
}
