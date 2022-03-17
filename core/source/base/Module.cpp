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
	suppressNextUpdate = false;
	DiscreteSpace::operator=(*sim);
	// Units
	native = sim->units;
	natural = tw::UnitConverter(tw::units::natural,native);
	atomic = tw::UnitConverter(tw::units::atomic,native);
	plasma = tw::UnitConverter(tw::units::plasma,native);
	cgs = tw::UnitConverter(tw::units::cgs,native);
	mks = tw::UnitConverter(tw::units::mks,native);
	directives.AttachUnits(native);
	// Any Module recognizes smoothing keys
	directives.Add("smoothing",new tw::input::Numbers<tw::Int>(&smoothing[1],3),false);
	directives.Add("compensation",new tw::input::Numbers<tw::Int>(&compensation[1],3),false);
}

Module::~Module()
{
	#ifdef USE_OPENCL

	if (programFilename!="")
		clReleaseProgram(program);

	#endif
}

void Module::AddSubmodule(Module* sub)
{
	submodule.push_back(sub);
	sub->super = this;
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
}

void Module::MoveWindow()
{
	corner.z += spacing.z;
	globalCorner.z += spacing.z;
}

void Module::ReadCheckpoint(std::ifstream& inFile)
{
	DiscreteSpace::ReadCheckpoint(inFile);
	inFile.read((char *)&dt,sizeof(dt));
	inFile.read((char *)&dth,sizeof(dth));
	dti = 1.0/dt;
}

void Module::WriteCheckpoint(std::ofstream& outFile)
{
	outFile << name << " ";
	DiscreteSpace::WriteCheckpoint(outFile);
	outFile.write((char *)&dt,sizeof(dt));
	outFile.write((char *)&dth,sizeof(dth));
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
		Profile *theProfile = dynamic_cast<Profile*>(tool);
		if (theProfile) profile.push_back(theProfile);

		Wave *theWave = dynamic_cast<Wave*>(tool);
		if (theWave) wave.push_back(theWave);

		Conductor *theConductor = dynamic_cast<Conductor*>(tool);
		if (theConductor) conductor.push_back(theConductor);
	}
}

void Module::ReadInputFileBlock(std::stringstream& inputString)
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

bool Module::Test(tw::Int& id)
{
	return false; // did not test, id unchanged means no next test
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

bool Module::AutoModuleType(tw::module_type theType)
{
	return theType==tw::module_type::kinetics || theType==tw::module_type::equilibriumGroup;
}

tw::module_type Module::RequiredSupermoduleType(const tw::module_type submoduleType)
{
	std::map<tw::module_type,tw::module_type> containmentMap =
	{
		{tw::module_type::species,tw::module_type::kinetics},
		{tw::module_type::chemical,tw::module_type::equilibriumGroup},
		{tw::module_type::equilibriumGroup,tw::module_type::sparcHydroManager}
	};
	if (containmentMap.find(submoduleType)==containmentMap.end())
		return tw::module_type::none;
	else
		return containmentMap[submoduleType];
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
		{"far field diagnostic",tw::module_type::farFieldDiagnostic},
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
		{"kinetics",tw::module_type::kinetics},
		{"population diagnostic",tw::module_type::populationDiagnostic}
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
		case tw::module_type::farFieldDiagnostic:
			ans = new FarFieldDiagnostic(name,sim);
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
		case tw::module_type::populationDiagnostic:
			ans = new PopulationDiagnostic(name,sim);
			break;
	}
	return ans;
}

bool Module::SetTestGrid(tw::module_type theType,tw::Int testId,Simulation *sim)
{
	switch (theType)
	{
		case tw::module_type::species:
			if (testId==1)
			{
				sim->Initialize(tw::idx4(1,1,2).array,tw::idx4(4,1,4).array,tw::idx4(1,1,0).array);
				sim->Resize(*sim,tw::vec3(0,0,0),tw::vec3(0.8,0.2,0.8),2);
				sim->SetupTimeInfo(0.1);
				return true;
			} else if (testId==2) {
				sim->Initialize(tw::idx4(1,1,2).array,tw::idx4(4,4,4).array,tw::idx4(1,1,0).array);
				sim->Resize(*sim,tw::vec3(0,0,0),tw::vec3(0.8,0.8,0.8),2);
				sim->SetupTimeInfo(0.1);
				return true;
			}
			return false;
		default:
			if (testId>1)
				return false;
			sim->Initialize(tw::idx4(1,1,2).array,tw::idx4(4,4,4).array,tw::idx4(1,1,0).array);
			sim->Resize(*sim,tw::vec3(0,0,0),tw::vec3(0.8,0.8,0.8),2);
			sim->SetupTimeInfo(0.1);
			return true;
	}
}
