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
	corner[3] += spacing[3];
	globalCorner[3] += spacing[3];
}

void Module::ReadCheckpoint(std::ifstream& inFile)
{
	DiscreteSpace::ReadCheckpoint(inFile);
	inFile.read((char *)&spacing,sizeof(spacing));
	inFile.read((char *)&freq,sizeof(freq));
}

void Module::WriteCheckpoint(std::ofstream& outFile)
{
	outFile << name << " ";
	DiscreteSpace::WriteCheckpoint(outFile);
	outFile.write((char *)&spacing,sizeof(spacing));
	outFile.write((char *)&freq,sizeof(freq));
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

/// @brief called if an assignment was not handled normally
/// @param curs on a directive, get, new, generate, or custom assignment
/// @param src source document
/// @returns whether any directive was handled
bool Module::ReadInputFileDirective(const TSTreeCursor *curs0,const std::string& src)
{
	std::string command = tw::input::node_kind(curs0);
	// Get an existing tool by searching for a name -- ``get = <name>``
	if (command=="get") {
		TSTreeCursor curs = ts_tree_cursor_copy(curs0);
		owner->ToolFromDirective(moduleTool,&curs,src);
		return true;
	}

	// Take care of nested declarations
	if (command=="new" || command=="generate") {
		TSTreeCursor curs = ts_tree_cursor_copy(curs0);
		owner->NestedDeclaration(&curs,src,this);
		return true;
	}

	return false;
}

/// @brief read all directives in the block
/// @param curs can be on block or on first child of block
/// @param src source document
void Module::ReadInputFileBlock(TSTreeCursor *curs,const std::string& src)
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

bool Module::ReadQuasitoolBlock(const TSTreeCursor *curs,const std::string& src)
{
	return false;
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

/// @brief weak matching of the input file key (legacy compatibility)
/// @param preamble data extracted while parsing object
/// @return the corresponding module_type enumeration
tw::module_type Module::CreateTypeFromInput(const tw::input::Preamble& preamble)
{
	// strategy is to lop off trailing words until we get a match
	std::map<std::string,tw::module_type> module_map = Module::Map();
	std::string key_now = preamble.obj_key;
	while (true) {
		if (module_map.find(key_now)!=module_map.end()) {
			return module_map[key_now];
		} else {
			auto p = key_now.rfind(' ');
			if (p != std::string::npos) {
				key_now = key_now.substr(0,p);
			} else {
				break;
			}
		}
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

bool Module::SetTestGrid(tw::module_type theType,tw::Int gridId,Simulation *sim)
{
	switch (theType)
	{
		case tw::module_type::species:
			if (gridId==1)
			{
				sim->Initialize(tw::idx4(1,1,1,2).array,tw::idx4(1,4,1,4).array,tw::idx4(0,1,1,0).array);
				sim->Resize(*sim,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.2,0.8),2);
				return true;
			} else if (gridId==2) {
				sim->Initialize(tw::idx4(1,1,1,2).array,tw::idx4(1,4,4,4).array,tw::idx4(0,1,1,0).array);
				sim->Resize(*sim,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.8,0.8),2);
				return true;
			}
			return false;
		default:
			if (gridId>1)
				return false;
			sim->Initialize(tw::idx4(1,1,1,2).array,tw::idx4(1,4,4,4).array,tw::idx4(0,1,1,0).array);
			sim->Resize(*sim,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.8,0.8),2);
			return true;
	}
}
