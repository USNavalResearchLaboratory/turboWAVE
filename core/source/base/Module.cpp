module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module twmodule;
import input;
import compute_tool;
import region;
import injection;
import diagnostics;
import discrete_space;
import metric_space;
import tw_iterator;

namespace tw
{
	export enum class module_type {	none,
					electrostatic,
					coulombSolver,directSolver,curvilinearDirectSolver,farFieldDiagnostic,
					qsLaser,pgcLaser,
					kinetics,species,fluidFields,equilibriumGroup,chemical,sparcHydroManager,
					boundElectrons,schroedinger,pauli,kleinGordon,dirac,populationDiagnostic};
	export enum class priority { diagnostic, source, field };
	inline std::map<tw::priority,tw::Int> priority_sort_map()
	{
		std::map<tw::priority,tw::Int> ans = {{priority::diagnostic,100},{priority::source,200},{priority::field,300}};
		return ans;
	}
}

export struct Module;

/// This class reads the grid block from the input file.  The grid block
/// contains data pertinent to both `Task` and `MetricSpace`, hence the
/// need for a broker class.
class GridReader
{
	tw::input::DirectiveReader directives;
	tw::grid::geometry geo;
	tw::Int req_dim[4],req_dom[4];
	tw::vec4 relativeRef,absoluteRef,spacing;
	bool adaptiveTimestep,adaptiveGrid,found;
	public:
	GridReader(tw::UnitConverter& uc);
	bool Read(const TSTreeCursor *curs,const std::string& src);
	void UpdateTask(Task& tsk);
	void UpdateSpace(MetricSpace& ms);
	tw::vec4 GlobalCorner();
	tw::vec4 GlobalSize();
	tw::idx4 GlobalDims() { return tw::idx4(req_dim[0],req_dim[1],req_dim[2],req_dim[3]); }
	tw::grid::geometry Geometry() { return geo; }
	bool FoundGrid() { return found; }
};

export struct Simulation: Task, MetricSpace, tw::input::Visitor
{
	tw::Float dtMin,dtMax,dtCritical,elapsedTime,elapsedTimeMax;
	tw::Float signalPosition,windowPosition,signalSpeed;
	tw::Float antiSignalPosition,antiWindowPosition;
	tw::Float unitDensityCGS;
	tw::units nativeUnits;

	bool neutralize,movingWindow;
	tw::Int outputLevel,errorCheckingLevel,success_count,failure_count;
	tw::Int stepNow;
	tw::Int lastTime;

	tw::Int dumpPeriod;

	tw::bc::par bc0[4],bc1[4];

	tw::Int inputFilePass;
	tw::input::DirectiveReader outerDirectives;
	GridReader *gridReader;
	std::string raw_src; // raw input file
	std::string src; // expanded input file

	std::vector<Region*> clippingRegion;
	std::vector<ComputeTool*> computeTool;
	std::vector<Module*> module;
	// Map of the most recently created module of a given type
	std::map<tw::module_type,Module*> module_map;

	Simulation(const std::string& unitTest,
		const std::string& inputFileName,
		const std::string& restartFileName,
		const std::string& platform,
		const std::string& device,
		const tw::Int& outputLevel,
		const tw::Int& errorCheckingLevel);
	virtual ~Simulation();
	void SetupIO();
	virtual void Run();
	virtual void Test();
	void PrepareSimulation();
	void FundamentalCycle();
	void MoveWindow();
	void AntiMoveWindow();
	void Diagnose();

	void UpdateTimestep(tw::Float dt0);
	bool IsFirstStep() {
		return stepNow == 0;
	}
	bool IsLastStep() {
		// Actual steps taken = stepsToTake+1
		// This allows us to write both the initial condition and the last available data.
		// If there is a restart the initial step is stepsToTake+1 and 1 less step is taken.
		return stepNow == dim[0];
	}
	Evolution GetEvo() {
		return Evolution {
			.dtMin = dtMin,
			.dtMax = dtMax,
			.dtCritical = dtCritical,
			.elapsedTime = elapsedTime,
			.elapsedTimeMax = elapsedTimeMax,
			.signalPosition = signalPosition,
			.windowPosition = windowPosition,
			.signalSpeed = signalSpeed,
			.antiSignalPosition = antiSignalPosition,
			.antiWindowPosition = antiWindowPosition
		};
	}

	bool MangleModuleName(std::string& name);
	bool CheckModule(const std::string& name);
	Module* GetModule(const std::string& name);

	bool MangleToolName(std::string& name);
	ComputeTool* CreateTool(const std::string& basename,tw::tool_type theType);
	ComputeTool* GetTool(const std::string& name,bool attaching);
	void ToolFromDirective(std::vector<ComputeTool*>& tool,TSTreeCursor *curs,const std::string& src);
	bool RemoveTool(ComputeTool *theTool);

	tw::input::navigation visit(TSTreeCursor *curs);
	void InputFileFirstPass();
	void NestedDeclaration(TSTreeCursor *curs,const std::string& src,Module *sup);
	Module* RecursiveAutoSuper(tw::module_type reqType,const std::string& basename);
	void ReadInputFile();
	void ReadCheckpoint(std::ifstream& inFile);
	void WriteCheckpoint(std::ofstream& outFile);
	void InteractiveCommand(const std::string& cmd,std::ostream *theStream);

	#ifdef USE_OPENCL
	void PrintGPUInformation();
	void CellUpdateProtocol(cl_kernel k)
	{
		MetricSpace::CellUpdateProtocol(k,commandQueue);
	}
	void ElementUpdateProtocol(cl_kernel k)
	{
		MetricSpace::ElementUpdateProtocol(k,commandQueue);
	}
	void LocalUpdateProtocol(cl_kernel k)
	{
		MetricSpace::LocalUpdateProtocol(k,commandQueue);
	}
	void PointUpdateProtocol(cl_kernel k)
	{
		MetricSpace::PointUpdateProtocol(k,commandQueue);
	}
	void StripUpdateProtocol(cl_kernel k,tw::Int axis,tw::Int stripArgument)
	{
		MetricSpace::StripUpdateProtocol(k,commandQueue,axis,stripArgument);
	}
	#endif


};

struct Module:DiscreteSpace
{
	std::string name;
	Simulation *owner;
	Module* super;
	std::vector<Module*> submodule;
	std::vector<ComputeTool*> moduleTool;
	tw::input::DirectiveReader directives;
	tw::UnitConverter native,natural,atomic,plasma,cgs,mks;

	tw::priority updateSequencePriority;
	tw::Int subSequencePriority;
	tw::Int smoothing[4],compensation[4];
	bool suppressNextUpdate;

	// Strongly typed ComputeTool lists that are provided for free
	std::vector<Profile*> profile;
	std::vector<Wave*> wave;
	std::vector<Conductor*> conductor;

	// OpenCL Support
	private:
	std::string programFilename;
	std::string buildLog;
	public:
	#ifdef USE_OPENCL
	cl_program program;
	#endif

	Module(const std::string& name,Simulation *sim);
	virtual ~Module();
	void AddSubmodule(Module* sub);
	virtual void PublishResource(void* resource,const std::string& description);
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void ExchangeResources() {;}
	virtual void Initialize();
	virtual void Reset() {;}
	virtual void Update() {;}
	virtual void MoveWindow();
	virtual void AntiMoveWindow() {;}
	virtual void AdaptGrid() {;}
	virtual tw::Float AdaptTimestep() { return 0.0;}

	void InitializeCLProgram(const std::string& filename);

	virtual void VerifyInput();
	virtual void ReadInputFileBlock(TSTreeCursor *curs,const std::string& src);
	virtual bool ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual bool ReadQuasitoolBlock(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void StartDiagnostics();
	virtual void Report(Diagnostic&);
	virtual void WarningMessage();
	virtual void StatusMessage(std::ostream *dest) {;}
	virtual bool Test(tw::Int& id);

	static std::map<std::string,tw::module_type> Map();
	static bool SingularType(tw::module_type theType);
	static bool AutoModuleType(tw::module_type theType);
	static tw::module_type RequiredSupermoduleType(tw::module_type submoduleType);
	static tw::module_type CreateTypeFromInput(const tw::input::Preamble& preamble);
	static bool SetTestGrid(tw::module_type theType,tw::Int gridId,Simulation *sim);
};

export struct ModuleComparator
{
	bool operator() (Module* const& m1,Module* const& m2)
	{
		tw::Int idx1 = tw::priority_sort_map()[m1->updateSequencePriority] + m1->subSequencePriority;
		tw::Int idx2 = tw::priority_sort_map()[m2->updateSequencePriority] + m2->subSequencePriority;
		return idx1 < idx2;
	}
};

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

void Module::WarningMessage()
{
	if (buildLog.size()>4)
	{
		std::println(std::cout,"{}: Build log for {} is not empty:",term::warning,programFilename);
		std::println(std::cout,"{}",buildLog);
	}
	if (owner->movingWindow)
		for (tw::Int i=0;i<profile.size();i++)
			if (profile[i]->theRgn->moveWithWindow==true && profile[i]->theRgn->name!="entire")
				std::println(std::cout,"{}: region {} in motion in module {}",term::warning,profile[i]->theRgn->name,name);
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
		{"maxwell",tw::module_type::directSolver},
		{"curvilinear maxwell",tw::module_type::curvilinearDirectSolver},
		{"coulomb gauge",tw::module_type::coulombSolver},
		{"far field diagnostic",tw::module_type::farFieldDiagnostic},
		{"electrostatic",tw::module_type::electrostatic},
		{"quasistatic",tw::module_type::qsLaser},
		{"pgc",tw::module_type::pgcLaser},
		{"bound",tw::module_type::boundElectrons},
		{"schroedinger",tw::module_type::schroedinger},
		{"pauli",tw::module_type::pauli},
		{"klein gordon",tw::module_type::kleinGordon},
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

/// @brief match input file key after normalizing white space
/// @param preamble data extracted while parsing object
/// @return the corresponding module_type enumeration, may return tw::module_type::none
tw::module_type Module::CreateTypeFromInput(const tw::input::Preamble& preamble)
{
	std::map<std::string,tw::module_type> module_map = Module::Map();
	std::stringstream raw_key(preamble.obj_key);
	std::string normalized,word;
	do {
		raw_key >> word;
		if (!normalized.empty()) {
			normalized += " ";
		}
		normalized += word;
	} while (!raw_key.eof());
	if (module_map.find(normalized)!=module_map.end()) {
		return module_map[normalized];
	} else {
		std::string messg;
		for (const auto& pair : module_map) {
			tw::input::BuildSimilar(messg,normalized,pair.first);
		}
		if (messg.size() > 0) {
			std::println(std::cout,"{}INFO{}: <{}> is similar to {}",term::cyan,term::reset_all,normalized,messg);
		}
		return tw::module_type::none;
	}
}

bool Module::SetTestGrid(tw::module_type theType,tw::Int gridId,Simulation *sim)
{
	const tw::Int cyclic[4] = {0,1,1,0};
	const tw::Int domains[4] = {1,1,1,2};
	const tw::Int gdim1d[4] = {1,1,1,4};
	const tw::Int gdim2d[4] = {1,4,1,4};
	const tw::Int gdim3d[4] = {1,4,4,4};
	sim->Initialize(domains,cyclic);
	switch (theType)
	{
		case tw::module_type::species:
			if (gridId==1)
			{
				sim->Resize(sim,gdim2d,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.2,0.8),2);
				return true;
			} else if (gridId==2) {
				sim->Resize(sim,gdim3d,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.8,0.8),2);
				return true;
			}
			return false;
		default:
			if (gridId>1)
				return false;
			sim->Resize(sim,gdim3d,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.8,0.8),2);
			return true;
	}
}
