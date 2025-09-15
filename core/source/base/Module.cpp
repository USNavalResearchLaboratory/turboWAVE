module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module twmodule;
import input;
import compute_tool;
import region;
import injection;
import diagnostics;
import static_space;
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
	tw::Int steps,req_dim[4],req_dom[4];
	tw::vec4 relativeRef,absoluteRef,spacing;
	bool adaptiveTimestep,adaptiveGrid,found;
	public:
	GridReader(tw::UnitConverter& uc);
	bool Read(const TSTreeCursor *curs,const std::string& src);
	void UpdateTask(Task& tsk);
	void UpdateSpace(MetricSpace& ms);
	tw::Int Steps() { return steps; }
	tw::vec4 GlobalCorner();
	tw::vec4 GlobalSize();
	tw::idx4 GlobalDims() { return tw::idx4(req_dim[0],req_dim[1],req_dim[2],req_dim[3]); }
	tw::grid::geometry Geometry() { return geo; }
	bool FoundGrid() { return found; }
};

export struct Simulation: Task, MetricSpace, tw::input::Visitor
{
	tw::Float unitDensityCGS;
	tw::units nativeUnits; ///< only used as a bridge from input file

	bool interactive,neutralize,movingWindow;
	tw::Int outputLevel,errorCheckingLevel,success_count,failure_count;
	tw::Int previous_timestamp;

	tw::Int dumpPeriod;

	tw::bc::par bc0[4],bc1[4];

	tw::Int inputFilePass;
	tw::input::DirectiveReader outerDirectives;
	std::unique_ptr<GridReader> gridReader;
	std::string raw_src; // raw input file
	std::string src; // expanded input file

	std::vector<Region*> clippingRegion;
	std::vector<SharedTool> tools;
	std::vector<Module*> module;
	// Map of the most recently created module of a given type
	std::map<tw::module_type,Module*> module_map;

	Simulation(const bool interactive,
		const std::string& unitTest,
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

	bool MangleModuleName(std::string& name);
	bool CheckModule(const std::string& name);
	Module* GetModule(const std::string& name);

	bool MangleToolName(std::string& name);
	/// Create a tool that can be shared, returns name in case it had to be mangled
	std::string CreateTool(const std::string& basename,tw::tool_type theType);
	/// Use a previously created tool, if the name was mangled, the mangled name must be passed in
	SharedTool UseTool(const std::string& name);
	/// Remove tool from the master list, only used to conserve memory during testing
	void RemoveTool(const SharedTool& tool);
	void ParseUse(std::vector<SharedTool>& tools,TSTreeCursor *curs,const std::string& src);

	tw::input::navigation visit(TSTreeCursor *curs);
	void InputFileFirstPass();
	void ParseNestedDeclaration(TSTreeCursor *curs,const std::string& src,Module *sup);
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

struct Module:StaticSpace,Testable
{
	std::string name;
	Simulation *owner;
	Module* super;
	std::vector<Module*> submodule;
	std::vector<SharedTool> tools;
	tw::input::DirectiveReader directives;
	tw::UnitConverter native,natural,atomic,plasma,cgs,mks;

	tw::priority updateSequencePriority;
	tw::Int subSequencePriority;
	tw::Int smoothing[4],compensation[4];
	bool suppressNextUpdate;

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
	virtual void Initialize() {;}
	virtual void Reset() {;}
	virtual void Update() {;}
	virtual void MoveWindow() {;}
	virtual void AntiMoveWindow() {;}
	virtual void AdaptGrid() {;}
	virtual tw::Float AdaptTimestep() { return 0.0;}

	void InitializeCLProgram(const std::string& filename);

	/// Carry out checks that are appropriate after the whole input file block has been read in.
	/// This includes searching the SharedTool list for what we need.
	virtual void VerifyInput() {;}
	virtual void ReadInputFileBlock(TSTreeCursor *curs,const std::string& src);
	virtual bool ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual bool ReadQuasitoolBlock(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void StartDiagnostics() {;}
	virtual void Report(Diagnostic&) {;}
	virtual void StatusMessage(std::ostream *dest) {;}

	static std::map<std::string,tw::module_type> Map();
	static bool SingularType(tw::module_type theType);
	static bool AutoModuleType(tw::module_type theType);
	static tw::module_type RequiredSupermoduleType(tw::module_type submoduleType);
	static tw::module_type CreateTypeFromInput(const tw::input::Preamble& preamble);
	/// setup grid, native units, or other environmental factors affecting a test, returns
	/// whether an environment was defined for the given enviro
	static bool SetTestEnvironment(tw::module_type theType,tw::Int enviro,Simulation *sim);
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
	StaticSpace::operator=(*sim);
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

void Module::ReadCheckpoint(std::ifstream& inFile)
{
	StaticSpace::ReadCheckpoint(inFile);
}

void Module::WriteCheckpoint(std::ofstream& outFile)
{
	outFile << name << " ";
	StaticSpace::WriteCheckpoint(outFile);
}

/// @brief called if an assignment was not handled normally
/// @param curs on a directive, use, new, generate, or custom assignment
/// @param src source document
/// @returns whether any directive was handled
bool Module::ReadInputFileDirective(const TSTreeCursor *curs0,const std::string& src)
{
	std::string command = tw::input::node_kind(curs0);
	// Get an existing tool by searching for a name -- ``use <name>``
	if (command=="use") {
		auto curs = tw::input::Cursor(curs0);
		owner->ParseUse(tools,curs.get(),src);
		return true;
	}

	// Take care of nested declarations
	if (command=="new" || command=="generate") {
		auto curs = tw::input::Cursor(curs0);
		owner->ParseNestedDeclaration(curs.get(),src,this);
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

bool Module::ReadQuasitoolBlock(const TSTreeCursor *curs,const std::string& src)
{
	return false;
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
		{"maxwell solver",tw::module_type::directSolver},
		{"curvilinear maxwell solver",tw::module_type::curvilinearDirectSolver},
		{"coulomb gauge solver",tw::module_type::coulombSolver},
		{"far field diagnostic",tw::module_type::farFieldDiagnostic},
		{"electrostatic solver",tw::module_type::electrostatic},
		{"quasistatic solver",tw::module_type::qsLaser},
		{"pgc laser solver",tw::module_type::pgcLaser},
		{"bound solver",tw::module_type::boundElectrons},
		{"schroedinger solver",tw::module_type::schroedinger},
		{"pauli solver",tw::module_type::pauli},
		{"klein gordon solver",tw::module_type::kleinGordon},
		{"dirac solver",tw::module_type::dirac},
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
		return tw::module_type::none;
	}
}

bool Module::SetTestEnvironment(tw::module_type theType,tw::Int enviro,Simulation *sim)
{
	const tw::node4 cyclic = {0,1,1,0};
	const tw::node4 domains = {1,1,1,2};
	const tw::node5 gdim1d = {1,1,1,4,1};
	const tw::node5 gdim2d = {1,4,1,4,1};
	const tw::node5 gdim3d = {1,4,4,4,1};
	const tw::node4 layers = {0,2,2,2};
	sim->Initialize(domains,cyclic);
	switch (theType)
	{
		case tw::module_type::species:
			if (enviro==1) {
				sim->AttachUnits(tw::units::plasma,sim->unitDensityCGS);
				sim->Resize(sim,gdim2d,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.2,0.8),std_packing,layers);
				return true;
			} else if (enviro==2) {
				sim->AttachUnits(tw::units::plasma,sim->unitDensityCGS);
				sim->Resize(sim,gdim3d,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.8,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
		case tw::module_type::pauli:
			return false;
		case tw::module_type::schroedinger:
			if (enviro==1) {
				sim->AttachUnits(tw::units::atomic,sim->unitDensityCGS);
				sim->Resize(sim,gdim3d,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.8,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
		case tw::module_type::dirac:
		case tw::module_type::kleinGordon:
			if (enviro==1) {
				sim->AttachUnits(tw::units::natural,sim->unitDensityCGS);
				sim->Resize(sim,gdim3d,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.8,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
		default:
			if (enviro==1) {
				sim->AttachUnits(tw::units::plasma,sim->unitDensityCGS);
				sim->Resize(sim,gdim3d,tw::vec4(0,0,0,0),tw::vec4(0.1,0.8,0.8,0.8),std_packing,layers);
				return true;
			} else {
				return false;
			}
	}
}
