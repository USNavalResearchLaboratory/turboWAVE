module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

module twmodule;
import base;
import preprocess;
import input;
import factory;
import logger;
#include "tw_logger.h"

GridReader::GridReader(tw::UnitConverter& uc)
{
	adaptiveGrid = adaptiveTimestep = found = false;
	geo = tw::grid::cartesian;
	absoluteRef = 0.0;

	directives.Add("steps",new tw::input::Int(&steps),true);
	directives.Add("origin",new tw::input::Numbers<tw::Float>(&relativeRef[0],4));
	directives.Add("shift",new tw::input::Numbers<tw::Float>(&absoluteRef[0],4),false);
	directives.Add("cell size",new tw::input::Numbers<tw::Float>(&spacing[0],4));
	directives.Add("adaptive timestep",new tw::input::Bool(&adaptiveTimestep),false);
	directives.Add("adaptive grid",new tw::input::Bool(&adaptiveGrid),false);
	directives.Add("dimensions",new tw::input::Numbers<tw::Int>(&req_dim[0],4));
	directives.Add("decomposition",new tw::input::Numbers<tw::Int>(&req_dom[0],4));
	std::map<std::string,tw::grid::geometry> geo_map = {{"cartesian",tw::grid::cartesian},{"cylindrical",tw::grid::cylindrical},{"spherical",tw::grid::spherical}};
	directives.Add("geometry",new tw::input::Enums<tw::grid::geometry>(geo_map,&geo),false);
	directives.AttachUnits(uc);
}

bool GridReader::Read(const TSTreeCursor *curs0,const std::string& src)
{
	auto curs = tw::input::Cursor(curs0);
	if (tw::input::node_kind(curs.get())=="new") {
		ts_tree_cursor_goto_first_child(curs.get()); // `new` token
		ts_tree_cursor_goto_next_sibling(curs.get());
		if (tw::input::node_text(curs.get(),src)=="grid") {
			ts_tree_cursor_goto_next_sibling(curs.get());
			if (tw::input::node_kind(curs.get())=="block") {
				ts_tree_cursor_goto_first_child(curs.get());
				tw::input::next_named_node(curs.get(),true);
				directives.ReadAll(curs.get(),src);
				directives.ThrowErrorIfMissingKeys("grid");
				if (req_dom[0]!=1)
					tw::input::ThrowParsingError(curs.get(),src,"Time decomposition must be 1");
				found = true;
				return true;
			} else {
				tw::input::ThrowParsingError(curs.get(),src,"Ill-formed grid block");
			}
		}
	}
	return false;
}

tw::vec4 GridReader::GlobalSize()
{
	return tw::vec4(
		spacing[0] * req_dim[0],
		spacing[1] * req_dim[1],
		spacing[2] * req_dim[2],
		spacing[3] * req_dim[3]
	);
}

tw::vec4 GridReader::GlobalCorner()
{
	return absoluteRef - relativeRef*GlobalSize();
}

void GridReader::UpdateTask(Task& tsk)
{
	for (tw::Int i=0;i<4;i++)
	{
		tsk.domains[i] = req_dom[i];
	}
}

void GridReader::UpdateSpace(MetricSpace& ms)
{
	// whatever not handled by Resize method
	ms.adaptiveGrid = adaptiveGrid;
	ms.adaptiveTimestep = adaptiveTimestep;
}

/// @brief handle retrieval of named tools
/// @param tool tool list owned by a module
/// @param curs should be on a `get` node
/// @param src source document
void Simulation::ToolFromDirective(std::vector<ComputeTool*>& tool,TSTreeCursor *curs,const std::string& src)
{
	logger::DEBUG("get a named tool");
	auto word = tw::input::next_named_node_text(curs,src);
	tw::input::StripQuotes(word);
	if (CheckModule(word))
		tw::input::ThrowParsingError(curs,src,"Tried to <get> module, but <get> can only be used for tools.");
	tool.push_back(GetTool(word,true));
}

/// @brief main visitor callback for AST walker
/// @param curs cursor, typically passed from walker
/// @param inputFilePass 0=syntax, 1=grid, 2=objects
/// @return action for the walker to take upon exit
tw::input::navigation Simulation::visit(TSTreeCursor *curs) {

	if (inputFilePass == 1) {
		std::flush(std::cout);
		if (tw::input::node_kind(curs) == "input_file") {
			return tw::input::navigation::gotoChild;
		} else if (tw::input::node_kind(curs) == "assignment") {
			// Process outer assignments
			outerDirectives.ReadNext(curs,src);
			return tw::input::navigation::gotoSibling;
		} else if (tw::input::node_kind(curs) == "new") {
			// if (!outerDirectives.TestKey("native units")) {
			// 	throw tw::FatalError("`new` encountered before `native units`");
			// }
			if (!outerDirectives.TestKey("unit density")) {
				throw tw::FatalError("`new` encountered before `unit density`");
			}
			ts_tree_cursor_goto_first_child(curs);
			ts_tree_cursor_goto_next_sibling(curs);
			std::string key = tw::input::node_text(curs,src);
			ts_tree_cursor_goto_parent(curs);
			if (key == "grid") {
				if (gridReader->FoundGrid()) {
					tw::input::ThrowParsingError(curs,src,"multiple grids");
				}
				gridReader->Read(curs,src);
				gridReader->UpdateTask(*this);
				gridReader->UpdateSpace(*this);
				return tw::input::navigation::gotoSibling;
			} else if (key == "warp") {
				auto curs1 = tw::input::Cursor(curs);
				tw::input::Preamble preamble = tw::input::GetPreamble(curs1.get(),src);
				MangleToolName(preamble.obj_name);
				std::println(std::cout,"Creating Tool <{}>...",preamble.obj_name);
				// Do not use CreateTool, do not want to increase refCount
				computeTool.push_back(factory::CreateToolFromType(preamble.obj_name,tw::tool_type::warp,this,this));
				computeTool.back()->ReadInputFileBlock(curs1.get(),src);
				computeTool.back()->Initialize(); // OK and necessary to init here
				AttachWarp(dynamic_cast<Warp*>(computeTool.back()));
				return tw::input::navigation::gotoSibling;
			} else {
				// not a grid or warp, so skip it
				return tw::input::navigation::gotoSibling;
			}
		} else {
			// not a basic `new` directive so skip it
			return tw::input::navigation::gotoSibling;
		}


	} else if (inputFilePass == 2) {
		std::flush(std::cout);
		TSNode node = ts_tree_cursor_current_node(curs);
		std::string typ = ts_node_type(node);
		if (tw::input::node_kind(curs) == "input_file") {
			return tw::input::navigation::gotoChild;
		} else if (typ == "comment") {
			return tw::input::navigation::gotoSibling;
		} else if (typ == "assignment") {
			// Process outer directives again, after first pass.
			// This is redundant, but helps minimize keyword collisions with object names.
			outerDirectives.ReadNext(curs,src);
			return tw::input::navigation::gotoSibling;
		} else if (typ == "new" || typ == "associative_new" || typ == "generate") {
			tw::input::Preamble preamble = tw::input::GetPreamble(curs,src);

			// Intercept items processed during the first pass
			if (preamble.obj_key=="grid" || preamble.obj_key=="warp")
				return tw::input::navigation::gotoParentSibling;

			// Intercept regions which are handled specially
			if (preamble.obj_key.substr(0,6)=="region") {
				logger::DEBUG("parsing region");
				std::string rgnType = tw::input::trim(preamble.obj_key.substr(6));
				clippingRegion.push_back(Region::CreateObjectFromString(clippingRegion,rgnType));
				clippingRegion.back()->directives.AttachUnits(units);
				clippingRegion.back()->name = preamble.obj_name;
				clippingRegion.back()->ReadInputFileBlock(curs,src);
				return tw::input::navigation::gotoSibling;
			}
			
			std::string similar_keys;
			if (!factory::VerifyKey(preamble.obj_key,similar_keys)) {
				tw::input::ThrowParsingError(curs, src, std::format("unknown key <{}>",preamble.obj_key),
					similar_keys.size() > 0 ? std::format("<{}> is similar to {}",preamble.obj_key,similar_keys) : "no similar keys found");
			}
			tw::tool_type whichTool = ComputeTool::CreateTypeFromInput(preamble);
			tw::module_type whichModule = Module::CreateTypeFromInput(preamble);

			// Install a pre or post declared tool
			if (whichTool!=tw::tool_type::none) {
				logger::DEBUG("parsing tool");
				if (preamble.attaching) {
					ComputeTool *tool = CreateTool(preamble.obj_name,whichTool);
					std::println(std::cout,"Attaching <{}> to <{}>",tool->name,preamble.owner_name);
					tool->ReadInputFileBlock(curs,src);
					GetModule(preamble.owner_name)->moduleTool.push_back(tool);
				} else {
					bool duplicate = MangleToolName(preamble.obj_name);
					std::println(std::cout,"Creating Tool <{}>...",preamble.obj_name);
					if (duplicate && errorCheckingLevel>0)
						tw::input::ThrowParsingError(curs,src,"duplicate tool name.");
					// Do not use CreateTool, do not want to increase refCount
					computeTool.push_back(factory::CreateToolFromType(preamble.obj_name,whichTool,this,this));
					computeTool.back()->ReadInputFileBlock(curs,src);
				}
				return tw::input::navigation::gotoSibling;
			}

			// Module Installation
			if (whichModule!=tw::module_type::none) {
				logger::DEBUG("parsing module");
				if (Module::SingularType(whichModule))
					if (module_map.find(whichModule)!=module_map.end())
						tw::input::ThrowParsingError(curs,src,"singular module type was created twice.  Check order of input file.");
				MangleModuleName(preamble.obj_name);
				std::println(std::cout,"Installing module <{}>...",preamble.obj_name);
				Module *sub = factory::CreateModuleFromType(preamble.obj_name,whichModule,this);
				module.push_back(sub);
				module_map[whichModule] = sub;
				sub->ReadInputFileBlock(curs,src); // important to note this can change module vector and map if there are nested declarations
				if (preamble.attaching && preamble.owner_name.size() > 0) {
					Module *super = GetModule(preamble.owner_name);
					std::println(std::cout,"Attaching <{}> to <{}>",preamble.obj_name,preamble.owner_name);
					super->AddSubmodule(sub);
				} else if (!preamble.attaching && preamble.owner_name.size() == 0) {
					// If not explicitly attaching, but supermodule is required, find or create one
					tw::module_type reqType = Module::RequiredSupermoduleType(whichModule);
					if (reqType!=tw::module_type::none) {
						Module *super = RecursiveAutoSuper(reqType,preamble.obj_name);
						super->AddSubmodule(sub);
					}
				} else {
					tw::input::ThrowParsingError(curs,src,"supermodule mismatch");
				}
				return tw::input::navigation::gotoSibling;
			}

			throw tw::FatalError("unknown key <" + preamble.obj_key + "> at " + tw::input::loc_str(curs));

		} else if (typ == "reaction" || typ == "collision" || typ == "excitation") {
			for (tw::Int i=0;i<module.size();i++)
				module[i]->ReadQuasitoolBlock(curs,src);
			return tw::input::navigation::gotoSibling;
		} else {
			tw::input::ThrowParsingError(curs,src,"unhandled directive");
			return tw::input::navigation::exit; // compiler doesn't know above line throws
		}


	} else {
		throw tw::FatalError("invalid number of input file passes");
	}
}

/// @brief to be called by supermodules that want to add submodules or tools, can be recursive
/// @param curs should be on the outer directive node
/// @param src source document
/// @param super module that is adding the item
void Simulation::NestedDeclaration(TSTreeCursor *curs,const std::string& src,Module *super)
{
	logger::DEBUG("handling nested declaration");
	tw::input::Preamble preamble = tw::input::GetPreamble(curs,src);
	if (preamble.attaching)
		tw::input::ThrowParsingError(curs,src,"keyword <for> is not allowed in a nested declaration.");

	std::string similar_keys;
	if (!factory::VerifyKey(preamble.obj_key,similar_keys)) {
		tw::input::ThrowParsingError(curs, src, std::format("unknown key <{}>",preamble.obj_key),
			similar_keys.size() > 0 ? std::format("<{}> is similar to {}",preamble.obj_key,similar_keys) : "no similar keys found");
	}

	tw::module_type whichModule = Module::CreateTypeFromInput(preamble);
	tw::tool_type whichTool = ComputeTool::CreateTypeFromInput(preamble);

	if (whichModule!=tw::module_type::none)
	{
		if (Module::SingularType(whichModule))
			if (module_map.find(whichModule)!=module_map.end())
				tw::input::ThrowParsingError(curs,src,"Singular module type was created twice.  Check order of input file.");
		MangleModuleName(preamble.obj_name);
		std::println(std::cout,"   Attaching nested module <{}>...",preamble.obj_name);
		Module *sub = factory::CreateModuleFromType(preamble.obj_name,whichModule,this);
		module.push_back(sub);
		module_map[whichModule] = sub;
		// The following may lead to a recursive call of this function.
		// If so the module vector and map can be modified.
		sub->ReadInputFileBlock(curs,src);
		super->AddSubmodule(sub);
		return;
	}

	if (whichTool!=tw::tool_type::none)
	{
		ComputeTool *tool = CreateTool(preamble.obj_name,whichTool);
		std::println(std::cout,"   Attaching nested tool <{}>...",tool->name);
		tool->ReadInputFileBlock(curs,src);
		super->moduleTool.push_back(tool);
		return;
	}

	tw::input::ThrowParsingError(curs, src, "unprocessed match in nested declaration (region?)");
}

Module* Simulation::RecursiveAutoSuper(tw::module_type reqType,const std::string& basename)
{
	logger::DEBUG("handling automatic supermodule");
	// If it already exists we are done
	if (module_map.find(reqType)!=module_map.end())
			return module_map[reqType];

	// Automatic creation of supermodule
	if (!Module::AutoModuleType(reqType))
		throw tw::FatalError("Module <"+basename+"> requires a supermodule that cannot be created automatically.");
	std::string super_module_name = basename + "_sup";
	MangleModuleName(super_module_name);
	std::println(std::cout,"Installing supermodule triggered by <{}>...",basename);
	Module *super = factory::CreateModuleFromType(super_module_name,reqType,this);
	module.push_back(super);
	module_map[reqType] = super;

	// Handle recursion
	tw::module_type superSuperType = Module::RequiredSupermoduleType(reqType);
	if (superSuperType!=tw::module_type::none)
	{
		Module *superSuper = RecursiveAutoSuper(superSuperType,super_module_name);
		superSuper->AddSubmodule(super);
	}
	return super;
}

/// The first pass through the input file is used to fully initialize the `Task` and
/// `MetricSpace` parent classes.
void Simulation::InputFileFirstPass()
{
	int numRanksProvided,worldRank;
	MPI_Comm_size(MPI_COMM_WORLD,&numRanksProvided);
	MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
	// world rank is suitable for reading task data from restart file
	// because this data is the same in every restart file

	std::stringstream fileName;
	tw::input::FileEnv fenv(inputFileName);
	fenv.OpenDeck(raw_src);
	src = tw::input::Preprocess(&fenv,raw_src);
	// std::print(std::cout,"{}",src);
	auto tree = tw::input::GetTree(src);

	AttachUnits(tw::input::GetNativeUnits(tree,src),tw::input::GetUnitDensityCGS(tree,src));

	outerDirectives.AttachUnits(units);
	outerDirectives.Reset();
	gridReader = new GridReader(units);

	inputFilePass = 1;
	tw::input::WalkTree(tree,this);

	outerDirectives.ThrowErrorIfMissingKeys("Simulation");
	if (!gridReader->FoundGrid())
		throw tw::FatalError("Grid directive was not found.");

	periodic[1] = bc0[1]==tw::bc::par::periodic ? 1 : 0;
	periodic[2] = bc0[2]==tw::bc::par::periodic ? 1 : 0;
	periodic[3] = bc0[3]==tw::bc::par::periodic ? 1 : 0;
	stepsToTake = gridReader->Steps();
	auto gdim_idx4 = gridReader->GlobalDims();
	tw::Int gdim[4] = { 1, gdim_idx4.array[1], gdim_idx4.array[2], gdim_idx4.array[3] };

	// Check integer viability
	int64_t totalCellsPerRank = int64_t(gdim[1])*int64_t(gdim[2])*int64_t(gdim[3])/int64_t(numRanksProvided);
	if (totalCellsPerRank>=pow(2,31) && sizeof(tw::Int)==4)
		throw tw::FatalError("You must recompile turboWAVE with 64 bit integers to handle this many grid cells.");

	// Verify and (if necessary) correct decomposition
	if (NumTasks() != numRanksProvided)
	{
		std::print(std::cout,"{}: Bad decomposition ",term::warning);
		tw::Int ax1=1,ax2=2,ax3=3; // to be sorted so ax1 is longest
		for (tw::Int i=1;i<=3;i++)
		{
			domains[i] = 1;
			if (gdim[i]>=gdim[1] && gdim[i]>=gdim[2] && gdim[i]>=gdim[3])
				ax1 = i;
		}
		for (tw::Int i=1;i<=3;i++)
			ax2 = i==ax1 ? ax2 : i;
		for (tw::Int i=1;i<=3;i++)
			ax3 = i==ax1 || i==ax2 ? ax3 : i;
		if (gdim[ax2]<gdim[ax3])
			std::swap(ax2,ax3);
		domains[ax1] = numRanksProvided;
		while (gdim[ax1]%(domains[ax1]*2)!=0 && domains[ax1]>0)
		{
			domains[ax1] /= 2;
			domains[ax2] *= 2;
		}
		std::println(std::cout,"(defaulting to {}x{}x{})",domains[1],domains[2],domains[3]);
	}
	std::println(std::cout,"{}-Way Decomposition",NumTasks());

	// Set up the domain decomposition
	for (tw::Int i=1;i<=3;i++) {
		if (gdim[i]%domains[i]!=0)
			throw tw::FatalError("global number of cells is not divisible by number of domains along axis");
	}
	Task::Initialize(domains,periodic);
	for (tw::Int i=1;i<=3;i++) {
		if (gdim[i]>1 && gdim[i]/domains[i]%2>0)
			throw tw::FatalError(std::format("local number of cells is not even along non-ignorable axis {}",i));
	}
	Resize(this,gdim,gridReader->GlobalCorner(),gridReader->GlobalSize(),2,gridReader->Geometry());
	delete gridReader;
	ts_tree_delete(tree);

	// Random numbers
	uniformDeviate = new UniformDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided));
	gaussianDeviate = new GaussianDeviate(1 + strip[0].Get_rank()*(MaxSeed()/numRanksProvided) + MaxSeed()/(2*numRanksProvided));
}

void Simulation::ReadInputFile()
{
	std::println(std::cout,"Reading Input File...\n");
	std::flush(std::cout);

	outerDirectives.Reset();
	auto tree = tw::input::GetTree(src);
	inputFilePass = 2;
	tw::input::WalkTree(tree,this);
	std::flush(std::cout);
	ts_tree_delete(tree);
	outerDirectives.ThrowErrorIfMissingKeys("Simulation");
	if (movingWindow && !outerDirectives.TestKey("window speed")) {
		solutionVelocity[3] = 1.0;
	}
}

