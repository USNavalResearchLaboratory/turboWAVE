module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"

module simulation;
import base;
import static_space;
import preprocess;
import input;
import factory;
import logger;

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
			try {
				gridReader->ReadOuter(curs,src);
			} catch (tw::FatalError& e) {
				if (std::strcmp(e.what(),"unknown key")) {
					outerDirectives.ReadNext(curs,src);
				} else {
					throw e;
				}
			}
			return tw::input::navigation::gotoSibling;
		} else if (tw::input::node_kind(curs) == "new") {
			ts_tree_cursor_goto_first_child(curs);
			ts_tree_cursor_goto_next_sibling(curs);
			std::string key = tw::input::node_text(curs,src);
			ts_tree_cursor_goto_parent(curs);
			if (key == "grid") {
				if (gridReader->FoundGrid()) {
					tw::input::ThrowParsingError(curs,src,"multiple grids");
				}
				gridReader->ReadGridBlock(curs,src);
				gridReader->UpdateTask(*task);
				gridReader->UpdateSpace(*space);
				return tw::input::navigation::gotoSibling;
			} else if (key == "warp") {
				auto curs1 = tw::input::Cursor(curs);
				tw::input::Preamble preamble = tw::input::GetPreamble(curs1.get(),src);
				logger::INFO(std::format("Creating Tool <{}>...",preamble.obj_name));
				auto new_tool = CreateTool(preamble.obj_name, tw::tool_type::warp);
				AddTool(new_tool);
				new_tool->ReadInputFileBlock(curs1.get(),src);
				new_tool->Initialize(); // OK and necessary to init here
				space->AttachWarp(std::dynamic_pointer_cast<warp_base>(new_tool));
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
			// We don't need the values, but must process them for parsing reasons.
			try {
				gridReader->ReadOuter(curs,src);
			} catch (tw::FatalError& e) {
				if (std::strcmp(e.what(),"unknown key")) {
					outerDirectives.ReadNext(curs,src);
				} else {
					throw e;
				}
			}
			return tw::input::navigation::gotoSibling;
		} else if (typ == "new" || typ == "associative_new" || typ == "generate") {
			tw::input::Preamble preamble = tw::input::GetPreamble(curs,src);

			// Intercept items processed during the first pass
			if (preamble.obj_key=="grid" || preamble.obj_key=="warp")
				return tw::input::navigation::gotoParentSibling;

			// TODO: deleted region processing, now we have to make sure regions are found as tools
			
			std::string similar_keys;
			if (!factory::VerifyKey(preamble.obj_key,similar_keys)) {
				tw::input::ThrowParsingError(curs, src, std::format("unknown key <{}>",preamble.obj_key),
					similar_keys.size() > 0 ? std::format("<{}> is similar to {}",preamble.obj_key,similar_keys) : "no similar keys found");
			}
			auto whichTool = ComputeTool::CreateTypeFromInput(preamble);

			// Install a pre or post declared tool
			if (whichTool!=tw::tool_type::none && !tw::IsDriver(whichTool)) {
				logger::DEBUG("parsing tool");
				if (preamble.attaching) {
					std::println(std::cout,"Attaching <{}> to <{}>",preamble.obj_name,preamble.owner_name);
					auto new_tool = CreateTool(preamble.obj_name,whichTool);
					AddTool(new_tool);
					new_tool->ReadInputFileBlock(curs,src);
					auto owner = FindDriver(preamble.owner_name,true);
					if (owner==NULL) {
						tw::input::ThrowParsingError(curs,src,"could not find requested owner");
					} else {
						owner->tools.push_back(new_tool);
					}
				} else {
					logger::INFO(std::format("Creating Tool <{}>...",preamble.obj_name));
					auto new_tool = CreateTool(preamble.obj_name,whichTool);
					AddTool(new_tool);
					if (new_tool->name != preamble.obj_name && errorCheckingLevel>0)
						tw::input::ThrowParsingError(curs,src,"duplicate tool name.");
					new_tool->ReadInputFileBlock(curs,src);
				}
				return tw::input::navigation::gotoSibling;
			}

			// Driver Installation
			if (whichTool!=tw::tool_type::none && tw::IsDriver(whichTool)) {
				logger::DEBUG("parsing driver");
				if (Module::SingularType(whichTool))
					if (most_recent.find(whichTool)!=most_recent.end())
						tw::input::ThrowParsingError(curs,src,"singular driver type was created twice.  Check order of input file.");
				std::println(std::cout,"Installing driver <{}>...",preamble.obj_name);
				Module *super;
				if (preamble.attaching && preamble.owner_name.size() > 0) {
					super = FindDriver(preamble.owner_name,true);
					std::println(std::cout,"Attaching <{}> to <{}>",preamble.obj_name,preamble.owner_name);
				} else if (!preamble.attaching && preamble.owner_name.size() == 0) {
					// If not explicitly attaching, but supermodule is required, find or create one
					tw::tool_type reqType = Module::RequiredSupermoduleType(whichTool);
					if (reqType!=tw::tool_type::none) {
						super = AutoCreateSupers(reqType,preamble.obj_name);
					} else {
						super = Root();
					}
				} else {
					tw::input::ThrowParsingError(curs,src,"supermodule mismatch");
				}
				Module *sub = super->CreateDriver(preamble.obj_name,whichTool);
				super->AddDriver(sub);
				most_recent[whichTool] = sub;
				// this can set off a chain of recursive calls that modify the tree
				sub->ReadInputFileBlock(curs,src);
				return tw::input::navigation::gotoSibling;
			}

			throw tw::FatalError("unknown key <" + preamble.obj_key + "> at " + tw::input::loc_str(curs));

		} else if (typ == "reaction" || typ == "collision" || typ == "excitation") {
			for (auto d : sub_drivers)
				d->ReadQuasitoolBlock(curs,src);
			return tw::input::navigation::gotoSibling;
		} else {
			tw::input::ThrowParsingError(curs,src,"unhandled directive");
			return tw::input::navigation::exit; // compiler doesn't know above line throws
		}


	} else {
		throw tw::FatalError("invalid number of input file passes");
	}
}

Module* Simulation::AutoCreateSupers(tw::tool_type reqType,const std::string& basename)
{
	logger::DEBUG("handling automatic super-driver");
	auto inv = ComputeTool::InvMap();

	// Build the abstract chain
	std::vector<tw::tool_type> chain;
	std::vector<std::string> super_names;
	chain.push_back(reqType);
	super_names.push_back(basename + "_sup");
	while (chain.back()!=tw::tool_type::none) {
		chain.push_back(Module::RequiredSupermoduleType(chain.back()));
		super_names.push_back(super_names.back() + "_sup");
	}
	chain.pop_back();
	super_names.pop_back();

	// Build the real chain to the extent necessary
	Module *curr = this;
	Module *super;
	for (int i = chain.size()-1; i>=0; i--) {
		logger::TRACE(std::format("handling super-driver {}: {}",i,inv[chain[i]]));
		auto it = curr->most_recent.find(chain[i]);
		if (it==curr->most_recent.end()) {
			if (!Module::AutoModuleType(chain[i])) {
				throw tw::FatalError(std::format("<{}> requires super-driver <{}> which cannot be created automatically.",basename,inv[chain[i]]));
			}
			logger::TRACE("(create new one)");
			super = curr->CreateDriver(super_names[i],chain[i]);
			curr->AddDriver(super);
			curr->most_recent[chain[i]] = super;
		} else {
			logger::TRACE("(use existing one)");
			super = (*it).second;
		}
		curr = super;
	}
	return curr;
}

/// The first pass through the input file is used to fully initialize the `Task` and
/// `MetricSpace` objects that will be accessible to every tool/driver.
void Simulation::InputFileFirstPass()
{
	int numRanksProvided,worldRank;
	MPI_Comm_size(MPI_COMM_WORLD,&numRanksProvided);
	MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
	// world rank is suitable for reading task data from restart file
	// because this data is the same in every restart file

	std::stringstream fileName;
	tw::input::FileEnv fenv(task->inputFileName);
	fenv.OpenDeck(raw_src);
	src = tw::input::Preprocess(&fenv,raw_src);
	// std::print(std::cout,"{}",src);
	auto tree = tw::input::GetTree(src);

	// before walking the tree we want to directly setup units
	auto tempNativeUnits = tw::input::GetNativeUnits(tree,src);
	auto tempUnitDensityCGS = tw::input::GetUnitDensityCGS(tree,src);
	space->AttachUnits(tempNativeUnits,tempUnitDensityCGS);
	outerDirectives.AttachUnits(tempNativeUnits,tempUnitDensityCGS);
	outerDirectives.Reset();
	gridReader = std::make_unique<GridReader>(tempNativeUnits,tempUnitDensityCGS);

	// now walk the tree
	inputFilePass = 1;
	tw::input::WalkTree(tree,this);

	outerDirectives.ThrowErrorIfMissingKeys("Simulation");
	if (!gridReader->FoundGrid())
		throw tw::FatalError("Grid directive was not found.");

	auto gdim_idx4 = gridReader->GlobalDims();
	tw::node5 gdim { 1, gdim_idx4.array[1], gdim_idx4.array[2], gdim_idx4.array[3], 1 };

	// Check integer viability
	int64_t totalCellsPerRank = int64_t(gdim[1])*int64_t(gdim[2])*int64_t(gdim[3])/int64_t(numRanksProvided);
	if (totalCellsPerRank>=std::pow(2,31) && sizeof(tw::Int)==4)
		throw tw::FatalError("You must recompile turboWAVE with 64 bit integers to handle this many grid cells.");

	// Verify and (if necessary) correct decomposition
	if (task->NumTasks() != numRanksProvided)
	{
		std::print(std::cout,"{}: Bad decomposition ",term::warning);
		tw::Int ax1=1,ax2=2,ax3=3; // to be sorted so ax1 is longest
		for (tw::Int i=1;i<=3;i++)
		{
			task->domains[i] = 1;
			if (gdim[i]>=gdim[1] && gdim[i]>=gdim[2] && gdim[i]>=gdim[3])
				ax1 = i;
		}
		for (tw::Int i=1;i<=3;i++)
			ax2 = i==ax1 ? ax2 : i;
		for (tw::Int i=1;i<=3;i++)
			ax3 = i==ax1 || i==ax2 ? ax3 : i;
		if (gdim[ax2]<gdim[ax3])
			std::swap(ax2,ax3);
		task->domains[ax1] = numRanksProvided;
		while (gdim[ax1]%(task->domains[ax1]*2)!=0 && task->domains[ax1]>0)
		{
			task->domains[ax1] /= 2;
			task->domains[ax2] *= 2;
		}
		std::println(std::cout,"(defaulting to {}x{}x{})",task->domains[1],task->domains[2],task->domains[3]);
	}
	std::println(std::cout,"{}-Way Decomposition",task->NumTasks());

	// Set up the domain decomposition
	for (auto i=1;i<=3;i++) {
		if (gdim[i]%task->domains[i]!=0)
			throw tw::FatalError("global number of cells is not divisible by number of domains along axis");
	}
	task->Initialize(task->domains,task->periodic);
	for (auto i=1;i<=3;i++) {
		if (gdim[i]>1 && gdim[i]/task->domains[i]%2>0)
			throw tw::FatalError(std::format("local number of cells is not even along non-ignorable axis {}",i));
	}
	space->Resize(task,gdim,gridReader->GlobalCorner(),gridReader->GlobalSize(),std_packing,std_layers,gridReader->Geometry());
	ts_tree_delete(tree);

	// Random numbers
	task->uniformDeviate = new UniformDeviate(1 + task->strip[0].Get_rank()*(MaxSeed()/numRanksProvided));
	task->gaussianDeviate = new GaussianDeviate(1 + task->strip[0].Get_rank()*(MaxSeed()/numRanksProvided) + MaxSeed()/(2*numRanksProvided));
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
}

