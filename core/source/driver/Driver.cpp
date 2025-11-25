module;

#ifndef USE_STD_MODULE
#include <ranges>
#endif
#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"

export module driver;
export import :tool;
export import :bounded;
export import :region;
export import :warp;
export import :engine;
export import :profile;
export import :diagnostic;

import logger;
import input;
import static_space;
import metric_space;
import tw_iterator;

namespace tw
{
	export enum class priority { diagnostic, source, field };
	inline std::map<tw::priority,tw::Int> priority_sort_map()
	{
		std::map<tw::priority,tw::Int> ans = {{priority::diagnostic,100},{priority::source,200},{priority::field,300}};
		return ans;
	}
}

export struct Driver:StaticSpace,Engine
{
	Driver* super;
	std::vector<Driver*> sub_drivers;
	std::vector<SharedTool> tools;
	/// tool_type is mapped to the most recently created driver of that type
	std::map<tw::tool_type,Driver*> most_recent;
	/// names of all tools or drivers owned by this driver, used for name mangling
	std::set<std::string> all_names;

	tw::priority updateSequencePriority;
	tw::Int subSequencePriority;
	tw::Int smoothing[4],compensation[4];
	bool suppressNextUpdate;

	Driver(const std::string& name,MetricSpace *ms,Task *tsk) : Engine(name,ms,tsk) {
		StaticSpace::operator=(*ms);
		super = NULL;
		updateSequencePriority = tw::priority::source;
		subSequencePriority = 0;
		for (tw::Int i=0;i<4;i++)
			smoothing[i] = compensation[i] = 0;
		suppressNextUpdate = false;
		// Any driver recognizes smoothing keys
		directives.Add("smoothing",new tw::input::Numbers<tw::Int>(&smoothing[1],3),false);
		directives.Add("compensation",new tw::input::Numbers<tw::Int>(&compensation[1],3),false);
	}
	/// make this pointer available to all siblings of this driver
	virtual void PublishResource(void* resource,const std::string& description) {
		for (auto driver : super->sub_drivers)
			driver->InspectResource(resource,description);
	}
	/// check the description to see if this pointer should be copied,
	/// returning whether or not it was copied
	virtual bool InspectResource(void* resource,const std::string& description) {
		return false;
	}
	/// Check inputs after input file has been fully processed, often tool pointers
	/// are checked for existence and created/copied as necessary.
	virtual void VerifyInput() {;}
	virtual void ExchangeResources() {;}
	virtual void Reset() {;}
	virtual void Update() {;}
	virtual void MoveWindow() {;}
	virtual void AntiMoveWindow() {;}
	virtual void AdaptGrid() {;}
	virtual tw::Float AdaptTimestep() { return 0.0; }
	virtual void StartDiagnostics() {;}
	virtual void Report(Diagnostic&) {;}
	virtual void StatusMessage(std::ostream *dest) {;}
	virtual bool ReadQuasitoolBlock(const TSTreeCursor *curs,const std::string& src) {
		return false;
	}

	Driver* Root() {
		tw::Int depth = 0;
		tw::Int max = 16;
		Driver *ans = this;
		while (ans->super != NULL) {
			logger::TRACE(std::format("go up to {}",ans->super->name));
			ans = ans->super;
			depth++;
			if (depth > max) {
				throw tw::FatalError("exceeded max depth (probably broken driver tree)");
			}
		}
		return ans;
	}
	/// @brief search for a named sub-driver on this driver's list, usually called on root
	/// @returns pointer to the driver, or NULL if driver not found
	Driver* FindDriver(const std::string& name,bool recursive) {
		for (auto check : this->sub_drivers) {
			if (name == check->name) {
				return check;
			} else if (recursive) {
				auto within = check->FindDriver(name,recursive);
				if (within!=NULL) {
					return within;
				}
			}
		}
		return NULL;
	}
	/// @brief search for a named tool on this driver's list, usually called on root
	/// @returns shared_ptr to the tool, or empty shared_ptr if tool not found
	SharedTool FindTool(const std::string& name) {
		for (auto check : this->tools) {
			if (name == check->name) {
				return check;
			}
		}
		SharedTool empty_tool;
		return empty_tool;
	}
	/// @brief parse a name and search for it on the driver's list, usually called on root
	/// @param curs should be on a `use` node
	/// @param src source document
	/// @returns the shared tool if found, otherwise throw error
	SharedTool ParseUse(TSTreeCursor *curs,const std::string& src)
	{
		logger::DEBUG("use a named tool");
		ts_tree_cursor_goto_first_child(curs);
		auto word = tw::input::next_named_node_text(curs,src);
		tw::input::StripQuotes(word);
		auto existing_tool = FindTool(word);
		if (existing_tool.use_count()==0) {
			tw::input::ThrowParsingError(curs,src,std::format("No tool <{}> contained in <{}>",name,this->name));
		}
		return existing_tool;
	}
	/// @brief to be called by supermodules that want to add submodules or tools, can be recursive
	/// @param curs should be on the outer directive node
	/// @param src source document
	void ParseNestedDeclaration(TSTreeCursor *curs,const std::string& src);
	/// @brief Create a tool or engine with name mangling
	/// @param name suggested name of tool
	/// @param whichTool type of tool or engine
	/// @returns the shared tool that was created
	SharedTool CreateTool(const std::string& name,tw::tool_type whichTool);
	/// @brief Create a driver with name mangling
	/// @param name suggested name of driver
	/// @param whichTool type of driver
	/// @returns pointer to new driver
	Driver* CreateDriver(const std::string& name,tw::tool_type whichDriver);
	/// @brief add a tool to this driver and to the root driver
	void AddTool(const SharedTool& tool) {
		Driver *root = Root();
		this->tools.push_back(tool);
		// always add shared tools to the root
		if (this!=root) {
			root->tools.push_back(tool);
		}
	}
	/// @brief add sub-driver to this driver
	void AddDriver(Driver *driver) {
		sub_drivers.push_back(driver);
		driver->super = this;
	}
	/// @brief called if an assignment was not handled normally
	/// @param curs on a directive, use, new, generate, or custom assignment
	/// @param src source document
	/// @returns whether any directive was handled
	virtual bool ReadInputFileDirective(const TSTreeCursor *curs0,const std::string& src)
	{
		std::string command = tw::input::node_kind(curs0);
		// Get an existing tool by searching for a name -- ``use <name>``
		if (command=="use") {
			auto curs = tw::input::Cursor(curs0);
			auto new_tool = this->Root()->ParseUse(curs.get(),src);
			this->tools.push_back(new_tool);
			return true;
		}

		// Take care of nested declarations
		if (command=="new" || command=="generate") {
			auto curs = tw::input::Cursor(curs0);
			this->ParseNestedDeclaration(curs.get(),src);
			return true;
		}

		return false;
	}
	void RecursiveTreeDisplay(std::ostream *theStream,int currLevel,int maxLevel) {
		if (currLevel > maxLevel) {
			throw tw::FatalError("reached max recursions (probably broken driver tree)");
		}
		std::string indentation = "";
		for (auto i=0; i< currLevel*2; i++) {
			indentation += " ";
		}
		*theStream << indentation << "Tools:" << std::endl;
		for (auto i=0; i<tools.size(); i++) {
			*theStream << indentation << "  " << i+1 << ". " << tools[i]->name << std::endl;			
		}
		// for (auto const& [i,t] : std::views::enumerate(tools)) {
		// 	*theStream << indentation << "  " << i+1 << ". " << t->name << std::endl;
		// }
		*theStream << indentation << "Drivers:" << std::endl;
		for (auto i=0; i<sub_drivers.size(); i++) {
			*theStream << indentation << "  " << i+1 << ". " << sub_drivers[i]->name << std::endl;
			sub_drivers[i]->RecursiveTreeDisplay(theStream,currLevel+1,maxLevel);
		}
		// for (auto const& [i,d] : std::views::enumerate(sub_drivers)) {
		// 	*theStream << indentation << "  " << i+1 << ". " << d->name << std::endl;
		// 	d->RecursiveTreeDisplay(theStream,currLevel+1,maxLevel);
		// }
	}
	static bool SingularType(tw::tool_type theType) {
		return theType==tw::tool_type::kinetics || theType==tw::tool_type::sparcHydroManager;
	}
	static bool AutoModuleType(tw::tool_type theType) {
		return theType==tw::tool_type::kinetics || theType==tw::tool_type::equilibriumGroup;
	}
	static tw::tool_type RequiredSupermoduleType(tw::tool_type submoduleType) {
		std::map<tw::tool_type,tw::tool_type> containmentMap =
		{
			{tw::tool_type::species,tw::tool_type::kinetics},
			{tw::tool_type::chemical,tw::tool_type::equilibriumGroup},
			{tw::tool_type::equilibriumGroup,tw::tool_type::sparcHydroManager}
		};
		if (containmentMap.find(submoduleType)==containmentMap.end())
			return tw::tool_type::none;
		else
			return containmentMap[submoduleType];
	}
};

export struct DriverComparator
{
	bool operator() (Driver* const& m1,Driver* const& m2)
	{
		tw::Int idx1 = tw::priority_sort_map()[m1->updateSequencePriority] + m1->subSequencePriority;
		tw::Int idx2 = tw::priority_sort_map()[m2->updateSequencePriority] + m2->subSequencePriority;
		return idx1 < idx2;
	}
};
