module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"

export module driver:engine;
import :tool;

import base;
import input;
import static_space;
import metric_space;
import logger;

/// This is something short of a driver that can contain other tools,
/// main use is in handling of regions.
export struct Engine : ComputeTool
{
	std::vector<std::string> element_names;
	/// This is a case where we have a tree without a clear ownership pattern, so shared_ptr is used.
	/// We are helped somewhat by the fact that there is no need to point back to the parent.
	std::vector<std::shared_ptr<Engine>> elements;

	std::string region_name;
	std::shared_ptr<Engine> region;

	Engine(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk) {
		directives.Add("elements",new tw::input::StringList(&element_names),false);
		directives.Add("clipping region",new tw::input::String(&region_name),false);
	}
	/// Usually called by the root driver on its engines.
	/// Drivers do not need this since they have access to the rest of the driver tree.
	virtual void AttachTools(const std::vector<SharedTool>& tools) {
		for (auto nm : element_names) {
			logger::TRACE(std::format("looking for {}",nm));
			auto old_count = elements.size();
			for (auto tool : tools) {
				if (tool->name == nm) {
					auto engine = std::dynamic_pointer_cast<Engine>(tool);
					if (engine.use_count()==0) {
						throw tw::FatalError(std::format("<{}> was not an engine",nm));
					} else {
						elements.push_back(engine);
					}
				}
			}
			if (old_count == elements.size()) {
				throw tw::FatalError(std::format("<{}> was missing",nm));
			}
		}
		if (region_name.size() > 0) {
			logger::TRACE(std::format("looking for region {}",region_name));
			for (auto tool : tools) {
				if (tool->name == region_name) {
					auto engine = std::dynamic_pointer_cast<Engine>(tool);
					if (engine.use_count()==0) {
						throw tw::FatalError(std::format("<{}> was not an engine",region_name));
					} else {
						region = engine;
					}
				}
			}
			if (region.use_count()==0) {
				throw tw::FatalError(std::format("<{}> was missing",region_name));
			}
		}
	}
};
