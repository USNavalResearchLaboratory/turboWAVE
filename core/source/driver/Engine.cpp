module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"

export module driver:engine;
import :tool;
import :region;

import base;
import input;
import static_space;
import metric_space;
import logger;

export struct Engine : ComputeTool
{
	// Allow for a region
	std::string region_name;
	std::shared_ptr<Region> theRgn;

	Engine(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk) {
		directives.Add("clipping region",new tw::input::String(&region_name),false);
	}
	/// Usually called by the root driver on its engines.
	/// Drivers do not need this since they have access to the rest of the driver tree.
	/// The default attaches a named region, either one named in the input file, or a new EntireRegion.
	virtual void AttachTools(const std::vector<SharedTool>& tools) {
		if (region_name.size()>0) {
			logger::TRACE(std::format("looking for {}",region_name));
			for (auto tool : tools) {
				if (tool->name == region_name) {
					theRgn = std::dynamic_pointer_cast<Region>(tool);
					if (theRgn.use_count()==0) {
						throw tw::FatalError(std::format("could not find region <{}>",region_name));
					}
					return;
				}
			}
		} else {
			logger::TRACE(std::format("create default entire"));
			theRgn = std::make_shared<EntireRegion>("entire",space,task);
			theRgn->Initialize();
		}
	}
};
