module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"

module driver;
import base;
import preprocess;
import input;
import factory;
import logger;

SharedTool Module::CreateTool(const std::string& name,tw::tool_type whichTool) {
	auto new_tool = factory::CreateToolFromType(name,whichTool,space,task);
	tw::input::MangleName(new_tool->name,this->all_names);
    return new_tool;
}

Module* Module::CreateDriver(const std::string& name,tw::tool_type whichDriver) {
	Module *sub = factory::CreateDriverFromType(name,whichDriver,space,task);
	tw::input::MangleName(sub->name,this->all_names);
	return sub;
}

void Module::ParseNestedDeclaration(TSTreeCursor *curs,const std::string& src)
{
	logger::DEBUG("handling nested declaration");
	tw::input::Preamble preamble = tw::input::GetPreamble(curs,src);
	// TODO: always allow `for`, it will apply to anything else on the same level
	// TODO: let's just have one declaration processor which is always nested starting from root
	if (preamble.attaching)
		tw::input::ThrowParsingError(curs,src,"keyword <for> is not allowed in a nested declaration.");

	std::string similar_keys;
	if (!factory::VerifyKey(preamble.obj_key,similar_keys)) {
		tw::input::ThrowParsingError(curs, src, std::format("unknown key <{}>",preamble.obj_key),
			similar_keys.size() > 0 ? std::format("<{}> is similar to {}",preamble.obj_key,similar_keys) : "no similar keys found");
	}

	tw::tool_type whichTool = ComputeTool::CreateTypeFromInput(preamble);
	try {
		if (Module::SingularType(whichTool))
			if (most_recent.find(whichTool)!=most_recent.end())
				tw::input::ThrowParsingError(curs,src,"Singular driver was created twice.  Check order of input file.");
		Module *sub = CreateDriver(preamble.obj_name,whichTool);
		AddDriver(sub);
		std::println(std::cout,"   Attaching nested module <{}>...",sub->name);
		most_recent[whichTool] = sub;
		// The following may lead to a recursive call of this function.
		// If so the tree and maps are being modified each time.
		sub->ReadInputFileBlock(curs,src);
		return;
	} catch (tw::FactoryError& e) {
		auto new_tool = CreateTool(preamble.obj_name,whichTool);
		AddTool(new_tool);
		std::println(std::cout,"   Attaching nested tool <{}>...",new_tool->name);
		new_tool->ReadInputFileBlock(curs,src);
		return;
	}

	tw::input::ThrowParsingError(curs, src, "unprocessed match in nested declaration");
}
