module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"

export module input;
export import base;
export import navigate;
export import assignment;
export import units;
import logger;

extern "C" {
	TSLanguage *tree_sitter_turbowave();
}

export namespace tw
{
	namespace input
	{
		/// @brief Assignments can be added to this object for automatic parsing
		class DirectiveReader
		{
			tw::Int maxKeyWords;
			std::vector<std::string> requiredKeys;
			std::map<std::string,tw::Int> keysFound;
			std::map<std::string,tw::input::Assignment*> dmap;
			tw::UnitConverter native;
		public:
			DirectiveReader();
			~DirectiveReader();
			void AttachUnits(const tw::UnitConverter& native);
			void Reset();
			void Add(const std::string& key,tw::input::Assignment *dir,bool required=true);
			bool ReadNext(const TSTreeCursor *curs,const std::string& src);
			void ReadAll(TSTreeCursor *curs,const std::string& src);
			bool TestKey(const std::string& test);
			void ThrowErrorIfMissingKeys(const std::string& src);
		};
		/// @brief encapsulates information about how an object is being created in the input file
		struct Preamble
		{
			bool attaching;
			std::string obj_key,obj_name,owner_name;
		};
		tw::Float GetUnitDensityCGS(TSTree *tree,const std::string& src);
		tw::units GetNativeUnits(TSTree *tree,const std::string& src);
		Preamble GetPreamble(TSTreeCursor *curs,const std::string& src);
		void BuildSimilar(std::string& messg,const std::string& wrong,const std::string& valid);
		std::string PythonRange(TSTreeCursor *curs,const std::string& src,tw::Float *v0,tw::Float *v1);
	}
}


////////////////////////
//                    //
//  Directive Reader  //
//                    //
////////////////////////

tw::input::DirectiveReader::DirectiveReader()
{
	maxKeyWords = 0;
}

tw::input::DirectiveReader::~DirectiveReader()
{
	for (auto pair : dmap)
		delete pair.second;
}

void tw::input::DirectiveReader::AttachUnits(const tw::UnitConverter& native)
{
	this->native = native;
}

void tw::input::DirectiveReader::Reset()
{
	keysFound.clear();
}

/// @brief This important function adds an assignment to the reader.  Modules and tools will use this
///        to setup automatic parsing of the input file.
/// @param key input file key triggering this assignment
/// @param assignment specific subclass of Assignment object
/// @param required is this assignment required
void tw::input::DirectiveReader::Add(const std::string& key,tw::input::Assignment *assignment,bool required)
{
	dmap[key] = assignment;
	// count number of words in this key
	std::string word;
	std::stringstream s(key);
	tw::Int words = 1;
	while (!s.eof())
	{
		s >> word;
		if (!s.eof())
			words++;
	}
	// update max words
	if (words>maxKeyWords)
		maxKeyWords = words;
	// update required keys
	if (required)
		requiredKeys.push_back(key);
}

/// @brief process an assignment
/// @param curs cursor on a directive
/// @param src text of the source document
/// @return false if directive not an assignment, or custom assignment
bool tw::input::DirectiveReader::ReadNext(const TSTreeCursor *curs0,const std::string& src)
{
	auto curs = tw::input::Cursor(curs0);
	if (tw::input::node_kind(curs.get()) == "assignment") {
		ts_tree_cursor_goto_first_child(curs.get());
		std::string key = input::node_text(curs.get(),src);
		logger::TRACE(std::format("handling key {}",key));
		if (dmap.find(key)!=dmap.end()) {
			keysFound[key] = 0;
			if (native.ne==0.0) {
				tw::input::ThrowParsingError(curs0,src,"units missing");
			}
			auto is_normal = dmap[key]->Read(curs.get(),src,key,native);
			return is_normal;
		} else {
			tw::input::ThrowParsingError(curs0,src,"unknown key");
		}
	}
	return false;
}

/// @brief read all assignments in this block
/// @param curs current cursor position is on first child in block
/// @param src text of the source document
/// @throw if non-assignment encountered
void tw::input::DirectiveReader::ReadAll(TSTreeCursor *curs,const std::string& src)
{
	do {
		if (!ReadNext(curs,src)) {
			auto kind = tw::input::node_kind(curs);
			if (kind != "comment" && kind != "}") {
				tw::input::ThrowParsingError(curs,src,"expected assignment, got " + kind);
			}
		}
	} while (ts_tree_cursor_goto_next_sibling(curs));
}

bool tw::input::DirectiveReader::TestKey(const std::string& test)
{
	return (keysFound.find(test)!=keysFound.end());
}

void tw::input::DirectiveReader::ThrowErrorIfMissingKeys(const std::string& src)
{
	for (auto s : requiredKeys)
		if (!TestKey(s))
			throw tw::FatalError("Missing required key <"+s+"> in <"+src+">");
}

tw::Float tw::input::GetUnitDensityCGS(TSTree *tree,const std::string& src) {
	auto curs = tw::input::Cursor(ts_tree_root_node(tree));
	if (ts_tree_cursor_goto_first_child(curs.get())) {
		do {
			if (tw::input::node_kind(curs.get()) == "assignment") {
				ts_tree_cursor_goto_first_child(curs.get());
				if (input::node_text(curs.get(),src) == "unit density") {
					ts_tree_cursor_goto_next_sibling(curs.get());
					ts_tree_cursor_goto_next_sibling(curs.get());
					if (tw::input::node_kind(curs.get()) == "decimal") {
						auto ans = std::stod(input::node_text(curs.get(),src));
						return ans;
					} else {
						throw tw::FatalError("expected raw value for unit density");
					}
				}
				ts_tree_cursor_goto_parent(curs.get());
			}
		} while (ts_tree_cursor_goto_next_sibling(curs.get()));
		throw tw::FatalError("could not find unit density assignment");
	}
	throw tw::FatalError("empty input file");
}

tw::units tw::input::GetNativeUnits(TSTree *tree,const std::string& src) {
	std::map<std::string,tw::units> m = tw::get_unit_map();
	auto curs = tw::input::Cursor(ts_tree_root_node(tree));
	if (ts_tree_cursor_goto_first_child(curs.get())) {
		do {
			if (tw::input::node_kind(curs.get()) == "assignment") {
				ts_tree_cursor_goto_first_child(curs.get());
				if (input::node_text(curs.get(),src) == "native units") {
					ts_tree_cursor_goto_next_sibling(curs.get());
					ts_tree_cursor_goto_next_sibling(curs.get());
					if (tw::input::node_kind(curs.get()) == "identifier") {
						std::string txt = input::node_text(curs.get(),src);
						if (m.find(txt) == m.end()) {
							throw tw::FatalError("Unknown system of units <" + txt + ">.");
						} else {
							return m[txt];
						}
					} else {
						throw tw::FatalError("expected identifier for native units");
					}
				}
				ts_tree_cursor_goto_parent(curs.get());
			}
		} while (ts_tree_cursor_goto_next_sibling(curs.get()));
		return tw::units::plasma;
	}
	throw tw::FatalError("empty input file");
}

/// @brief read python style range
/// @param curs on range node
/// @param src source document
/// @param v0 start of range, blank resolves to 0
/// @param v1 end of range, blank resolves to big_pos
/// @returns name of quantity, e.g., `Te(5:10)` returns "Te"
std::string tw::input::PythonRange(TSTreeCursor *curs,const std::string& src,tw::Float *v0,tw::Float *v1)
{
	if (tw::input::node_kind(curs) != "range")
		tw::input::ThrowParsingError(curs,src,"expected range node");
	
	ts_tree_cursor_goto_first_child(curs);
	std::string ans = tw::input::node_text(curs,src);
	ts_tree_cursor_goto_next_sibling(curs); // `(`
	ts_tree_cursor_goto_next_sibling(curs);
	if (tw::input::node_text(curs,src)==":")
		*v0 = 0.0;
	else
		*v0 = std::stod(tw::input::node_text(curs,src));
	if (tw::input::next_named_node(curs,false))
		*v1 = std::stod(tw::input::node_text(curs,src));
	else
		*v1 = tw::big_pos;
	return ans;
}

/// @brief get opening parts of any new or generate directive
/// @param curs cursor, should be on parent of the preamble items
/// @param src source document
/// @return the preamble data, cursor is left on the block node
tw::input::Preamble tw::input::GetPreamble(TSTreeCursor *curs,const std::string& src) {
	bool found_string = false;
	tw::input::Preamble ans;
	std::string root = tw::input::node_kind(curs);
	ts_tree_cursor_goto_first_child(curs);
	ans.obj_key = next_named_node_text(curs,src);
	ans.obj_name = ans.obj_key;
	ans.owner_name = "";
	ts_tree_cursor_goto_next_sibling(curs);
	if (tw::input::node_kind(curs) == "string_literal") {
		found_string = true;
		ans.obj_name = node_text(curs,src);
		ts_tree_cursor_goto_next_sibling(curs);
	}
	if (root == "new") {
		ans.attaching = false;
	} else if (root == "associative_new") {
		ans.attaching = true;
		ts_tree_cursor_goto_next_sibling(curs);
		ans.owner_name = node_text(curs,src);
	} else if (root == "generate") {
		if (!found_string) {
			ts_tree_cursor_goto_previous_sibling(curs);
			ThrowParsingError(curs, src, "bad generate, is the name quoted?");
		}
		ans.attaching = true;
		ans.owner_name = ans.obj_name;
		ans.obj_name = ans.obj_key;
	}
	while (tw::input::node_kind(curs) != "block") {
		if (!ts_tree_cursor_goto_next_sibling(curs)) {
			throw tw::FatalError("missing block");
		}
	}
	tw::input::StripQuotes(ans.obj_key);
	tw::input::StripQuotes(ans.obj_name);
	tw::input::StripQuotes(ans.owner_name);
	return ans;
}

/// @brief build a message to display with valid keys that are similar to an erroneous one
/// @param messg a message string that will be appended after each call
/// @param wrong a string that failed to match
/// @param valid one of the valid strings we are looking for
void tw::input::BuildSimilar(std::string& messg,const std::string& wrong,const std::string& valid) {
	// levenshtein distance (this is an AI program)
    const size_t len1 = wrong.size(), len2 = valid.size();
    std::vector<std::vector<int>> d(len1 + 1, std::vector<int>(len2 + 1));

    for (size_t i = 0; i <= len1; ++i) d[i][0] = i;
    for (size_t i = 0; i <= len2; ++i) d[0][i] = i;

    for (size_t i = 1; i <= len1; ++i) {
        for (size_t j = 1; j <= len2; ++j) {
            int cost = (wrong[i - 1] == valid[j - 1]) ? 0 : 1;
            d[i][j] = std::min({d[i - 1][j] + 1,          // deletion
                               d[i][j - 1] + 1,          // insertion
                               d[i - 1][j - 1] + cost}); // substitution
        }
    }
	// d[len1][len2] contains the number of edits needed to transform
	// one string into the other
    if (d[len1][len2] < std::abs(int(len1-len2))+1) {
		messg += valid;
		messg += ", ";
		return;
	}
	if (d[len1][len2] < 5*std::max(len1,len2)/10) {
		messg += valid;
		messg += ", ";
		return;
	}
}