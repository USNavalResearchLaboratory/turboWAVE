module;

#include "meta_base.h"

export module input;

export import assignment;

export namespace tw
{
	namespace input
	{
		enum class navigation {
			gotoSelf,
			gotoChild,
			gotoSibling,
			gotoParentSibling,
			descend,
			exit,
			abort
		};
		class Visitor {
			public:
			virtual navigation visit(TSTreeCursor *curs) = 0;
			virtual navigation descend(TSTreeCursor *curs) = 0;
		};
		class MacroVisitor : public Visitor {
			std::string src;
			std::map<std::string,std::string> dict;
			std::string expandedDeck;
			bool didChange;
			public:
			void init() { didChange = false; }
			bool done() { return didChange == false; }
			std::string get() { return expandedDeck; }
			virtual navigation visit(TSTreeCursor *curs);
			virtual navigation descend(TSTreeCursor *curs);
		};
		inline std::string trim(const std::string& str) {
			const std::string whitespace = " \t\n\r\f\v";
			size_t first = str.find_first_not_of(whitespace);
			if (first == std::string::npos) {
				return "";
			}
			size_t last = str.find_last_not_of(whitespace);
			return str.substr(first, (last - first + 1));
		}
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
		/// @brief encapsulates finding and opening files
		struct FileEnv
		{
			std::string inputFileName;
			std::vector<std::string> searchPaths;
			FileEnv(const std::string& inputFileName);
			bool OpenDeck(std::string& contents) const;
			bool FindAndOpen(const std::string& fileName,std::string& contents) const;
			bool FindAndOpen(const std::string& fileName,std::stringstream& contents) const;
		};

		TSTree* GetTree(const std::string& src);
		void WalkTree(const TSTree *tree,Visitor *visitor);

		tw::Float GetUnitDensityCGS(TSTree *tree,const std::string& src);
		tw::units GetNativeUnits(TSTree *tree,const std::string& src);
		std::string MacroSubstitution(const std::string& src);
		Preamble GetPreamble(TSTreeCursor *curs,const std::string& src);
		std::string PythonRange(TSTreeCursor *curs,const std::string& src,tw::Float *v0,tw::Float *v1);
		void ThrowParsingError(const TSTreeCursor *curs,const std::string& src,const std::string& messg);
	}
}

extern "C" {
	TSLanguage *tree_sitter_turbowave();
}

void tw::input::ThrowParsingError(const TSTreeCursor *curs,const std::string& src,const std::string& messg) {
	std::string xmessg("While parsing ");
	xmessg += tw::input::node_kind(curs) + " = " + tw::input::node_text(curs,src);
	xmessg += "\n";
	xmessg += "on line " + std::to_string(ts_node_start_point(ts_tree_cursor_current_node(curs)).row);
	xmessg += ":\n" + messg;
	throw tw::FatalError(xmessg);
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
	TSTreeCursor curs = ts_tree_cursor_copy(curs0);
	if (tw::input::node_kind(&curs) == "assignment") {
		ts_tree_cursor_goto_first_child(&curs);
		std::string key = input::node_text(&curs,src);
		if (dmap.find(key)!=dmap.end()) {
			keysFound[key] = 0;
			if (native.ne==0.0)
				throw tw::FatalError("DirectiveReader class is missing units while processing key <"+key+">.");
			return dmap[key]->Read(&curs,src,key,native);
		} else {
			throw tw::FatalError("Unknown key <"+key+">.");
		}
	} else {
		return false;
	}
}

/// @brief read all assignments in this block
/// @param curs current cursor position is on first child in block
/// @param src text of the source document
/// @throw if non-assignment encountered
void tw::input::DirectiveReader::ReadAll(TSTreeCursor *curs,const std::string& src)
{
	do {
		if (!ReadNext(curs,src)) {
			throw tw::FatalError("non-assignment encountered");
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


////////////////////////////////////
//                                //
//  Pre-Processing of Input File  //
//                                //
////////////////////////////////////

tw::input::FileEnv::FileEnv(const std::string& inputFilename)
{
	// Purpose of this class is mainly to encapsulate file searches
	this->inputFileName = inputFilename;
	searchPaths.push_back("");
	std::size_t sep_pos = inputFileName.find_last_of("/\\");
	searchPaths.push_back(inputFileName.substr(0,sep_pos+1)); // keep separator so we have the right one later
}

bool tw::input::FileEnv::OpenDeck(std::string& contents) const
{
	std::ifstream inFile;
	std::stringstream temp;
	inFile.open(inputFileName);
	if (inFile.good())
	{
		temp << inFile.rdbuf();
		contents = temp.str();
		inFile.close();
		return true;
	}
	return false;
}

bool tw::input::FileEnv::FindAndOpen(const std::string& fileName,std::string& contents) const
{
	std::ifstream *inFile;
	std::stringstream temp;
	for (auto it=searchPaths.begin();it!=searchPaths.end();++it)
	{
		inFile = new std::ifstream(*it + fileName);
		if (inFile->good())
		{
			temp << inFile->rdbuf();
			contents = temp.str();
			inFile->close();
			delete inFile;
			return true;
		}
		delete inFile;
	}
	return false;
}

bool tw::input::FileEnv::FindAndOpen(const std::string& fileName,std::stringstream& contents) const
{
	std::ifstream *inFile;
	for (auto it=searchPaths.begin();it!=searchPaths.end();++it)
	{
		inFile = new std::ifstream(*it + fileName);
		if (inFile->good())
		{
			contents << inFile->rdbuf();
			inFile->close();
			delete inFile;
			return true;
		}
		delete inFile;
	}
	return false;
}

tw::Float tw::input::GetUnitDensityCGS(TSTree *tree,const std::string& src) {
	TSTreeCursor curs = ts_tree_cursor_new(ts_tree_root_node(tree));
	if (ts_tree_cursor_goto_first_child(&curs)) {
		do {
			if (tw::input::node_kind(&curs) == "assignment") {
				ts_tree_cursor_goto_first_child(&curs);
				if (input::node_text(&curs,src) == "unit density") {
					ts_tree_cursor_goto_next_sibling(&curs);
					ts_tree_cursor_goto_next_sibling(&curs);
					if (tw::input::node_kind(&curs) == "decimal") {
						return std::stod(input::node_text(&curs,src));
					} else {
						throw tw::FatalError("expected raw value for unit density");
					}
				}
				ts_tree_cursor_goto_parent(&curs);
			}
		} while (ts_tree_cursor_goto_next_sibling(&curs));
		throw tw::FatalError("could not find unit density assignment");
	}
	throw tw::FatalError("empty input file");
}

tw::units tw::input::GetNativeUnits(TSTree *tree,const std::string& src) {
	std::map<std::string,tw::units> m = tw::get_unit_map();
	TSTreeCursor curs = ts_tree_cursor_new(ts_tree_root_node(tree));
	if (ts_tree_cursor_goto_first_child(&curs)) {
		do {
			if (tw::input::node_kind(&curs) == "assignment") {
				ts_tree_cursor_goto_first_child(&curs);
				if (input::node_text(&curs,src) == "native units") {
					ts_tree_cursor_goto_next_sibling(&curs);
					ts_tree_cursor_goto_next_sibling(&curs);
					if (tw::input::node_kind(&curs) == "identifier") {
						std::string txt = input::node_text(&curs,src);
						if (m.find(txt) == m.end())
							throw tw::FatalError("Unknown system of units <" + txt + ">.");
						else
							return m[txt];
					} else {
						throw tw::FatalError("expected identifier for native units");
					}
				}
			}
		} while (ts_tree_cursor_goto_next_sibling(&curs));
		return tw::units::plasma;
	}
	throw tw::FatalError("empty input file");
}

void tw::input::WalkTree(const TSTree *tree,Visitor *visitor) {
	TSTreeCursor curs = ts_tree_cursor_new(ts_tree_root_node(tree));
	navigation choice = navigation::gotoSelf;
	while (choice!=navigation::exit && choice!=navigation::abort) {
		if (choice == navigation::gotoSelf) {
			choice = visitor->visit(&curs);
		} else if (choice == navigation::descend) {
			choice = visitor->descend(&curs);
		} else if (choice == navigation::gotoChild && ts_tree_cursor_goto_first_child(&curs)) {
			choice = visitor->visit(&curs);
		} else if (choice == navigation::gotoParentSibling && ts_tree_cursor_goto_parent(&curs) && ts_tree_cursor_goto_next_sibling(&curs)) {
			choice = visitor->visit(&curs);
		} else if (choice == navigation::gotoSibling && ts_tree_cursor_goto_next_sibling(&curs)) {
			choice = visitor->visit(&curs);
		} else if (ts_tree_cursor_goto_next_sibling(&curs)) {
			choice = visitor->visit(&curs);
		} else if (ts_tree_cursor_goto_parent(&curs)) {
			choice = navigation::gotoSibling;
		} else {
			choice = navigation::exit;
		}
	}
}

tw::input::navigation tw::input::MacroVisitor::visit(TSTreeCursor *curs) {
	TSNode node = ts_tree_cursor_current_node(curs);
	if (tw::input::node_kind(curs)=="define") {
		std::string key = tw::input::node_text(ts_node_child(node,1),src);
		// value node is hidden and may resolve to a sequence of any length,
		// so we simply go to end of the parent.
		const int s = ts_node_start_byte(ts_node_child(node,2));
		const int e = ts_node_end_byte(node);
		std::string val = src.substr(s,e);
		if (key==trim(val))
			throw tw::FatalError("Macro cannot refer to itself ("+key+")");
		if (dict.find(key)==dict.end())
			dict[key] = val;
		else
			throw tw::FatalError("Macro "+key+" was already used.");
		expandedDeck += tw::input::node_text(node,src);
		expandedDeck += " ";
		return tw::input::navigation::gotoSibling;
	}
	if (tw::input::node_kind(curs)=="define_ref") {
		std::string key = tw::input::node_text(node,src);
		if (dict.find(key)==dict.end()) {
			throw tw::FatalError("Macro "+key+" not defined.");
		}
		expandedDeck += trim(dict[key]);
		expandedDeck += " ";
		return tw::input::navigation::gotoSibling;
	}
	if (ts_node_child_count(node)==0) {
		expandedDeck += tw::input::node_text(node,src);
	}
	return tw::input::navigation::gotoChild;
}

tw::input::navigation tw::input::MacroVisitor::descend(TSTreeCursor *curs) {
	return tw::input::navigation::exit;
}

std::string tw::input::MacroSubstitution(const std::string& src) {
	tw::input::MacroVisitor visitor;
	TSParser *parser = ts_parser_new();
	ts_parser_set_language(parser,tree_sitter_turbowave());
	std::string new_src(src);

	// Do multiple sweeps to account for macros within macros.
	do {
		visitor.init();
		TSTree *tree = ts_parser_parse_string(parser,NULL,new_src.c_str(),new_src.length());
		tw::input::WalkTree(tree,&visitor);
		ts_tree_delete(tree);
		new_src = visitor.get();
	} while (!visitor.done());

	return visitor.get();
}

TSTree* tw::input::GetTree(const std::string& src) {
	TSParser * parser = ts_parser_new();
	ts_parser_set_language(parser,tree_sitter_turbowave());
	return ts_parser_parse_string(parser,NULL,src.c_str(),src.length());
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
		throw tw::FatalError("expected range node");
	
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
	tw::input::Preamble ans;
	std::string root = tw::input::node_kind(curs);
	ts_tree_cursor_goto_first_child(curs);
	ans.obj_key = next_named_node_text(curs,src);
	ans.obj_name = ans.obj_key;
	ans.owner_name = "";
	ts_tree_cursor_goto_next_sibling(curs);
	if (tw::input::node_kind(curs) == "string_literal") {
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
