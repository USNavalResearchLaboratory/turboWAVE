module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module navigate;
import base;

extern "C" {
	TSLanguage *tree_sitter_turbowave();
}

export namespace tw
{
	namespace input
	{
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
			virtual navigation descend(TSTreeCursor *curs) { return navigation::exit; };
		};
		/// @brief wraps TSTreeCursor in a way that avoids leaking cursor copies
		class Cursor {
			TSTreeCursor curs;
			public:
			Cursor(TSNode root) {
				curs = ts_tree_cursor_new(root);
			}
			Cursor(const TSTreeCursor *src) {
				curs = ts_tree_cursor_copy(src);
			}
			~Cursor() {
				ts_tree_cursor_delete(&curs);
			}
			TSTreeCursor *get() {
				return &curs;
			}
		};
		TSTree* GetTree(const std::string& src);
		void WalkTree(const TSTree *tree,Visitor *visitor);
		void ThrowParsingError(const TSTreeCursor *curs,const std::string& src,const std::string& messg,const std::string& info="");

		/// @brief goto the next named node, optionally accepting the current node
		/// @param curs tree cursor
		/// @param include_curr whether to accept the current node
		/// @return was a named node found
		inline std::string loc_str(const TSTreeCursor *curs) {
			TSNode n = ts_tree_cursor_current_node(curs);
			TSPoint p = ts_node_start_point(n);
			return std::format("{},{}",p.row,p.column);
		}
		inline std::string node_text(TSNode node,const std::string& src) {
			const int s = ts_node_start_byte(node);
			const int e = ts_node_end_byte(node);
			return src.substr(s,e-s);
		}
		inline std::string node_text(const TSTreeCursor *curs,const std::string &src) {
			const int s = ts_node_start_byte(ts_tree_cursor_current_node(curs));
			const int e = ts_node_end_byte(ts_tree_cursor_current_node(curs));
			return src.substr(s,e-s);
		}
		inline std::string node_kind(const TSTreeCursor *curs) {
			return ts_node_type(ts_tree_cursor_current_node(curs));
		}
		inline bool next_named_node(TSTreeCursor *curs,bool include_curr) {
			if (include_curr) {
				if (ts_node_is_named(ts_tree_cursor_current_node(curs))) {
					return true;
				}
			}
			while (ts_tree_cursor_goto_next_sibling(curs)) {
				if (ts_node_is_named(ts_tree_cursor_current_node(curs))) {
					return true;
				}
			}
			return false;
		}
		/// @brief advance to the next named sibling and return the text, or throw error
		inline std::string next_named_node_text(TSTreeCursor *curs,const std::string &src) {
			if (next_named_node(curs,false)) {
				return node_text(curs,src);
			} else {
				ThrowParsingError(curs,src,"something missing");
				return "";
			}
		}
		/// @brief advance to the next named sibling and return the type, or throw error
		inline std::string next_named_node_kind(TSTreeCursor *curs) {
			if (next_named_node(curs,false)) {
				return node_kind(curs);
			} else {
				throw tw::FatalError("missing syntax node");
			}
		}
		std::string trim(const std::string& str) {
			const std::string whitespace = " \t\n\r\f\v";
			size_t first = str.find_first_not_of(whitespace);
			if (first == std::string::npos) {
				return "";
			}
			size_t last = str.find_last_not_of(whitespace);
			return str.substr(first, (last - first + 1));
		}
		void StripQuotes(std::string& str)
		{
			bool needed = false;
			if (str.length() == 0) {
				return;
			}
			if (str.front()=='\'' && str.back()=='\'')
				needed = true;
			if (str.front()=='\"' && str.back()=='\"')
				needed = true;
			if (needed)
			{
				str.pop_back();
				str = str.substr(1);
			}
		}
	}
}

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

TSTree* tw::input::GetTree(const std::string& src) {
	TSParser * parser = ts_parser_new();
	ts_parser_set_language(parser,tree_sitter_turbowave());
	auto tree = ts_parser_parse_string(parser,NULL,src.c_str(),src.length());
	ts_parser_delete(parser);
	return tree;
}

void tw::input::WalkTree(const TSTree *tree,Visitor *visitor) {
	auto curs = tw::input::Cursor(ts_tree_root_node(tree));
	navigation choice = navigation::gotoSelf;
	while (choice!=navigation::exit && choice!=navigation::abort) {
		if (choice == navigation::gotoSelf) {
			choice = visitor->visit(curs.get());
		} else if (choice == navigation::descend) {
			choice = visitor->descend(curs.get());
		} else if (choice == navigation::gotoChild && ts_tree_cursor_goto_first_child(curs.get())) {
			choice = visitor->visit(curs.get());
		} else if (choice == navigation::gotoParentSibling && ts_tree_cursor_goto_parent(curs.get()) && ts_tree_cursor_goto_next_sibling(curs.get())) {
			choice = visitor->visit(curs.get());
		} else if (choice == navigation::gotoSibling && ts_tree_cursor_goto_next_sibling(curs.get())) {
			choice = visitor->visit(curs.get());
		} else if (ts_tree_cursor_goto_next_sibling(curs.get())) {
			choice = visitor->visit(curs.get());
		} else if (ts_tree_cursor_goto_parent(curs.get())) {
			choice = navigation::gotoSibling;
		} else {
			choice = navigation::exit;
		}
	}
}

void tw::input::ThrowParsingError(const TSTreeCursor *curs,const std::string& src,const std::string& messg,const std::string& info) {
	std::stringstream ss(src);
	std::string line;
	auto linenum = ts_node_start_point(ts_tree_cursor_current_node(curs)).row;
	for (auto i = 0; i<=linenum; i++) {
		std::getline(ss,line);
	}
	auto xmessg = std::format("{} while parsing line\n>>>{}{}{}<<<",messg,term::cyan,line,term::reset_color);
	if (info.size() > 0) {
		xmessg += std::format("\n{}INFO{}: {}",term::cyan,term::reset_color,info);
	}
	throw tw::FatalError(xmessg);
}
