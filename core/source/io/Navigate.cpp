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
		TSTree* GetTree(const std::string& src);
		void WalkTree(const TSTree *tree,Visitor *visitor);

		/// @brief goto the next named node, optionally accepting the current node
		/// @param curs tree cursor
		/// @param include_curr whether to accept the current node
		/// @return was a named node found
		inline std::string loc_str(const TSTreeCursor *curs) {
			TSNode n = ts_tree_cursor_current_node(curs);
			TSPoint p = ts_node_start_point(n);
			return std::format("{},{}",p.row,p.column);
		}
		inline std::string missing(const TSTreeCursor *curs) {
			return std::format("something missing at {}",loc_str(curs));
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
				throw tw::FatalError(missing(curs));
			}
		}
		/// @brief advance to the next named sibling and return the type, or throw error
		inline std::string next_named_node_kind(TSTreeCursor *curs) {
			if (next_named_node(curs,false)) {
				return node_kind(curs);
			} else {
				throw tw::FatalError(missing(curs));
			}
		}
    }
}

TSTree* tw::input::GetTree(const std::string& src) {
	TSParser * parser = ts_parser_new();
	ts_parser_set_language(parser,tree_sitter_turbowave());
	return ts_parser_parse_string(parser,NULL,src.c_str(),src.length());
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

