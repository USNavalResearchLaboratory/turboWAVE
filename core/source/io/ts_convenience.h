namespace tw
{
	namespace input
	{
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