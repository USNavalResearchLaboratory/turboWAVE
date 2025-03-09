module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module preprocess;
import base;
import navigate;

export namespace tw
{
	namespace input
	{
        typedef std::map<std::string,std::string> Subs;
        class MacroVisitor: public Visitor {
			std::string src, expansion;
			Subs subs;
			bool didChange;
            std::vector<bool> active; // stack for conditionals
            uint32_t last_byte_pos;
			public:
            /// @brief call once before first pass
            /// @param initial_src starting source string
            /// @param map substitutions
            void init(const std::string &initial_src) {
                active.clear();
                subs.clear();
                active.push_back(true);
                this->expansion = initial_src;
            }
			/// @brief call after init, and again before each pass
			void reset() {
                active.clear();
                active.push_back(true);
                last_byte_pos = 0;
                this->src = expansion;
                expansion = "";
                didChange = false;
            }
			/// @brief flush the result string
			/// @return true if there were no more substitutions
			bool done() {
                expansion += src.substr(last_byte_pos);
                return didChange == false;
            }
			std::string* get_src() { return &src; }
			std::string get_exp() { return expansion; }
			virtual navigation visit(TSTreeCursor *curs);
        };
        class IncludeVisitor: public Visitor {
			std::string src, expansion;
            tw::input::FileEnv *env;
            uint32_t last_byte_pos;
			public:
            /// @brief call once before first pass
            /// @param initial_src starting source string
            /// @param map substitutions
            void init(tw::input::FileEnv *env, const std::string& initial_src) {
                last_byte_pos = 0;
                this->env = env;
                this->src = initial_src;
                expansion = "";
            }
            void done() {
                expansion += src.substr(last_byte_pos);
            }
			std::string* get_src() { return &src; }
			std::string get_exp() { return expansion; }
			virtual navigation visit(TSTreeCursor *curs);
        };
        /// @brief parse define node and update substitution map
        /// @param curs should be on define node
        void MacroGather(const TSTreeCursor& curs,Subs& map,const std::string& src);
        /// @brief repeatedly substitute macros until no more references
        /// @param parser tree-sitter parser for TW input files
        /// @param map substitutions
        /// @param src starting source string
        /// @return expanded source string
        std::string MacroSubstitution(const std::string& src);
        /// @brief handle includes and macros recursively
        /// @param src the main source string
        /// @return the expanded source string
        /// @throws syntax error
        std::string Preprocess(tw::input::FileEnv *env,const std::string& src);
    }
}

void tw::input::MacroGather(const TSTreeCursor& curs,tw::input::Subs& map,const std::string& src) {
	if (tw::input::node_kind(&curs)=="define") {
        TSNode node = ts_tree_cursor_current_node(&curs);
        TSNode key_node = ts_node_child(node,1); // node 0 is `#define`
        TSNode val_node = ts_node_child(node,2); // there could be subsequent nodes
		std::string key = tw::input::node_text(key_node,src);
		if (map.find(key)!=map.end()) {
			throw tw::FatalError("Macro "+key+" was already used.");
        }
        if (ts_node_is_null(val_node)) {
            map[key] = std::string("");
        } else {
            // the value in this model is just the text from start of first value node to end of parent
            const int s = ts_node_start_byte(val_node);
            const int e = ts_node_end_byte(node);
            std::string val = trim(src.substr(s,e-s));
            if (key==val)
                throw tw::FatalError("Macro cannot refer to itself ("+key+")");
            map[key] = val;
        }
	} else {
        tw::input::ThrowParsingError(&curs,src,"internal error");
    }
}

tw::input::navigation tw::input::MacroVisitor::visit(TSTreeCursor *curs) {

    if (active.size()==0) {
        tw::input::ThrowParsingError(curs,src,"conditional state is undefined");
    }

	TSNode node = ts_tree_cursor_current_node(curs);

    // if syntax error, hunt down the error and go out
    if (ts_node_has_error(node)) {
        if (ts_node_is_error(node)) {
            tw::input::ThrowParsingError(curs,src,"syntax error");
        } else if (ts_node_is_missing(node)) {
            tw::input::ThrowParsingError(curs,src,"missing " + tw::input::node_kind(curs));
        }
        return tw::input::navigation::gotoChild;
    }

    auto kind = tw::input::node_kind(curs);

    if (active.back()) {
        if (kind=="ifxdef") {
            expansion += src.substr(last_byte_pos,ts_node_start_byte(node) - last_byte_pos);
            ts_tree_cursor_goto_first_child(curs);
            auto which_if = tw::input::node_text(curs,src);
            ts_tree_cursor_goto_next_sibling(curs);
            auto key = tw::input::node_text(curs,src);
            if (which_if=="#ifdef") {
                active.push_back(subs.find(key)!=subs.end());
            } else {
                active.push_back(subs.find(key)==subs.end());
            }
            last_byte_pos = ts_node_end_byte(ts_tree_cursor_current_node(curs));
            return tw::input::navigation::gotoSibling;
        } else if (kind=="else_block") {
            // skip the whole block here and now
            expansion += src.substr(last_byte_pos,ts_node_start_byte(node) - last_byte_pos);
            last_byte_pos = ts_node_end_byte(node);
            return tw::input::navigation::gotoSibling;
        } else if (kind=="#endif") {
            expansion += src.substr(last_byte_pos,ts_node_start_byte(node) - last_byte_pos);
            last_byte_pos = ts_node_end_byte(node);
            active.pop_back();
            return tw::input::navigation::gotoSibling;
        } else if (kind=="define_ref") {
            std::string key = tw::input::node_text(node,src);
            if (subs.find(key)==subs.end()) {
                tw::input::ThrowParsingError(curs,src,"Macro "+key+" not defined.");
            }
            expansion += src.substr(last_byte_pos,ts_node_start_byte(node) - last_byte_pos);
            expansion += subs[key];
            last_byte_pos = ts_node_end_byte(node);
            didChange = true;
            return tw::input::navigation::gotoSibling;
        } else if (kind=="define") {
            tw::input::MacroGather(*curs,subs,src);
            // In this model it makes most sense to strip defines as we encounter them.
            // Otherwise, in the expanded document, we would be facing ambiguous scope.
            expansion += src.substr(last_byte_pos,ts_node_start_byte(node) - last_byte_pos);
            expansion += "\n";
            last_byte_pos = ts_node_end_byte(node);
            return tw::input::navigation::gotoSibling;
        }
    } else {
        if (kind=="else_block") {
            active.pop_back();
            active.push_back(true);
            ts_tree_cursor_goto_first_child(curs); // anonymous #else
            last_byte_pos = ts_node_end_byte(ts_tree_cursor_current_node(curs));
            return tw::input::navigation::gotoSibling;
        } else if (kind=="#endif") {
            last_byte_pos = ts_node_end_byte(node);
            active.pop_back();
            return tw::input::navigation::gotoSibling;
        } else {
            last_byte_pos = ts_node_end_byte(node);
        }
    }

	return tw::input::navigation::gotoChild;
}

/// @brief Handle conditionals and macros
/// @param src the original source string
/// @return the expanded source string
std::string tw::input::MacroSubstitution(const std::string& src) {

    const size_t max_iter = 100;
    size_t iterations = 0;
    tw::input::MacroVisitor visitor;
    visitor.init(src);

	// Do multiple sweeps to account for macros within macros.
	do {
        if (iterations > max_iter) {
            throw tw::FatalError("exceeded maximum depth during macro substitution");
        }
		visitor.reset();
		auto tree = tw::input::GetTree(*visitor.get_src());
		tw::input::WalkTree(tree,&visitor);
		ts_tree_delete(tree);
        iterations++;
	} while (!visitor.done());

	return visitor.get_exp();
}

tw::input::navigation tw::input::IncludeVisitor::visit(TSTreeCursor *curs) {

	TSNode node = ts_tree_cursor_current_node(curs);

    // assume syntax has been verified by macro pass

	if (tw::input::node_kind(curs)=="include") {
        auto fname_node = ts_node_child(node,1);
        if (!ts_node_is_null(fname_node)) {
            std::string include_file;
            std::string fname = tw::input::node_text(fname_node,src);
            tw::input::StripQuotes(fname);
            if (env->FindAndOpen(fname,include_file)) {
                expansion += src.substr(last_byte_pos,ts_node_start_byte(node) - last_byte_pos);
                expansion += tw::input::Preprocess(env,include_file);
                last_byte_pos = ts_node_end_byte(node);
            } else {
                tw::input::ThrowParsingError(curs,src,"could not open file");
            }
            return tw::input::navigation::gotoSibling;
        }
        tw::input::ThrowParsingError(curs,src,"broken tree");
	}
	return tw::input::navigation::gotoChild;
}

std::string tw::input::Preprocess(tw::input::FileEnv *env,const std::string& src) {
    // In order to keep macros scoped to files, it is important that before an include
    // is actually inserted, all its macro substitutions are fully resolved.
    // TODO: let child scopes access the parent scope
    std::string expansion = tw::input::MacroSubstitution(src);
    tw::input::IncludeVisitor visitor;
    visitor.init(env,expansion);
    auto tree = tw::input::GetTree(*visitor.get_src());
    tw::input::WalkTree(tree,&visitor);
    visitor.done();
    ts_tree_delete(tree);
    return visitor.get_exp();
}
