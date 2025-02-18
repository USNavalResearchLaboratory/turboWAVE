#include <meta_base.h>

// We are using explicit instantiation approach to keep some templates here.

template <class T>
bool tw::input::Numbers<T>::Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
{
    next_named_node(curs,false);
    if (num>1) {
        if (tw::input::node_kind(curs) == "tuple") {
            ts_tree_cursor_goto_first_child(curs);
        } else {
            throw tw::FatalError("expected tuple at " + loc_str(curs));
        }
    }
    for (tw::Int i=0;i<num;i++) {
        if (next_named_node(curs,true)) {
            dat[i] = tw::dnum(curs,src) >> native;
        } else {
            throw tw::FatalError(missing(curs));
        }
        ts_tree_cursor_goto_next_sibling(curs);
    }
    return true;
}

template struct tw::input::Numbers<tw::Float>;
template struct tw::input::Numbers<tw::Int>;

bool tw::input::String::Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
{
    if (next_named_node(curs,false)) {
        if (tw::input::node_kind(curs) == "string_literal") {
            *str = node_text(curs,src);
            StripQuotes(*str);
        } else {
            throw tw::FatalError("expected string at " + loc_str(curs));
        }
    } else {
        throw tw::FatalError(missing(curs));
    }
    return true;
}

bool tw::input::Bool::Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
{
    std::string word = next_named_node_text(curs,src);
    std::multimap<std::string,bool> m = {{"on",true},{"true",true},{"yes",true},{"off",false},{"false",false},{"no",false}};
    auto p = m.find(word);
    if (p==m.end())
        throw tw::FatalError("Invalid boolean value <"+word+"> while reading key <"+key+">");
    *dat = p->second;
    return true;
}
