module;

#include "tree_sitter/api.h"
#include "tw_includes.h"

export module assignment;
import base;
import tensor;
import navigate;
import units;

export namespace tw
{
	namespace input
	{
		/// @brief Assignment is a key element in automatic parsing if modules and tools.
        /// By adding assignments to a `DirectiveReader` various data forms can be parsed and error checked for free.
        /// Internally this works via a hash table fo input file keys mapped to `Assignment` instances.
        struct Assignment
        {
            virtual ~Assignment()
            {
            }
            /// @brief parse and capture data if this is an assignment
            /// @param curs on first child of assignment node
            /// @param src source document
            /// @param key the text of the current node
            /// @param native converter for native units
            /// @return was the assignment processed, e.g., will be false for custom directive
            virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
            {
                return true;
            }
        };
        struct Custom : Assignment
        {
            virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
            {
                return false;
            }
        };
        /// @brief capture literal-string or identifier, optionally casts any node to a string
        struct String : Assignment
        {
            std::string *str;
            bool allow_cast;
            String(std::string *s,bool allow_cast = false) { str=s; this->allow_cast = allow_cast; }
            virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native);
        };
        /// @brief this is a tuple supporting raw or dimensional numbers
        /// @tparam T numerical type
        template <class T>
        struct Numbers : Assignment
        {
            T *dat;
            tw::Int num;
            Numbers(T *d,tw::Int n) { num=n; dat=d; }
            virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native);
        };
        struct Int : Numbers<tw::Int>
        {
            Int(tw::Int *i) : Numbers<tw::Int>(i,1) { ; }
        };
        struct Float : Numbers<tw::Float>
        {
            Float(tw::Float *f) : Numbers<tw::Float>(f,1) { ; }
        };
        struct Complex : Numbers<tw::Float>
        {
            Complex(tw::Complex *c) : Numbers<tw::Float>((tw::Float*)c,2) { ; }
        };
        struct Vec3 : Numbers<tw::Float>
        {
            Vec3(tw::vec3 *v) : Numbers<tw::Float>(&v->x,3) { ; }
        };
        struct Vec4 : Numbers<tw::Float>
        {
            Vec4(tw::vec4 *v) : Numbers<tw::Float>(&((*v)[0]),4) { ; }
        };
        /// @brief variable length array of numbers
        /// @tparam T the container type including its inner type
        template <class T>
        struct NumberList : Assignment
        {
            T *dat;
            NumberList(T *d) { dat=d; }
            virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native);
        };
        /// @brief variable length array of strings
        /// @tparam T the container type including its inner type (should be std::string)
        template <class T>
        struct StringList : Assignment {
            T *dat;
            bool allow_cast;
            StringList(T *d, bool allow_cast = false) { dat=d; this->allow_cast = allow_cast; }
            virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native);
        };
        template <class T>
        struct Enums : Assignment
        {
            T *dat1,*dat2,*dat3;
            std::map<std::string,T> *emap; // has to be pointer for polymorphic destruction
            Enums(const std::map<std::string,T>& m,T* d1,T* d2=NULL,T* d3=NULL) {
                // must specify at least 1 data element, may specify up to 3.
                dat1=d1; dat2=d2; dat3=d3;
                emap = new std::map<std::string,T>(m);
            }
            virtual ~Enums() {
                delete emap;
            }
            virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native);
        };
        struct Bool : Assignment
        {
            bool *dat;
            Bool(bool *d) { dat=d; }
            virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native);
        };
    }
}

template <class T>
bool tw::input::NumberList<T>::Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
{
    using F = typename T::value_type;
    std::vector<F> temp;
    if (next_named_node(curs,false)) {
        if (node_kind(curs)=="list") {
            if (ts_tree_cursor_goto_first_child(curs)) {
                next_named_node(curs,true);
                do {
                    auto x = tw::dnum(curs,src) >> native;
                    temp.push_back(x);    
                } while (next_named_node(curs,false));
                ts_tree_cursor_goto_parent(curs);
            } else {
                tw::input::ThrowParsingError(curs,src,"empty list");
            }
        } else {
            tw::input::ThrowParsingError(curs,src,"expected list");
        }
    } else {
        tw::input::ThrowParsingError(curs,src,"something missing");
    }
    dat->resize(temp.size());
    for (tw::Int i=0;i<temp.size();i++)
        (*dat)[i] = temp[i];
    return true;
}

template <class T>
bool tw::input::StringList<T>::Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
{
    using F = typename T::value_type;
    std::vector<F> temp;
    if (next_named_node(curs,false)) {
        if (node_kind(curs)=="list") {
            if (ts_tree_cursor_goto_first_child(curs)) {
                next_named_node(curs,true);
                do {
                    if (node_kind(curs)=="string_literal") {
                        ts_tree_cursor_goto_first_child(curs);
                        ts_tree_cursor_goto_next_sibling(curs);
                        temp.push_back(node_text(curs,src));
                        ts_tree_cursor_goto_parent(curs);
                    } else if (node_kind(curs)=="identifier" || allow_cast) {
                        temp.push_back(node_text(curs,src));
                    } else {
                        tw::input::ThrowParsingError(curs,src,"expected string or identifier");
                    }
                } while (next_named_node(curs,false));
                ts_tree_cursor_goto_parent(curs);
            } else {
                tw::input::ThrowParsingError(curs,src,"empty list");
            }
        } else {
            tw::input::ThrowParsingError(curs,src,"expected list");
        }
    } else {
        tw::input::ThrowParsingError(curs,src,"something missing");
    }
    dat->resize(temp.size());
    for (tw::Int i=0;i<temp.size();i++)
        (*dat)[i] = temp[i];
    return true;
}

template <class T>
bool tw::input::Enums<T>::Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
{
    if (next_named_node_kind(curs) == "tuple") {
        ts_tree_cursor_goto_first_child(curs);
    }
    next_named_node(curs,true);
    auto read_one = [&] (T *d)
    {
        if (d!=NULL)
        {
            std::string word = node_text(curs,src);
            if (emap->find(word)==emap->end()) {
                tw::input::ThrowParsingError(curs,src,"invalid type");
            }
            *d = (*emap)[word];
            next_named_node(curs,false);
        }
    };
    read_one(dat1);
    read_one(dat2);
    read_one(dat3);
    return true;
}

// We are using explicit instantiation approach to keep some templates here.
// TODO: is this still needed now that we are a module?

template <class T>
bool tw::input::Numbers<T>::Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
{
    next_named_node(curs,false);
    if (num>1) {
        if (tw::input::node_kind(curs) == "tuple") {
            ts_tree_cursor_goto_first_child(curs);
        } else {
            tw::input::ThrowParsingError(curs,src,"expected tuple");
        }
    }
    for (tw::Int i=0;i<num;i++) {
        if (next_named_node(curs,true)) {
            dat[i] = tw::dnum(curs,src) >> native;
        } else {
            tw::input::ThrowParsingError(curs,src,"something missing");
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
        if (node_kind(curs)=="string_literal") {
            ts_tree_cursor_goto_first_child(curs);
            ts_tree_cursor_goto_next_sibling(curs);
            *str = node_text(curs,src);
            ts_tree_cursor_goto_parent(curs);
        } else if (node_kind(curs)=="identifier" || allow_cast) {
            *str = node_text(curs,src);
        } else {
            tw::input::ThrowParsingError(curs,src,"expected string or identifier");
        }
    } else {
        tw::input::ThrowParsingError(curs,src,"something missing");
    }
    return true;
}

bool tw::input::Bool::Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
{
    std::string word = next_named_node_text(curs,src);
    std::multimap<std::string,bool> m = {{"on",true},{"true",true},{"yes",true},{"off",false},{"false",false},{"no",false}};
    auto p = m.find(word);
    if (p==m.end()) {
        tw::input::ThrowParsingError(curs,src,"expected boolean");
    }
    *dat = p->second;
    return true;
}
