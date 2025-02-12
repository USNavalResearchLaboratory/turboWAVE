namespace tw
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
		inline std::string loc_str(const TSTreeCursor *curs) {
			TSNode n = ts_tree_cursor_current_node(curs);
			TSPoint p = ts_node_start_point(n);
			return std::to_string(p.row) + "," + std::to_string(p.column);
		}
		inline std::string missing(const TSTreeCursor *curs) {
			return "something missing at " + loc_str(curs);
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
		/// @brief goto the next named node, optionally accepting the current node
		/// @param curs tree cursor
		/// @param include_curr whether to accept the current node
		/// @return was a named node found
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
		void StripQuotes(std::string& str);
		/// @brief Directive is a base class for several child classes.
		/// Directive and its derivatives provide mainly the virtual function Read(...).
		/// Read(...) is invoked by the DirectiveReader class, which maintains a hash table of string-keys and pointers to Directive instances.
		/// This hash table is the basis for reading data from the input file.
		struct Directive
		{
			virtual ~Directive()
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
		struct Custom : Directive
		{
			virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
			{
				return false;
			}
		};
		struct String : Directive
		{
			std::string *str;
			String(std::string *s) { str=s; }
			virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
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
		};
		template <class T> // any type supported by input streams
		struct Numbers : Directive
		{
			T *dat;
			tw::Int num;
			Numbers(T *d,tw::Int n) { num=n; dat=d; }
			virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
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
		template <class T>
		struct List : Directive
		{
			// T is a standard container, or any class with resize, subscript, and value_type.
			T *dat;
			List(T *d) { dat=d; }
			virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
			{
				std::vector<typename T::value_type> temp;
				std::vector<typename T::value_type> dummy(1);
				if (next_named_node(curs,false)) {
					if (ts_node_type(ts_tree_cursor_current_node(curs))=="list") {
						if (ts_tree_cursor_goto_first_child(curs)) {
							next_named_node(curs,true);
							do {
								dummy[0] = tw::dnum(curs,src) >> native;
								temp.push_back(dummy[0]);
							} while (next_named_node(curs,false));
							ts_tree_cursor_goto_parent(curs);
						} else {
							throw tw::FatalError("empty list at " + loc_str(curs));
						}
					} else {
						throw tw::FatalError("expected list at " + loc_str(curs));
					}
				} else {
					throw tw::FatalError(missing(curs));
				}
				dat->resize(temp.size());
				for (tw::Int i=0;i<temp.size();i++)
					(*dat)[i] = temp[i];
				return true;
			}
		};
		template <class T>
		struct Enums : Directive
		{
			T *dat1,*dat2,*dat3;
			std::map<std::string,T> *emap; // has to be pointer for polymorphic destruction
			Enums(const std::map<std::string,T>& m,T* d1,T* d2=NULL,T* d3=NULL)
			{
				// must specify at least 1 data element, may specify up to 3.
				dat1=d1; dat2=d2; dat3=d3;
				emap = new std::map<std::string,T>(m);
			}
			virtual ~Enums()
			{
				delete emap;
			}
			virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
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
						if (emap->find(word)==emap->end())
							throw tw::FatalError("Invalid type <"+word+"> while reading key <"+key+">");
						*d = (*emap)[word];
						next_named_node(curs,false);
					}
				};
				read_one(dat1);
				read_one(dat2);
				read_one(dat3);
				return true;
			}
		};
		struct Bool : Directive
		{
			bool *dat;
			Bool(bool *d) { dat=d; }
			virtual bool Read(TSTreeCursor *curs,const std::string& src,const std::string& key,const tw::UnitConverter& native)
			{
				std::string word = next_named_node_text(curs,src);
				std::multimap<std::string,bool> m = {{"on",true},{"true",true},{"yes",true},{"off",false},{"false",false},{"no",false}};
				auto p = m.find(word);
				if (p==m.end())
					throw tw::FatalError("Invalid boolean value <"+word+"> while reading key <"+key+">");
				*dat = p->second;
				return true;
			}
		};
		class DirectiveReader
		{
			tw::Int maxKeyWords;
			std::vector<std::string> requiredKeys;
			std::map<std::string,tw::Int> keysFound;
			std::map<std::string,tw::input::Directive*> dmap;
			tw::UnitConverter native;
		public:
			DirectiveReader();
			~DirectiveReader();
			void AttachUnits(const tw::UnitConverter& native);
			void Reset();
			void Add(const std::string& key,tw::input::Directive *dir,bool required=true);
			bool ReadNext(const TSTreeCursor *curs,const std::string& src);
			void ReadAll(TSTreeCursor *curs,const std::string& src);
			bool TestKey(const std::string& test);
			void ThrowErrorIfMissingKeys(const std::string& src);
		};
		struct Preamble
		{
			// encapsulates information about how an object is being created in the input file
			bool attaching;
			std::string obj_key,obj_name,owner_name;
		};
		struct FileEnv
		{
			// encapsulates finding and opening files
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
