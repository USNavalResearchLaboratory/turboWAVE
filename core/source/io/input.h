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
		void StripQuotes(std::string& str);
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
