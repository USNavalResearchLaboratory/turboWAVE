namespace tw
{
	namespace input
	{
		struct Directive
		{
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				std::string word;
				in >> word;
				if (word!="=")
					throw tw::FatalError("Expected equals sign after key: "+key);
			}
		};
		struct Custom : Directive
		{
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				// Do not consume the equals sign or any other token.  Do nothing.
			}
		};
		template <class T> // any type supported by input streams
		struct Numbers : Directive
		{
			T *dat;
			tw::Int num;
			Numbers(T *d,tw::Int n) { num=n; dat=d; }
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				Directive::Read(in,key);
				for (tw::Int i=0;i<num;i++)
					in >> dat[i];
			}
		};
		template <class T> // any type with push_back method
		struct List : Directive
		{
			T *dat;
			List(T *d) { dat=d; }
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				Directive::Read(in,key);
				std::string word;
				in >> word;
				if (word!="{")
					throw tw::FatalError("Expected curly-brace at start of list.");
				do
				{
					in >> word;
					if (word!="}")
						dat->push_back(std::stod(word));
				} while (word!="}");
			}
		};
		template <class T>
		struct Enums : Directive
		{
			T *dat1,*dat2,*dat3;
			std::map<std::string,T> emap;
			Enums(std::map<std::string,T>& m,T* d1,T* d2=NULL,T* d3=NULL)
			{
				// must specify at least 1 data element, may specify up to 3.
				emap=m; dat1=d1; dat2=d2; dat3=d3;
			}
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				Directive::Read(in,key);
				std::string word;
				auto read_one = [&] (T *d)
				{
					if (d!=NULL)
					{
						in >> word;
						if (emap.find(word)==emap.end())
							throw tw::FatalError("Invalid value: "+word);
						*d = emap[word];
					}
				};
				read_one(dat1);
				read_one(dat2);
				read_one(dat3);
			}
		};
		struct Bool : Directive
		{
			bool *dat;
			Bool(bool *d) { dat=d; }
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				Directive::Read(in,key);
				std::string word;
				std::multimap<std::string,bool> m = {{"on",true},{"true",true},{"yes",true},{"off",false},{"false",false},{"no",false}};
				in >> word;
				auto p = m.find(word);
				if (p==m.end())
					throw tw::FatalError("Invalid boolean value: "+word);
				*dat = p->second;
			}
		};
		class DirectiveReader
		{
			tw::Int maxKeyWords;
			std::map<std::string,tw::Int> keysFound;
			std::map<std::string,tw::input::Directive*> dmap;
		public:
			DirectiveReader();
			~DirectiveReader();
			void Add(const std::string& key,tw::input::Directive *dir);
			std::string ReadNext(std::stringstream& in);
			void ReadAll(std::stringstream& in);
			bool TestKey(const std::string& test);
		};

		tw::Float GetUnitDensityCGS(std::stringstream& in);
		tw::Int IncludeFiles(std::stringstream& in,std::stringstream& out);
		void StripComments(std::ifstream& inputFile,std::stringstream& out);
		void StripDecorations(std::stringstream& in,std::stringstream& out);
		void InsertWhitespace(std::stringstream& in,std::stringstream& out);
		void UserMacros(std::stringstream& in,std::stringstream& out);
		void UnitMacros(std::stringstream& in,std::stringstream& out);
		void PreprocessInputFile(std::ifstream& inputFile,std::stringstream& out);

		void PythonRange(std::string& source,tw::Float *v0,tw::Float *v1);

		void ReadRect(Region *ans,std::stringstream& source);

		Profile* GetProfile(Simulation* sim,const std::string& name,const std::string& profileType);

		template <class T>
		void ReadArray(std::valarray<T>& data,std::stringstream& inputString);

		tw_boundary_spec ConvertBoundaryString(std::string& theString);

		void ReadBoundaryTerm(tw_boundary_spec *low,tw_boundary_spec *high,std::stringstream& theString,const std::string& command);

		void NormalizeInput(const UnitConverter& uc,std::string& in_out);

		bool GetQuotedString(std::string& str);

		std::vector<std::string> EnterInputFileBlock(std::stringstream& inputString,const std::string& end_tokens);

		void ExitInputFileBlock(std::stringstream& inputString,bool alreadyEntered);

		std::string GetPhrase(const std::vector<std::string>& words,tw::Int num_words);
	}
}

// Read curly braced sequence of numbers into "theArray"
// Expects equals sign and opening brace still on the input stringstream

template <class T>
void tw::input::ReadArray(std::valarray<T>& theArray,std::stringstream& inputString)
{
	std::vector<T> temp;
	std::string word;
	inputString >> word >> word;
	do
	{
		inputString >> word;
		if (word!="}")
			temp.push_back(std::stod(word,NULL));
	} while (word!="}");

	theArray.resize(temp.size());
	for (tw::Int i=0;i<temp.size();i++)
		theArray[i] = temp[i];
}
