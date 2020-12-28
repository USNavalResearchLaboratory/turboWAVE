namespace tw
{
	namespace input
	{
		void StripQuotes(std::string& str);
		struct Directive
		{
			// Directive is a base class for several child classes.
			// Directive and its derivatives provide mainly the virtual function Read(...).
			// Read(...) is invoked by the DirectiveReader class, which maintains a hash table of string-keys and pointers to Directive instances.
			// This hash table is the basis for reading data from the input file.
			virtual ~Directive()
			{
			}
			virtual void Read(std::stringstream& in,const std::string& key,const UnitConverter& uc)
			{
				std::string word;
				in >> word;
				if (word!="=")
					throw tw::FatalError("Expected equals sign after key <"+key+">");
			}
		};
		struct Custom : Directive
		{
			virtual void Read(std::stringstream& in,const std::string& key,const UnitConverter& uc)
			{
				// Do not consume the equals sign or any other token.  Do nothing.
			}
		};
		struct String : Directive
		{
			std::string *str;
			String(std::string *s) { str=s; }
			virtual void Read(std::stringstream& in,const std::string& key,const UnitConverter& uc)
			{
				Directive::Read(in,key,uc);
				in >> *str;
				StripQuotes(*str);
			}
		};
		template <class T> // any type supported by input streams
		struct Numbers : Directive
		{
			T *dat;
			tw::Int num;
			Numbers(T *d,tw::Int n) { num=n; dat=d; }
			virtual void Read(std::stringstream& in,const std::string& key,const UnitConverter& uc)
			{
				Directive::Read(in,key,uc);
				for (tw::Int i=0;i<num;i++)
				{
					if (std::is_floating_point<T>::value)
					{
						tw::dnum dimensional_number;
						in >> dimensional_number;
						dat[i] = uc.ConvertToNative(dimensional_number);
					}
					else
					{
						if (!(in >> dat[i]))
							throw tw::FatalError("Invalid data type for key <"+key+">");
					}
				}
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
			// The contained type has to support input streams.
			T *dat;
			List(T *d) { dat=d; }
			virtual void Read(std::stringstream& in,const std::string& key,const UnitConverter& uc)
			{
				Directive::Read(in,key,uc);
				std::string word;
				in >> word;
				if (word!="{")
					throw tw::FatalError("Expected curly-brace at start of list.");
				// We are using streams to convert strings to data polymorphically.
				std::vector<typename T::value_type> temp;
				std::vector<typename T::value_type> dummy(1);
				while (!in.eof() && word!="}")
				{
					in >> word;
					if (!in.eof() && word!="}")
					{
						in.seekg(-word.size(),std::ios::cur);
						if (std::is_floating_point<typename T::value_type>::value)
						{
							tw::dnum dimensional_number;
							in >> dimensional_number;
							dummy[0] = uc.ConvertToNative(dimensional_number);
						}
						else
						{
							if (!(in >> dummy[0]))
								throw tw::FatalError("Invalid data type for key <"+key+">");
						}
						temp.push_back(dummy[0]);
					}
				}
				dat->resize(temp.size());
				for (tw::Int i=0;i<temp.size();i++)
					(*dat)[i] = temp[i];
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
			virtual void Read(std::stringstream& in,const std::string& key,const UnitConverter& uc)
			{
				Directive::Read(in,key,uc);
				std::string word;
				auto read_one = [&] (T *d)
				{
					if (d!=NULL)
					{
						in >> word;
						if (emap->find(word)==emap->end())
							throw tw::FatalError("Invalid type <"+word+"> while reading key <"+key+">");
						*d = (*emap)[word];
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
			virtual void Read(std::stringstream& in,const std::string& key,const UnitConverter& uc)
			{
				Directive::Read(in,key,uc);
				std::string word;
				std::multimap<std::string,bool> m = {{"on",true},{"true",true},{"yes",true},{"off",false},{"false",false},{"no",false}};
				in >> word;
				auto p = m.find(word);
				if (p==m.end())
					throw tw::FatalError("Invalid boolean value <"+word+"> while reading key <"+key+">");
				*dat = p->second;
			}
		};
		class DirectiveReader
		{
			tw::Int maxKeyWords;
			std::vector<std::string> requiredKeys;
			std::map<std::string,tw::Int> keysFound;
			std::map<std::string,tw::input::Directive*> dmap;
			UnitConverter *uc;
		public:
			DirectiveReader();
			~DirectiveReader();
			void AttachUnits(UnitConverter *uc);
			void Reset();
			void Add(const std::string& key,tw::input::Directive *dir,bool required=true);
			std::string ReadNext(std::stringstream& in);
			void ReadAll(std::stringstream& in);
			bool TestKey(const std::string& test);
			void ThrowErrorIfMissingKeys(const std::string& src);
		};
		struct Preamble
		{
			// encapsulates information about how an object is being created in the input file
			bool attaching;
			std::string str,obj_name,owner_name,end_token,err_prefix;
			std::vector<std::string> words;
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

		tw::Float GetUnitDensityCGS(const std::string& in);
		tw::units GetNativeUnits(const std::string& in);

		std::string include_regex();
		std::string define_key_regex();
		std::string define_regex();
		std::string ifdef_regex();
		std::string ifndef_regex();
		std::string else_regex();
		std::string endif_regex();
		void StripComments(std::string& in_out);
		void PreprocessorSyntaxCheck(const std::string& in);
		void IncludeFiles(const FileEnv& file_env,std::string& in_out);
		void AddMacro(const std::string& line,std::map<std::string,std::string>& macros);
		void EnterConditional(std::string& line,std::stringstream& in,std::stringstream& out,std::map<std::string,std::string>& macros);
		std::map<std::string,std::string> StripConditionalCode(std::string& in_out);
		void MacroSubstitution(std::string& in_out,const std::map<std::string,std::string>& macros);
		void PreprocessString(std::string& in_out);
		void PreprocessInputFile(const FileEnv& file_env,std::stringstream& out);

		Preamble EnterInputFileBlock(const std::string& com,std::stringstream& inputString,const std::string& end_tokens);
		void ExitInputFileBlock(std::stringstream& inputString,bool alreadyEntered);

		std::string GetPhrase(const std::vector<std::string>& words,tw::Int num_words);
		void PopExpectedWord(std::stringstream& inputString,const std::string& word,const std::string& obj);
		void PythonRange(std::string& source,tw::Float *v0,tw::Float *v1);
	}
}
