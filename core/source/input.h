namespace tw
{
	namespace input
	{
		// Directive is a base class for several child classes.
		// Directive and its derivatives provide mainly the virtual function Read(...).
		// Read(...) is invoked by the DirectiveReader class, which maintains a hash table of string-keys and pointers to Directive instances.
		// This hash table is the basis for reading data from the input file.
		void StripQuotes(std::string& str);
		struct Directive
		{
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				std::string word;
				in >> word;
				if (word!="=")
					throw tw::FatalError("Expected equals sign after key <"+key+">");
			}
		};
		struct Custom : Directive
		{
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				// Do not consume the equals sign or any other token.  Do nothing.
			}
		};
		struct String : Directive
		{
			std::string *str;
			String(std::string *s) { str=s; }
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				Directive::Read(in,key);
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
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				Directive::Read(in,key);
				for (tw::Int i=0;i<num;i++)
					if (!(in >> dat[i]))
						throw tw::FatalError("Invalid number while reading data for key <"+key+">");
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
			// T has to support the resize method, and be subscriptable.
			// The contained type has to support input streams.
			T *dat;
			List(T *d) { dat=d; }
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				Directive::Read(in,key);
				std::string word;
				in >> word;
				if (word!="{")
					throw tw::FatalError("Expected curly-brace at start of list.");
				// First put the string-data on a temporary stream and count.
				// We are using streams to convert strings to data polymorphically.
				std::stringstream temp;
				tw::Int count = 0;
				do
				{
					in >> word;
					if (word!="}")
					{
						count++;
						temp << word << " ";
					}
				} while (word!="}");
				temp.seekg(0,temp.beg);
				dat->resize(count);
				for (tw::Int i=0;i<count;i++)
					if (!(temp >> (*dat)[i]))
						throw tw::FatalError("Invalid data type for key <"+key+">");
			}
		};
		template <class T>
		struct Enums : Directive
		{
			T *dat1,*dat2,*dat3;
			std::map<std::string,T> emap;
			Enums(const std::map<std::string,T>& m,T* d1,T* d2=NULL,T* d3=NULL) : emap(m)
			{
				// must specify at least 1 data element, may specify up to 3.
				dat1=d1; dat2=d2; dat3=d3;
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
							throw tw::FatalError("Invalid type <"+word+"> while reading key <"+key+">");
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
		public:
			DirectiveReader();
			~DirectiveReader();
			void Reset();
			void Add(const std::string& key,tw::input::Directive *dir,bool required=true);
			std::string ReadNext(std::stringstream& in);
			void ReadAll(std::stringstream& in);
			bool TestKey(const std::string& test);
			void ThrowErrorIfMissingKeys(const std::string& src);
		};
		struct Preamble
		{
			bool attaching;
			std::string str,obj_name,owner_name,end_token,err_prefix;
			std::vector<std::string> words;
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

		void NormalizeInput(const UnitConverter& uc,std::string& in_out);

		Preamble EnterInputFileBlock(const std::string& com,std::stringstream& inputString,const std::string& end_tokens);

		void ExitInputFileBlock(std::stringstream& inputString,bool alreadyEntered);

		std::string GetPhrase(const std::vector<std::string>& words,tw::Int num_words);

		void PopExpectedWord(std::stringstream& inputString,const std::string& word,const std::string& obj);
	}
}
