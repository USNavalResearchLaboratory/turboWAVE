enum tw_dimensions
{
	time_dim,
	length_dim,
	number_dim,
	mass_dim,
	energy_dim,
	momentum_dim,
	angular_momentum_dim,
	density_dim,
	power_dim,
	fluence_dim,
	intensity_dim,
	energy_density_dim,
	power_density_dim,
	charge_dim,
	current_dim,
	current_density_dim,
	charge_density_dim,
	electric_field_dim,
	magnetic_field_dim,
	scalar_potential_dim,
	diffusivity_dim,
	conductivity_dim,
	rate_coefficient_2_dim,
	rate_coefficient_3_dim,
	mobility_dim,
	temperature_dim,
	cross_section_dim,
	susceptibility_dim,
	susceptibility_2_dim,
	susceptibility_3_dim
};

struct UnitConverter
{
	tw::Float n1,wp;
	tw::Float c,qe,me,eps0,kB,hbar;

	UnitConverter(tw::Float unitDensityCGS);
	tw::Float MKSValue(tw_dimensions dim) const;
	tw::Float CGSValue(tw_dimensions dim) const;

	tw::Float SimToMKS(tw_dimensions dim,tw::Float val) const
	{
		return val*MKSValue(dim);
	}
	tw::Float MKSToSim(tw_dimensions dim,tw::Float val) const
	{
		return val/MKSValue(dim);
	}
	tw::Float SimToCGS(tw_dimensions dim,tw::Float val) const
	{
		return val*CGSValue(dim);
	}
	tw::Float CGSToSim(tw_dimensions dim,tw::Float val) const
	{
		return val/CGSValue(dim);
	}
	tw::Float sim_to_eV(tw::Float val) const
	{
		return val*MKSValue(energy_dim)/qe;
	}
	tw::Float eV_to_sim(tw::Float val) const
	{
		return val*qe/MKSValue(energy_dim);
	}
};

namespace tw
{
	namespace input
	{
		// Directive is a base class for several child classes.
		// Directive and its derivatives provide mainly the virtual function Read(...).
		// Read(...) is invoked by the DirectiveReader class, which maintains a hash table of string-keys and pointers to Directive instances.
		// This hash table is the basis for reading data from the input file.
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
		struct Vec3 : Numbers<tw::Float>
		{
			Vec3(tw::vec3 *v) : Numbers<tw::Float>(&v->x,3) { ; }
		};
		template <class T>
		struct List : Directive
		{
			// This is designed to work only for T being a container of tw::Float types.
			// The container has to support the resize method, and be subscriptable.
			T *dat;
			List(T *d) { dat=d; }
			virtual void Read(std::stringstream& in,const std::string& key)
			{
				std::vector<tw::Float> temp;
				Directive::Read(in,key);
				std::string word;
				in >> word;
				if (word!="{")
					throw tw::FatalError("Expected curly-brace at start of list.");
				do
				{
					in >> word;
					if (word!="}")
						temp.push_back(std::stod(word));
				} while (word!="}");
				dat->resize(temp.size());
				for (tw::Int i=0;i<temp.size();i++)
					(*dat)[i] = temp[i];
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
