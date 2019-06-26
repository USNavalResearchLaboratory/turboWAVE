namespace tw
{
	namespace input
	{
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

		void ExitInputFileBlock(std::stringstream& inputString);

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
