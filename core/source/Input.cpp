#include "meta_base.h"
#include <regex>

////////////////////////
//                    //
//  Directive Reader  //
//                    //
////////////////////////

tw::input::DirectiveReader::DirectiveReader()
{
	maxKeyWords = 0;
	uc = NULL;
}

tw::input::DirectiveReader::~DirectiveReader()
{
	for (auto pair : dmap)
		delete pair.second;
}

void tw::input::DirectiveReader::AttachUnits(UnitConverter *uc)
{
	this->uc = uc;
}

void tw::input::DirectiveReader::Reset()
{
	keysFound.clear();
}

void tw::input::DirectiveReader::Add(const std::string& key,tw::input::Directive *dir,bool required)
{
	dmap[key] = dir;
	// count number of words in this key
	std::string word;
	std::stringstream s(key);
	tw::Int words = 1;
	while (!s.eof())
	{
		s >> word;
		if (!s.eof())
			words++;
	}
	// update max words
	if (words>maxKeyWords)
		maxKeyWords = words;
	// update required keys
	if (required)
		requiredKeys.push_back(key);
}

std::string tw::input::DirectiveReader::ReadNext(std::stringstream& in)
{
	std::string word,key;
	in >> key;
	if (in.eof())
		return "tw::EOF";
	if (key=="}" || key=="new" || key=="get" || key=="generate" || key=="open")
		return key;
	tw::Int word_count = 0;
	do
	{
		word_count++;
		if (dmap.find(key)!=dmap.end()) // key has been found
		{
			keysFound[key] = 0;
			if (uc==NULL)
				throw tw::FatalError("DirectiveReader class is missing units while processing key <"+key+">.");
			dmap[key]->Read(in,key,*uc);
			return key;
		}
		if (word_count<maxKeyWords)
		{
			in >> word;
			key += " " + word;
		}
	} while(word_count<maxKeyWords);
	throw tw::FatalError("Unexpected directive: <"+key+">.");
	return key;
}

void tw::input::DirectiveReader::ReadAll(std::stringstream& in)
{
	std::string s;
	do
	{
		s = ReadNext(in);
	} while(s!="}");
}

bool tw::input::DirectiveReader::TestKey(const std::string& test)
{
	return (keysFound.find(test)!=keysFound.end());
}

void tw::input::DirectiveReader::ThrowErrorIfMissingKeys(const std::string& src)
{
	for (auto s : requiredKeys)
		if (!TestKey(s))
			throw tw::FatalError("Missing required key <"+s+"> in <"+src+">");
}


///////////////////////////////////////
//                                   //
//  Stream Processing of Input File  //
//                                   //
///////////////////////////////////////

tw::input::FileEnv::FileEnv(const std::string& inputFilename)
{
	// Purpose of this class is mainly to encapsulate file searches
	this->inputFileName = inputFilename;
	searchPaths.push_back("");
	std::size_t sep_pos = inputFileName.find_last_of("/\\");
	searchPaths.push_back(inputFileName.substr(0,sep_pos+1)); // keep separator so we have the right one later
}

bool tw::input::FileEnv::OpenDeck(std::ifstream& inFile) const
{
	inFile.open(inputFileName);
	return inFile.good();
}

bool tw::input::FileEnv::FindAndOpen(const std::string& fileName,std::ifstream& inFile) const
{
	std::ifstream *temp;
	for (auto it=searchPaths.begin();it!=searchPaths.end();++it)
	{
		temp = new std::ifstream(*it + fileName);
		if (temp->good())
		{
			delete temp;
			inFile.open(*it + fileName);
			return inFile.good();
		}
		delete temp;
	}
	return false;
}

tw::Float tw::input::GetUnitDensityCGS(const std::string& in)
{
	// Following regex has 6 capture groups.
	// Group 2 is the mantissa, with 3,4,5 being alternatives within it.
	// Group 6 is either the exponent or white space.
	std::regex ex(R"((\bunit\s+density\s*=\s*)(([+-]?\d*\.\d+)|([+-]?\d+\.\d*)|([+-]?\d+))([eE][+-]?\d+|\s))");
	std::smatch match;
	if (std::regex_search(in,match,ex))
		return std::stod(match[2].str()+match[6].str(),NULL);
	else
		throw tw::FatalError("The <unit density> parameter is missing or ill-formed.");
}

tw::units tw::input::GetNativeUnits(const std::string& in)
{
	std::map<std::string,tw::units> m = tw::get_unit_map();
	std::regex ex(R"((\bnative\s+units\s*=\s*)(\w+))");
	std::smatch match;
	if (!std::regex_search(in,match,ex))
		return tw::units::plasma;
	else
		if (m.find(match[2].str())==m.end())
			throw tw::FatalError("Unkown system of units <"+match[2].str()+">.");
	return m[match[2].str()];
}

void tw::input::StripComments(std::ifstream& inFile,std::stringstream& out)
{
	bool ignoreUntilMark,ignoreUntilEOL;
	char charNow,nextChar;

	ignoreUntilMark = false;
	ignoreUntilEOL = false;
	while (inFile.get(charNow))
	{
		if (charNow=='/' || charNow=='*' || charNow=='\n' || charNow=='\r')
		{
			if (charNow=='/')
			{
				inFile.get(nextChar);
				if (nextChar=='/')
					ignoreUntilEOL = true;
				if (nextChar=='*' && !ignoreUntilEOL)
					ignoreUntilMark = true;
				if (nextChar!='/' && nextChar!='*')
				{
					if (!ignoreUntilMark && !ignoreUntilEOL)
						out << charNow;
					inFile.putback(nextChar);
				}
			}
			if (charNow=='*')
			{
				inFile.get(nextChar);
				if (nextChar=='/')
					ignoreUntilMark = false;
				else
				{
					if (!ignoreUntilMark && !ignoreUntilEOL)
						out << charNow;
					inFile.putback(nextChar);
				}
			}
			if (charNow=='\n' || charNow=='\r')
			{
				if (!ignoreUntilMark)
					out << charNow;
				ignoreUntilEOL = false;
			}
		}
		else
		{
			if (!ignoreUntilMark && !ignoreUntilEOL)
			{
				out << charNow;
			}
		}
	}
}

void tw::input::StripDecorations(std::stringstream& in,std::stringstream& out)
{
	// Replace commas and parenthesis with spaces
	char charNow;
	while (in.get(charNow))
	{
		if (charNow==',' || charNow=='(' || charNow==')')
			out << ' ';
		else
			out << charNow;
	}
}

void tw::input::InsertWhitespace(std::stringstream& in,std::stringstream& out)
{
	// Insert spaces around braces and equal signs
	char charNow;
	while (in.get(charNow))
	{
		if (charNow=='{' || charNow=='}' || charNow=='=')
			out << ' ';
		out << charNow;
		if (charNow=='{' || charNow=='}' || charNow=='=')
			out << ' ';
	}
}

void tw::input::UserMacros(std::stringstream& in,std::stringstream& out)
{
	std::map<std::string,std::string> macros;
	std::map<std::string,std::string>::iterator it;
	std::string word,key,val,stripped;
	std::stringstream temp;
	// First eat and save the key/value pairs
	// Currently this process destroys line feeds
	while (!in.eof())
	{
		in >> word;
		if (!in.eof())
		{
			if (word=="#define")
			{
				in >> key;
				it = macros.find(key);
				if (it==macros.end())
				{
					val = "";
					char next;
					do
					{
						next = in.get();
						if (next!='\n' && next!='\r' && !in.eof())
							val += next;
					} while(next!='\n' && next!='\r' && !in.eof());
					// eliminate leading white space
					std::string::size_type i = val.find_first_not_of(" \t");
					if (i!=std::string::npos)
						val = val.substr(i,std::string::npos);
					macros[key] = val;
				}
				else
					throw tw::FatalError("Macro "+key+" was already used.");
			}
			else
			{
				temp << word << " ";
			}
		}
	}
	// Now replace all the key matches with the values.
	// The key is only replaced if it is enclosed by white space.
	// N.b. by this time delimiters are replaced with white space.
	stripped = temp.str();
	for (auto pair : macros)
	{
		// first escape any special regex characters in the key
		std::regex special(R"([.^$|()\[\]{}*+?\\])");
		std::string escaped = std::regex_replace(pair.first,special,R"(\$&)");
		// now make the replacement
		std::regex ex(R"(\s)"+escaped+R"(\s)");
		stripped = std::regex_replace(stripped,ex," "+pair.second+" ");
	}
	// Finally put the result on the output stream
	out.str(stripped);
}

tw::Int tw::input::IncludeFiles(const FileEnv& file_env,std::stringstream& in,std::stringstream& out)
{
	tw::Int count = 0;
	std::string line_str,word;
	while (!in.eof())
	{
		std::getline(in,line_str);
		std::stringstream line(line_str);
		while (!line.eof())
		{
			if (line >> word)
			{
				if (word=="#include")
				{
					line >> word;
					StripQuotes(word);
					std::ifstream includedFile;
					if (!file_env.FindAndOpen(word,includedFile))
						throw tw::FatalError("couldn't open " + word);
					StripComments(includedFile,out);
					count++;
				}
				else
				{
					out << word << " ";
				}
			}
		}
		out << std::endl;
	}
	return count;
}

void tw::input::PreprocessInputFile(const FileEnv& file_env,std::stringstream& out)
{
	std::ifstream inFile;
	std::stringstream temp;
	std::string word;
	if (!file_env.OpenDeck(inFile))
		throw tw::FatalError("couldn't open input file " + file_env.inputFileName);

	auto reset_in = [&] (std::stringstream& ss)
	{
		ss.clear();
		ss.seekg(0,ss.beg);
	};

	auto reset_in_out = [&] (std::stringstream& i,std::stringstream& o)
	{
		i.clear();
		i.str(o.str());
		i.seekg(0,i.beg);
		o.clear();
		o.str("");
		o.seekp(0,o.beg);
	};

	// Handle included files, strip comments
	tw::input::StripComments(inFile,temp);
	reset_in(temp);
	while (tw::input::IncludeFiles(file_env,temp,out))
		reset_in_out(temp,out);
	// std::cout << out.str();
	// std::cout << std::endl << std::endl;

	// Clean formatting
	reset_in_out(temp,out);
	tw::input::StripDecorations(temp,out);
	// std::cout << out.str();
	// std::cout << std::endl << std::endl;
	reset_in_out(temp,out);
	tw::input::InsertWhitespace(temp,out);
	// std::cout << out.str();
	// std::cout << std::endl << std::endl;

	// User macro substitution
	reset_in_out(temp,out);
	tw::input::UserMacros(temp,out);
	// std::cout << out.str();
	// std::cout << std::endl << std::endl;
	inFile.close();
}


// Read a python.numpy style range, but it is a floating point range
// thus, blank resolves to 0 or big_pos

void tw::input::PythonRange(std::string& source,tw::Float *v0,tw::Float *v1)
{
	size_t colonPos;
	std::string sub;
	colonPos = source.find(':');

	if (colonPos==std::string::npos)
		throw tw::FatalError("Missing colon while reading python style range.");
	if (colonPos>0)
	{
		sub = source.substr(0,colonPos);
		try {
			*v0 = std::stod(sub,NULL); }
		catch (std::invalid_argument) {
			throw tw::FatalError("Invalid number while reading python style range."); }
	}
	else
		*v0 = 0.0;

	if (source.length() > colonPos+1)
	{
		sub = source.substr(colonPos+1,source.length()-colonPos-1);
		try {
			*v1 = std::stod(sub,NULL); }
		catch (std::invalid_argument) {
			throw tw::FatalError("Invalid number while reading python style range."); }
	}
	else
		*v1 = tw::big_pos;
}

void tw::input::StripQuotes(std::string& str)
{
	bool needed = false;
	if (str.front()=='\'' && str.back()=='\'')
		needed = true;
	if (str.front()=='\"' && str.back()=='\"')
		needed = true;
	if (needed)
	{
		str.pop_back();
		str = str.substr(1);
	}
}

tw::input::Preamble tw::input::EnterInputFileBlock(const std::string& com,std::stringstream& inputString,const std::string& end_tokens)
{
	tw::input::Preamble ans;
	std::string word;
	// Get the word list
	do
	{
		inputString >> word;
		if (end_tokens.find(word)==std::string::npos)
			ans.words.push_back(std::string(word));
	} while (end_tokens.find(word)==std::string::npos && !inputString.eof());
	ans.end_token = std::string(word);
	// Form the contiguous string representation
	ans.str = "";
	for (auto s : ans.words)
	{
		if (ans.str!="")
			ans.str = ans.str + " " + s;
		else
			ans.str = s;
	}
	ans.err_prefix = "While processing <"+ans.str+">: ";
	if (inputString.eof())
		throw tw::FatalError(ans.err_prefix+"encountered EOF.");
	// Check to see if the for keyword is present and in the right position
	if (ans.words.size()>2)
		ans.attaching = (ans.words[ans.words.size()-2]=="for" ? true : false);
	else
		ans.attaching = false;
	bool keywordFound = false;
	for (auto s : ans.words)
		if (s=="for")
			keywordFound = true;
	if (keywordFound && !ans.attaching)
		throw tw::FatalError(ans.err_prefix+"misplaced keyword <for>.");
	// Apply modifiers based on the command type
	if (com=="generate")
		ans.attaching = true;
	// Get the object names
	if (ans.attaching)
	{
		std::string option1,option2,testString;
		if (keywordFound)
		{
			// option1 is the word before <for>, which is taken if it is quoted
			// option2 is the concatenation all the words before <for>
			option1 = std::string(*(ans.words.end()-3));
			option2 = "";
			for (auto iter=ans.words.begin();iter<ans.words.end()-2;++iter)
				option2.append(*iter+"_");
			option2.pop_back();
		}
		else
		{
			// This is a generate block.
			// option1 is the second to last word, which is taken if it is quoted
			// option2 is the concatenation of all the words before the last word
			option1 = std::string(*(ans.words.end()-2));
			option2 = "";
			for (auto iter=ans.words.begin();iter<ans.words.end()-1;++iter)
				option2.append(*iter+"_");
			option2.pop_back();
		}
		testString = std::string(option1);
		tw::input::StripQuotes(testString);
		if (testString!=option1)
			ans.obj_name = std::string(option1);
		else
			ans.obj_name = std::string(option2);
		ans.owner_name = std::string(ans.words.back());
	}
	else
	{
		ans.obj_name = std::string(ans.words.back());
		ans.owner_name = "";
	}
	tw::input::StripQuotes(ans.obj_name);
	tw::input::StripQuotes(ans.owner_name);
	return ans;
}

std::string tw::input::GetPhrase(const std::vector<std::string>& words,tw::Int num_words)
{
	std::string ans("");
	tw::Int num = num_words<=words.size() ? num_words : words.size();
	for (tw::Int i=0;i<num;i++)
		ans += words[i] + " ";
	ans.pop_back();
	return ans;
}

void tw::input::ExitInputFileBlock(std::stringstream& inputString,bool alreadyEntered)
{
	std::string word;
	tw::Int leftCount=alreadyEntered?1:0;
	tw::Int rightCount=0;
	do
	{
		inputString >> word;
		if (word=="{")
			leftCount++;
		if (word=="}")
			rightCount++;
	} while (leftCount==0 || leftCount>rightCount);
}

void tw::input::PopExpectedWord(std::stringstream& inputString,const std::string& word,const std::string& obj)
{
	std::string found;
	inputString >> found;
	if (found!=word)
		throw tw::FatalError("Encountered <"+found+"> instead of <"+word+"> while reading <"+obj+">.");
}
