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


////////////////////////////////////
//                                //
//  Pre-Processing of Input File  //
//                                //
////////////////////////////////////

tw::input::FileEnv::FileEnv(const std::string& inputFilename)
{
	// Purpose of this class is mainly to encapsulate file searches
	this->inputFileName = inputFilename;
	searchPaths.push_back("");
	std::size_t sep_pos = inputFileName.find_last_of("/\\");
	searchPaths.push_back(inputFileName.substr(0,sep_pos+1)); // keep separator so we have the right one later
}

bool tw::input::FileEnv::OpenDeck(std::string& contents) const
{
	std::ifstream inFile;
	std::stringstream temp;
	inFile.open(inputFileName);
	if (inFile.good())
	{
		temp << inFile.rdbuf();
		contents = temp.str();
		inFile.close();
		return true;
	}
	return false;
}

bool tw::input::FileEnv::FindAndOpen(const std::string& fileName,std::string& contents) const
{
	std::ifstream *inFile;
	std::stringstream temp;
	for (auto it=searchPaths.begin();it!=searchPaths.end();++it)
	{
		inFile = new std::ifstream(*it + fileName);
		if (inFile->good())
		{
			temp << inFile->rdbuf();
			contents = temp.str();
			inFile->close();
			delete inFile;
			return true;
		}
		delete inFile;
	}
	return false;
}

bool tw::input::FileEnv::FindAndOpen(const std::string& fileName,std::stringstream& contents) const
{
	std::ifstream *inFile;
	for (auto it=searchPaths.begin();it!=searchPaths.end();++it)
	{
		inFile = new std::ifstream(*it + fileName);
		if (inFile->good())
		{
			contents << inFile->rdbuf();
			inFile->close();
			delete inFile;
			return true;
		}
		delete inFile;
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

void tw::input::StripComments(std::string& in_out)
{
	// Strip comments : does not preserve number of lines
	// The following blindly matches any comment, whether quoted or not.
	// TW has always used the blind version, continue with that for now.
	std::regex ex_blind(R"(\/\/.*?(\n|\r\n)|\/\*(.|\s)*?\*\/)");
	try
	{
		in_out = std::regex_replace(in_out,ex_blind,std::string("\n"));
	}
	catch (...)
	{
		throw tw::FatalError("Error trying to strip comments.  Perhaps C-style comment wasn't closed.");
	}
}

// The following regex are intended to be processed line by line.
std::string tw::input::include_regex()
{
	return R"(([ \t]*#include[ \t]+)([^"']\S+[^"' \t]|"\S*"|'\S*')([ \t]*))";
}
std::string tw::input::define_key_regex()
{
	return R"(([^-+,=(){} \t][^,=(){} \t]+))";
}
std::string tw::input::define_regex()
{
	// group 2 is the key, group 5 is the value (if present)
	return R"(([ \t]*#define[ \t]+))" + define_key_regex() + R"((([ \t]*)(.*)))";
}
std::string tw::input::ifdef_regex()
{
	return R"(([ \t]*#ifdef[ \t]+))" + define_key_regex() + R"(([ \t]*))";
}
std::string tw::input::ifndef_regex()
{
	return R"(([ \t]*#ifndef[ \t]+))" + define_key_regex() + R"(([ \t]*))";
}
std::string tw::input::else_regex()
{
	return R"([ \t]*#else[ \t]*)";
}
std::string tw::input::endif_regex()
{
	return R"([ \t]*#endif[ \t]*)";
}

void tw::input::PreprocessorSyntaxCheck(const std::string& in)
{
	std::regex appearance(R"(.*(#include|#define|#ifdef|#ifndef|#else|#endif).*)");
	std::smatch match;
	std::stringstream temp(in);
	std::string line;
	tw::Int if_depth=0,else_depth=0;
	while (!temp.eof())
	{
		std::getline(temp,line);
		if (!temp.eof())
		{
			if (std::regex_search(line,appearance))
			{
				// First just check syntax of this line only
				bool good_inc = std::regex_match(line,std::regex(include_regex()));
				bool good_def = std::regex_match(line,std::regex(define_regex()));
				bool good_ifdef = std::regex_match(line,std::regex(ifdef_regex()));
				bool good_ifndef = std::regex_match(line,std::regex(ifndef_regex()));
				bool good_else = std::regex_match(line,std::regex(else_regex()));
				bool good_endif = std::regex_match(line,std::regex(endif_regex()));
				// Then if everything is OK check nesting, otherwise throw error
				if (good_inc || good_def || good_ifdef || good_ifndef || good_else || good_endif)
				{
					if (good_ifdef || good_ifndef)
						if_depth++;
					if (good_else)
						else_depth++;
					if (good_endif)
					{
						if_depth--;
						if (else_depth)
							else_depth--;
					}
					if (if_depth<0)
						throw tw::FatalError("There is an extra #endif");
					if (else_depth>if_depth)
						throw tw::FatalError("There is an extra #else");
				}
				else
					throw tw::FatalError("Preprocessor directive badly formed: "+line);
			}
		}
	}
	if (if_depth)
		throw tw::FatalError("Missing " + std::to_string(if_depth) + " #endif " + (if_depth==1?"level":"levels"));
}

void tw::input::AddMacro(const std::string& line,std::map<std::string,std::string>& macros)
{
	std::string key,val;
	std::smatch match;
	if (std::regex_match(line,match,std::regex(define_regex())))
	{
		key = match[2].str();
		val = match[5].str();
		if (key==val)
			throw tw::FatalError("Macro cannot refer to itself ("+key+")");
		if (macros.find(key)==macros.end())
			macros[key] = val;
		else
			throw tw::FatalError("Macro "+key+" was already used.");
	}
}

void tw::input::EnterConditional(std::string& line,std::stringstream& in,std::stringstream& out,std::map<std::string,std::string>& macros)
{
	// Keep or clear lines until matching #endif is encountered.
	// Continue to build macros as we go.
	// Call this recursively if further conditionals encountered.
	bool keeping,is_ifdef,is_ifndef,is_else,is_endif;
	tw::Int depth=1;
	std::smatch match;
	out << std::endl; // clear the conditional that got us here
	if (std::regex_match(line,match,std::regex(ifdef_regex())))
		keeping = macros.find(match[2].str())!=macros.end();
	if (std::regex_match(line,match,std::regex(ifndef_regex()))) // can't use else, must generate the match
		keeping = macros.find(match[2].str())==macros.end();
	do
	{
		std::getline(in,line);
		is_ifdef = std::regex_match(line,std::regex(ifdef_regex()));
		is_ifndef = std::regex_match(line,std::regex(ifndef_regex()));
		is_else = std::regex_match(line,std::regex(else_regex()));
		is_endif = std::regex_match(line,std::regex(endif_regex()));
		if (is_else)
		{
			keeping = !keeping;
			out << std::endl;
		}
		else if (is_endif)
		{
			depth--;
			out << std::endl;
		}
		else
		{
			if (keeping)
			{
				AddMacro(line,macros);
				if (is_ifdef || is_ifndef)
					EnterConditional(line,in,out,macros);
				else
					out << line << std::endl;
			}
			else
			{
				if (is_ifdef || is_ifndef)
					depth++;
				out << std::endl;
			}
		}
	} while (depth);
}

std::map<std::string,std::string> tw::input::StripConditionalCode(std::string& in_out)
{
	// Handle #ifdef and #ifndef blocks
	// Assume syntax has already been checked
	// This routine also loads the map of #define constants
	std::map<std::string,std::string> macros;
	std::string line;
	std::stringstream in(in_out),out;

	while (!in.eof())
	{
		std::getline(in,line);
		if (!in.eof())
		{
			AddMacro(line,macros);
			if (std::regex_match(line,std::regex(ifdef_regex())) || std::regex_match(line,std::regex(ifndef_regex())))
				EnterConditional(line,in,out,macros);
			else
				out << line << std::endl;
		}
	}
	in_out = out.str();
	return macros;
}

void tw::input::MacroSubstitution(std::string& in_out,const std::map<std::string,std::string>& macros)
{
	// assumes syntax has already been checked.
	// assumes conditionals have already been processed.
	std::string line;
	std::stringstream in(in_out),temp;
	bool found;

	// First eat #define lines
	while (!in.eof())
	{
		std::getline(in,line);
		if (!in.eof())
		{
			if (std::regex_match(line,std::regex(define_regex())))
				temp << std::endl;
			else
				temp << line << std::endl;
		}
	}
	in_out = temp.str();

	// Now replace all the key matches with the values.
	// Do multiple sweeps to account for macros within macros.
	do
	{
		found = false;
		for (auto pair : macros)
		{
			// first escape any special regex characters in the key
			std::regex special(R"([.^$|()\[\]{}*+?\\])");
			std::string escaped = std::regex_replace(pair.first,special,R"(\$&)");
			// now make the replacement (could be simpler if we had lookbehind)
			std::regex ex(R"(([-+,=(){}\s]|^))" + escaped + R"((?=[,=(){}\s]|$))");
			// must beware that pair.second very often starts with a digit!
			found = found || std::regex_search(in_out,ex);
			in_out = std::regex_replace(in_out,ex,R"($01)"+pair.second);
		}
	} while (found);
}

void tw::input::IncludeFiles(const FileEnv& file_env,std::string& in_out)
{
	tw::Int count;
	std::string line,filename;
	std::smatch match;
	std::stringstream in,out,temp;
	do
	{
		count = 0;
		in.clear();
		in.str(in_out);
		in.seekg(0,in.beg);
		out.clear();
		out.seekp(0,out.beg);
		while (!in.eof())
		{
			std::getline(in,line);
			if (!in.eof())
			{
				if (std::regex_match(line,match,std::regex(include_regex())))
				{
					count++;
					filename = match[2].str();
					StripQuotes(filename);
					if (!file_env.FindAndOpen(filename,line))
						throw tw::FatalError("couldn't open " + filename);
					PreprocessString(line);
				}
				out << line << std::endl;
			}
		}
		in_out = out.str();
	} while (count);
}

void tw::input::PreprocessString(std::string& in_out)
{
	std::map<std::string,std::string> macros; // treat as local to a file
	tw::input::StripComments(in_out);
	tw::input::PreprocessorSyntaxCheck(in_out);
	macros = tw::input::StripConditionalCode(in_out);
	tw::input::MacroSubstitution(in_out,macros);
	// std::cout << in_out << std::endl << std::endl;
}

void tw::input::PreprocessInputFile(const FileEnv& file_env,std::stringstream& out)
{
	std::string deck;

	if (!file_env.OpenDeck(deck))
		throw tw::FatalError("couldn't open input file " + file_env.inputFileName);

	tw::input::PreprocessString(deck);
	tw::input::IncludeFiles(file_env,deck);
	deck = std::regex_replace(deck,std::regex(R"([\,\(\)])")," ");
	deck = std::regex_replace(deck,std::regex(R"([\=\{\}])")," $& ");

	out.str(deck);
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
		if (s=="new" || s=="generate" || s=="get")
			throw tw::FatalError(ans.err_prefix+"misplaced keyword <"+s+">.");
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
