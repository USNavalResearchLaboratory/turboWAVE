#include "meta_base.h"

/////////////////////////////
//                         //
//  Unit Conversion Tools  //
//                         //
/////////////////////////////

UnitConverter::UnitConverter(tw::Float unitDensityCGS)
{
	n1 = unitDensityCGS*1e6;
	wp = sqrt(n1*sqr(mks::qe)/(mks::eps0*mks::me));
	c = mks::c;
	qe = mks::qe;
	me = mks::me;
	eps0 = mks::eps0;
	kB = mks::kB;
	hbar = mks::hbar;
}

tw::Float UnitConverter::MKSValue(tw_dimensions dim) const
{
	switch (dim)
	{
		case time_dim:
			return 1.0/wp;
		case length_dim:
			return c/wp;
		case number_dim:
			return n1*cub(c/wp);
		case mass_dim:
			return me;
		case energy_dim:
			return me*c*c;
		case momentum_dim:
			return me*c;
		case angular_momentum_dim:
			return me*c*c/wp;
		case density_dim:
			return n1;
		case power_dim:
			return me*c*c*wp;
		case fluence_dim:
			return sqr(me*c*wp/qe)*c*eps0/wp; // E^2/eta0/wp
		case intensity_dim:
			return sqr(me*c*wp/qe)*c*eps0; // E^2/eta0
		case energy_density_dim:
			return me*c*c*n1;
		case power_density_dim:
			return me*c*c*n1*wp;
		case charge_dim:
			return qe;
		case current_dim:
			return n1*qe*c*sqr(c/wp);
		case current_density_dim:
			return n1*qe*c;
		case charge_density_dim:
			return n1*qe;
		case electric_field_dim:
			return me*c*wp/qe;
		case magnetic_field_dim:
			return me*wp/qe;
		case scalar_potential_dim:
			return me*c*c/qe;
		case diffusivity_dim:
			return c*c/wp;
		case conductivity_dim:
			return (n1*qe*c)/(me*c*wp/qe); // j/E
		case rate_coefficient_2_dim:
			return wp/n1;
		case rate_coefficient_3_dim:
			return wp/(n1*n1);
		case mobility_dim:
			return c/(me*c*wp/qe); // c/E
		case temperature_dim:
			return me*c*c/kB;
		case cross_section_dim:
			return wp/(n1*c);
		case susceptibility_dim:
			return 1.0;
		case susceptibility_2_dim:
			return 1.0/(me*c*wp/qe); // 1/E
		case susceptibility_3_dim:
			return 1.0/sqr(me*c*wp/qe); // 1/E^2
		default:
			return 0.0;
	}
}

tw::Float UnitConverter::CGSValue(tw_dimensions dim) const
{
	switch (dim)
	{
		case time_dim:
			return 1.0/wp;
		case length_dim:
			return 1e2*c/wp;
		case number_dim:
			return n1*cub(c/wp);
		case mass_dim:
			return 1e3*me;
		case energy_dim:
			return 1e7*me*c*c;
		case momentum_dim:
			return 1e5*me*c;
		case angular_momentum_dim:
			return 1e7*me*c*c/wp;
		case density_dim:
			return 1e-6*n1;
		case power_dim:
			return 1e7*me*c*c*wp;
		case fluence_dim:
			return 1e3*sqr(me*c*wp/qe)*c*eps0/wp; // E^2/eta0/wp
		case intensity_dim:
			return 1e3*sqr(me*c*wp/qe)*c*eps0; // E^2/eta0
		case energy_density_dim:
			return 1e1*me*c*c*n1;
		case power_density_dim:
			return 1e1*me*c*c*n1*wp;
		case charge_dim:
			return 3e9*qe;
		case current_dim:
			return 3e9*n1*qe*c*sqr(c/wp);
		case current_density_dim:
			return 3e5*n1*qe*c;
		case charge_density_dim:
			return 3e3*n1*qe;
		case electric_field_dim:
			return 0.333333e-4*me*c*wp/qe;
		case magnetic_field_dim:
			return 1e4*me*wp/qe;
		case scalar_potential_dim:
			return 0.333333e-2*me*c*c/qe;
		case diffusivity_dim:
			return 1e4*c*c/wp;
		case conductivity_dim:
			return 9e9*(n1*qe*c)/(me*c*wp/qe); // j/E
		case rate_coefficient_2_dim:
			return 1e6*wp/n1;
		case rate_coefficient_3_dim:
			return 1e12*wp/(n1*n1);
		case mobility_dim:
			return 1e2*c/(0.333333e-4*me*c*wp/qe); // c/E
		case temperature_dim:
			return me*c*c/kB;
		case cross_section_dim:
			return 1e4*wp/(n1*c);
		case susceptibility_dim:
			return 1.0/(4.0*pi);
		case susceptibility_2_dim:
			return 1.0/(4.0*pi*0.333333e-4*me*c*wp/qe); // 1/4*pi*E
		case susceptibility_3_dim:
			return 1.0/sqr(4.0*pi*0.333333e-4*me*c*wp/qe); // 1/4*pi*E^2
		default:
			return 0.0;
	}
}


////////////////////////
//                    //
//  Directive Reader  //
//                    //
////////////////////////

tw::input::DirectiveReader::DirectiveReader()
{
	maxKeyWords = 0;
}

tw::input::DirectiveReader::~DirectiveReader()
{
	for (auto it=dmap.begin();it!=dmap.end();++it)
		delete it->second;
}

void tw::input::DirectiveReader::Add(const std::string& key,tw::input::Directive *dir)
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
}

std::string tw::input::DirectiveReader::ReadNext(std::stringstream& in)
{
	std::string word,key;
	in >> key;
	if (key=="}" || key=="new" || key=="get")
		return key;
	tw::Int word_count = 0;
	do
	{
		word_count++;
		if (dmap.find(key)!=dmap.end()) // key has been found
		{
			keysFound[key] = 0;
			dmap[key]->Read(in,key);
			return key;
		}
		if (word_count<maxKeyWords)
		{
			in >> word;
			key += " " + word;
		}
	} while(word_count<maxKeyWords);
	throw tw::FatalError("Unexpected directive: "+key);
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


///////////////////////////////////////
//                                   //
//  Stream Processing of Input File  //
//                                   //
///////////////////////////////////////

tw::Float tw::input::GetUnitDensityCGS(std::stringstream& in)
{
	tw::Float ans = 0.0;
	std::string word;
	while (!in.eof())
	{
		in >> word;
		if (word=="unit")
		{
			in >> word;
			if (word=="density")
				in >> word >> ans;
		}
	}
	return ans;
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
	std::string word,key,val;
	while (!in.eof())
	{
		in >> word;
		if (!in.eof())
		{
			if (word=="#define")
			{
				in >> key >> val;
				it = macros.find(key);
				if (it==macros.end())
					macros[key] = val;
				else
					throw tw::FatalError("Macro "+key+" was already used.");
			}
			else
			{
				char unaryOp = ' ';
				key = word;
				it = macros.find(key);
				if (it==macros.end())
				{
					// If not a direct macro, check to see if there is a unary operator.
					// At present negation is the only one.
					if (word[0]=='-')
					{
						unaryOp = '-';
						key = word.substr(1,std::string::npos);
						it = macros.find(key);
					}
				}
				if (it==macros.end())
					out << word << " ";
				else
					out << unaryOp << macros[key] << " ";
			}
		}
	}
}

void tw::input::UnitMacros(std::stringstream& in,std::stringstream& out)
{
	bool unitMacrosPresent = false;
	std::string word;

	while (!in.eof())
	{
		in >> word;
		if (!in.eof())
			if (word[0]=='%')
				unitMacrosPresent = true;
	}

	in.clear();
	in.seekg(0,in.beg);
	tw::Float unitDensityCGS = tw::input::GetUnitDensityCGS(in);
	if (unitDensityCGS==0.0 && unitMacrosPresent)
		throw tw::FatalError("Unit density directive is required.");
	UnitConverter uc(unitDensityCGS);

	in.clear();
	in.seekg(0,in.beg);
	while (!in.eof())
	{
		in >> word;
		if (!in.eof())
		{
			if (word[0]=='%')
				tw::input::NormalizeInput(uc,word);
			out << word << " ";
		}
	}
}

tw::Int tw::input::IncludeFiles(std::stringstream& in,std::stringstream& out)
{
	tw::Int count = 0;
	std::ifstream *includedFile;
	std::string word;
	while (!in.eof())
	{
		in >> word;
		if (!in.eof())
		{
			if (word=="#include")
			{
				in >> word;
				includedFile = new std::ifstream(word.c_str());
				if (includedFile->rdstate() & std::ios::failbit)
					throw tw::FatalError("couldn't open " + word);
				StripComments(*includedFile,out);
				delete includedFile;
				count++;
			}
			else
			{
				out << word << " ";
			}
		}
	}
	return count;
}

void tw::input::PreprocessInputFile(std::ifstream& inFile,std::stringstream& out)
{
	std::stringstream temp;
	std::string word;
	tw::Float unitDensityCGS;

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
	while (tw::input::IncludeFiles(temp,out))
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

	// Unit conversion macro substitution
	reset_in_out(temp,out);
	tw::input::UnitMacros(temp,out);
	// std::cout << out.str();
	// std::cout << std::endl << std::endl;
}


// Read a python.numpy style range, but it is a floating point range
// thus, blank resolves to 0 or big_pos

void tw::input::PythonRange(std::string& source,tw::Float *v0,tw::Float *v1)
{
	size_t colonPos;
	std::string sub;
	colonPos = source.find(':');

	if (colonPos>0)
	{
		sub = source.substr(0,colonPos);
		*v0 = std::stod(sub,NULL);
	}
	else
		*v0 = 0.0;

	if (source.length() > colonPos+1)
	{
		sub = source.substr(colonPos+1,source.length()-colonPos-1);
		*v1 = std::stod(sub,NULL);
	}
	else
		*v1 = tw::big_pos;
}

// Read boundary conditions

tw_boundary_spec tw::input::ConvertBoundaryString(std::string& theString)
{
	if (theString=="periodic" || theString=="cyclic")
		return cyclic;
	if (theString=="reflective" || theString=="reflecting")
		return reflecting;
	if (theString=="absorbing" || theString=="open")
		return absorbing;
	if (theString=="emitting" || theString=="emissive")
		return emitting;
	if (theString=="axisymmetric" || theString=="axisymmetry")
		return axisymmetric;

	return cyclic;
}

void tw::input::ReadBoundaryTerm(tw_boundary_spec *low,tw_boundary_spec *high,std::stringstream& theString,const std::string& command)
{
	std::string word;
	tw::Int axis = 0;
	if (command=="xboundary") axis = 1;
	if (command=="yboundary") axis = 2;
	if (command=="zboundary") axis = 3;
	if (axis>0)
	{
		theString >> word >> word;
		low[axis] = ConvertBoundaryString(word);
		theString >> word;
		high[axis] = ConvertBoundaryString(word);
	}
}

void tw::input::NormalizeInput(const UnitConverter& uc,std::string& in_out)
{
	// To be used in preprocessing input file
	// un-normalized quantities are signaled with percent,e.g. %1.0m
	// signals to replace the string with normalized length corresponding to 1 meter
	std::size_t endpos;
	std::string units;
	std::stringstream temp;
	tw::Float qty,nqty=tw::small_pos;
	if (in_out[0]=='%')
	{
		try { qty = std::stod(in_out.substr(1),&endpos); }
		catch (std::invalid_argument) { throw tw::FatalError("Invalid unit conversion macro : " + in_out); }
		units = in_out.substr(1+endpos);

		// Angle Conversions
		if (units=="rad")
			nqty = qty;
		if (units=="mrad")
			nqty = 1e-3*qty;
		if (units=="urad")
			nqty = 1e-6*qty;
		if (units=="deg")
			nqty = qty*pi/180.0;

		// Length Conversions
		if (units=="um")
			nqty = uc.MKSToSim(length_dim,1e-6*qty);
		if (units=="mm")
			nqty = uc.MKSToSim(length_dim,1e-3*qty);
		if (units=="cm")
			nqty = uc.MKSToSim(length_dim,1e-2*qty);
		if (units=="m")
			nqty = uc.MKSToSim(length_dim,qty);

		// Time Conversions
		if (units=="fs")
			nqty = uc.MKSToSim(time_dim,1e-15*qty);
		if (units=="ps")
			nqty = uc.MKSToSim(time_dim,1e-12*qty);
		if (units=="ns")
			nqty = uc.MKSToSim(time_dim,1e-9*qty);
		if (units=="us")
			nqty = uc.MKSToSim(time_dim,1e-6*qty);
		if (units=="s")
			nqty = uc.MKSToSim(time_dim,qty);

		// Number Density Conversions
		if (units=="m-3")
			nqty = uc.MKSToSim(density_dim,qty);
		if (units=="cm-3")
			nqty = uc.CGSToSim(density_dim,qty);

		// Energy Density Conversions
		if (units=="Jm3")
			nqty = uc.MKSToSim(energy_density_dim,qty);
		if (units=="Jcm3")
			nqty = uc.MKSToSim(energy_density_dim,1e6*qty);

		// Temperature Conversions
		if (units=="eV")
			nqty = uc.eV_to_sim(qty);
		if (units=="K")
			nqty = uc.MKSToSim(temperature_dim,qty);

		// Cross Section Conversions
		if (units=="cm2")
			nqty = uc.CGSToSim(cross_section_dim,qty);
		if (units=="m2")
			nqty = uc.MKSToSim(cross_section_dim,qty);

		// Diffusivity Conversions
		if (units=="cm2s")
			nqty = uc.CGSToSim(diffusivity_dim,qty);
		if (units=="m2s")
			nqty = uc.MKSToSim(diffusivity_dim,qty);

		if (nqty==tw::small_pos)
			throw tw::FatalError("Unrecognized Units " + in_out);
		else
		{
			temp << nqty; // don't use std::to_string; it rounds small numbers to 0.
			in_out = temp.str();
		}
	}
}

bool tw::input::GetQuotedString(std::string& str)
{
	bool isLeftQuote = str.front()=='\'' || str.front()=='\"';
	bool isRightQuote = str.back()=='\'' || str.back()=='\"';
	if (isLeftQuote && isRightQuote)
	{
		str.pop_back();
		str = str.substr(1);
		return true;
	}
	return false;
}

std::vector<std::string> tw::input::EnterInputFileBlock(std::stringstream& inputString,const std::string& end_tokens)
{
	std::string word;
	std::vector<std::string> preamble;
	do
	{
		inputString >> word;
		if (end_tokens.find(word)==std::string::npos)
			preamble.push_back(std::string(word));
	} while (end_tokens.find(word)==std::string::npos && !inputString.eof());
	if (inputString.eof())
		throw tw::FatalError("Encountered EOF while processing key <"+preamble[0]+">");
	return preamble;
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
