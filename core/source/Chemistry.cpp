#include "meta_base.h"
#include "computeTool.h"
#include "physics.h"
#include "chemistry.h"

tw::Float PrimitiveReaction::PrimitiveRate(tw::Float T)
{
	tw::Float rate;

	if (T<T0 || T>T1)
		return 0.0;

	if (c1==0.0) // use janev form
	{
		const tw::Float logT_eV = log(tw::small_pos + fabs(T)*unit_T_eV);
		rate = 0.0;
		for (tw::Int s=0;s<9;s++)
			rate += b[s]*pow(logT_eV,tw::Float(s));
		rate = exp(rate) / unit_rate_cgs;
	}
	else // use arrhenius form
		rate = c1 * pow(T,c2) * exp(-c3/T);

	return rate;
}

void PrimitiveReaction::ReadRate(std::stringstream& inputString,tw::Int numBodies,UnitConverter& uc)
{
	unit_T_eV = uc.sim_to_eV(1.0);
	unit_rate_cgs = std::pow(uc.CGSValue(density_dim),tw::Float(1-numBodies)) / uc.CGSValue(time_dim);
	std::string word;
	inputString >> word;
	if (word!="rate" && word!="janev_rate")
		throw tw::FatalError("Found "+word+" instead of <rate> or <janev_rate>.");
	if (word=="rate")
	{
		tw::input::PopExpectedWord(inputString,"=","reaction");
		if (!(inputString >> c1 >> c2 >> c3))
			throw tw::FatalError("Invalid number while reading rate coefficients for reaction.");
		c1 *= pow(unit_T_eV,c2);
		c1 /= unit_rate_cgs;
		c3 = uc.eV_to_sim(c3);
	}
	if (word=="janev_rate")
	{
		c1 = 0.0;
		tw::input::PopExpectedWord(inputString,"=","reaction");
		for (tw::Int i=0;i<9;i++)
			if (!(inputString >> b[i]))
				throw tw::FatalError("Invalid number while reading rate coefficients for reaction.");
		// janev coefficients are left in eV-cgs
	}
}

void Reaction::ReadInputFile(std::stringstream& inputString,tw::Float unitDensityCGS)
{
	UnitConverter uc(unitDensityCGS);
	std::string word;
	tw::Float sign,energy;
	bool rhs;
	numBodies = 0;

	tw::input::PopExpectedWord(inputString,"{","reaction");

	// Read in reactants and products
	do
	{
		sub.push_back(new SubReaction);
		sub.back()->heat = 0.0;
		sub.back()->vheat = 0.0;
		rhs = false;
		sign = 1.0;
		do
		{
			inputString >> word;
			if (word=="->")
				rhs = true;
			if (word=="+")
				sign = 1.0;
			if (word=="-")
				sign = -1.0;
			if (word!="->" && word!="+" && word!="-" && word!=":" && word!="}")
			{
				try
				{
					size_t pos;
					energy = uc.eV_to_sim(std::stod(word,&pos));
					if (pos<word.size()-1 || (pos==word.size()-1 && word.back()!='v'))
						throw tw::FatalError("Heat of reaction is ill-formed.");
					if (word.back()=='v')
						sub.back()->vheat = energy*sign;
					else
						sub.back()->heat = energy*sign;
				}
				catch (std::invalid_argument)
				{
					if (sign==-1.0)
						throw tw::FatalError("Chemical equation had misplaced minus sign.");
					if (rhs)
					{
						tw::input::StripQuotes(word);
						sub.back()->product_names.push_back(word);
					}
					else
					{
						tw::input::StripQuotes(word);
						sub.back()->reactant_names.push_back(word);
						numBodies++;
					}
				}
			}
		} while (word!=":" && word!="}");
	} while (word!="}");

	ReadRate(inputString,numBodies,uc);

	// Get temperature range and catalyst
	inputString >> catalyst_name >> word;
	tw::input::PythonRange(word,&T0,&T1);
	T0 = uc.eV_to_sim(T0);
	T1 = uc.eV_to_sim(T1);
}

void Excitation::ReadInputFile(std::stringstream& inputString,tw::Float unitDensityCGS)
{
	std::string word;
	UnitConverter uc(unitDensityCGS);

	inputString >> name1;
	tw::input::StripQuotes(name1);
	tw::input::PopExpectedWord(inputString,"->","excitation");
	inputString >> name2;
	tw::input::StripQuotes(name2);
	tw::input::PopExpectedWord(inputString,"level","excitation");
	tw::input::PopExpectedWord(inputString,"=","excitation");
	inputString >> level;

	T0 = 0.0;
	T1 = tw::big_pos;
	ReadRate(inputString,2,uc);
}

void Collision::ReadInputFile(std::stringstream& inputString,tw::Float unitDensityCGS)
{
	// hard sphere example
	// new collision = e <-> N , cross section = 1e-15

	// coulomb collision example
	// new collision = e <-> N[+] , coulomb

	// metallic collision example
	// new collision = e <-> Cu[+] , metallic , ks = 2.4 , fermi_energy_ev = 7.0 , ref_density = 1e23

	tw::Int i;
	std::string word,species;
	UnitConverter uc(unitDensityCGS);

	inputString >> name1;
	tw::input::StripQuotes(name1);
	tw::input::PopExpectedWord(inputString,"<->","collision");
	inputString >> name2;
	tw::input::StripQuotes(name2);

	inputString >> word;

	if (word!="cross" && word!="coulomb" && word!="metallic")
		throw tw::FatalError("Unrecognized collision type <"+word+">.");

	if (word=="cross")
	{
		type = sparc::hard_sphere;
		tw::input::PopExpectedWord(inputString,"section","hard sphere collision");
		tw::input::PopExpectedWord(inputString,"=","hard sphere collision");
		if (!(inputString >> crossSection))
			throw tw::FatalError("Invalid number encountered while reading collision cross section.");
		crossSection = uc.CGSToSim(cross_section_dim,crossSection);
	}

	if (word=="coulomb")
	{
		type = sparc::coulomb;
	}

	if (word=="metallic")
	{
		type = sparc::metallic;
		tw::input::PopExpectedWord(inputString,"ks","metallic collision");
		tw::input::PopExpectedWord(inputString,"=","metallic collision");
		if (!(inputString >> ks))
			throw tw::FatalError("Invalid number encountered while reading metallic collision.");
		tw::input::PopExpectedWord(inputString,"fermi_energy_ev","metallic collision");
		tw::input::PopExpectedWord(inputString,"=","metallic collision");
		if (!(inputString >> T_ref))
			throw tw::FatalError("Invalid number encountered while reading metallic collision.");
		tw::input::PopExpectedWord(inputString,"ref_density","metallic collision");
		tw::input::PopExpectedWord(inputString,"=","metallic collision");
		if (!(inputString >> n_ref))
			throw tw::FatalError("Invalid number encountered while reading metallic collision.");
		T_ref = uc.eV_to_sim(T_ref);
		n_ref = uc.CGSToSim(density_dim,n_ref);
	}
}
