#include "meta_base.h"
#include "computeTool.h"
#include "physics.h"
#include "chemistry.h"

// Subreaction class serves as an element of the Reaction quasitool
// No need to write chemical names into restart file.  We have the indices into the state array.

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

void Reaction::ReadInputFile(std::stringstream& inputString,tw::Float unitDensityCGS)
{
	std::string word;
	tw::Float sign,energy;
	bool rhs = false;

	numBodies = 0;

	inputString >> word >> word; // take off "{" and read next word (n.b. "=" terminates the preamble)

	UnitConverter uc(unitDensityCGS);

	// Read in reactants and products
	do
	{
		sub.push_back(new SubReaction);
		sub.back()->heat = 0.0;
		sub.back()->vheat = 0.0;
		rhs = false;
		do
		{
			if (word=="->")
			{
				rhs = true;
				inputString >> word;
			}
			if (word=="+")
			{
				sign = 1.0;
				inputString >> word;
			}
			if (word=="-")
			{
				sign = -1.0;
				inputString >> word;
			}
			if (word==":")
				inputString >> word;

			// see if this is a heat of reaction
			// if no number is found, encode this result by setting energy=0
			// WARNING: this fails if a species name begins with inf or nan
			try { energy = uc.eV_to_sim(std::stod(word,NULL)); }
			catch (std::invalid_argument) { energy = 0.0; }
			if (energy!=0.0) // it is
			{
				if (word.back()=='v')
					sub.back()->vheat = energy*sign;
				else
					sub.back()->heat = energy*sign;
			}
			else
			{
				if (rhs)
				{
					sub.back()->product_names.push_back(word);
				}
				else
				{
					sub.back()->reactant_names.push_back(word);
					numBodies++;
				}
			}
			inputString >> word;
		} while (word!=":" && word!="}");
	} while (word!="}");

	// Get reaction rate data and normalize
	unit_T_eV = uc.sim_to_eV(1.0);
	unit_rate_cgs = std::pow(uc.CGSValue(density_dim),tw::Float(1-numBodies)) / uc.CGSValue(time_dim);
	inputString >> word;
	if (word=="rate")
	{
		inputString >> word >> c1 >> c2 >> c3;
		c1 *= pow(uc.eV_to_sim(1.0),-c2);
		c1 /= unit_rate_cgs;
		c3 = uc.eV_to_sim(c3);
	}
	if (word=="janev_rate")
	{
		c1 = 0.0;
		inputString >> word;
		for (tw::Int i=0;i<9;i++)
			inputString >> b[i];
		// janev coefficients are left in eV-cgs
	}

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
	inputString >> word >> name2; // take off "->" and get name2
	inputString >> word >> word >> level;

	T0 = 0.0;
	T1 = tw::big_pos;
	unit_T_eV = uc.sim_to_eV(1.0);
	unit_rate_cgs = uc.SimToCGS(rate_coefficient_2_dim,1.0);
	inputString >> word;
	if (word=="rate")
	{
		inputString >> word >> c1 >> c2 >> c3;
		c1 *= pow(uc.eV_to_sim(1.0),-c2);
		c1 /= unit_rate_cgs;
		c3 = uc.eV_to_sim(c3);
	}
	if (word=="janev_rate")
	{
		c1 = 0.0;
		inputString >> word;
		for (tw::Int i=0;i<9;i++)
			inputString >> b[i];
		// janev coefficients are left in eV-cgs
	}
}

void Collision::ReadInputFile(std::stringstream& inputString,tw::Float unitDensityCGS)
{
	// hard sphere example
	// new collision = e <-> N , cross section = 5.0

	// coulomb collision example
	// new collision = e <-> N[+] , coulomb

	// metallic collision example
	// new collision = e <-> Cu[+] , metallic , ks = 2.4 , fermi_energy_ev = 7.0 , ref_density = 3000

	tw::Int i;
	std::string word,species;
	UnitConverter uc(unitDensityCGS);

	inputString >> name1;
	inputString >> word >> name2; // take off "->" and get name2

	inputString >> word;

	if (word=="cross")
	{
		type = sparc::hard_sphere;
		inputString >> word >> word >> crossSection;
	}

	if (word=="coulomb")
	{
		type = sparc::coulomb;
	}

	if (word=="metallic")
	{
		type = sparc::metallic;
		inputString >> word >> word >> ks;
		inputString >> word >> word >> T_ref;
		inputString >> word >> word >> n_ref;
		T_ref *= uc.eV_to_sim(1.0);
	}
}
