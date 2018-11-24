#include "sim.h"

// Subreaction class serves as an element of the Reaction quasitool
// No need to write chemical names into restart file.  We have the indices into the state array.

void SubReaction::ReadData(std::ifstream& inFile)
{
	tw::Int i,num;
	sparc::hydro_set htemp;
	sparc::material mtemp;

	reactants.clear();
	mat_r.clear();
	inFile.read((char *)&num,sizeof(tw::Int));
	for (int i=0;i<num;i++)
	{
		inFile.read((char *)&htemp,sizeof(sparc::hydro_set));
		reactants.push_back(htemp);
		inFile.read((char *)&mtemp,sizeof(sparc::material));
		mat_r.push_back(mtemp);
	}

	products.clear();
	mat_p.clear();
	inFile.read((char *)&num,sizeof(tw::Int));
	for (int i=0;i<num;i++)
	{
		inFile.read((char *)&htemp,sizeof(sparc::hydro_set));
		products.push_back(htemp);
		inFile.read((char *)&mtemp,sizeof(sparc::material));
		mat_p.push_back(mtemp);
	}

	inFile.read((char *)&heat,sizeof(tw::Float));
	inFile.read((char *)&vheat,sizeof(tw::Float));
}

void SubReaction::WriteData(std::ofstream& outFile)
{
	tw::Int i,num;

	num = reactants.size();
	outFile.write((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		outFile.write((char*)&reactants[i],sizeof(sparc::hydro_set));
		outFile.write((char*)&mat_r[i],sizeof(sparc::material));
	}

	num = products.size();
	outFile.write((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		outFile.write((char*)&products[i],sizeof(sparc::hydro_set));
		outFile.write((char*)&mat_p[i],sizeof(sparc::material));
	}

	outFile.write((char *)&heat,sizeof(tw::Float));
	outFile.write((char *)&vheat,sizeof(tw::Float));
}

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

void PrimitiveReaction::ReadData(std::ifstream& inFile)
{
	inFile.read((char *)&unit_T_eV,sizeof(tw::Float));
	inFile.read((char *)&unit_rate_cgs,sizeof(tw::Float));
	inFile.read((char *)&c1,sizeof(tw::Float));
	inFile.read((char *)&c2,sizeof(tw::Float));
	inFile.read((char *)&c3,sizeof(tw::Float));
	inFile.read((char *)&T0,sizeof(tw::Float));
	inFile.read((char *)&T1,sizeof(tw::Float));
	inFile.read((char *)b,sizeof(tw::Float)*9);
}

void PrimitiveReaction::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&unit_T_eV,sizeof(tw::Float));
	outFile.write((char *)&unit_rate_cgs,sizeof(tw::Float));
	outFile.write((char *)&c1,sizeof(tw::Float));
	outFile.write((char *)&c2,sizeof(tw::Float));
	outFile.write((char *)&c3,sizeof(tw::Float));
	outFile.write((char *)&T0,sizeof(tw::Float));
	outFile.write((char *)&T1,sizeof(tw::Float));
	outFile.write((char *)b,sizeof(tw::Float)*9);
}

void Reaction::ReadInputFile(std::stringstream& inputString,tw::Float unitDensityCGS)
{
	std::string word;
	tw::Float sign,energy;
	bool rhs = false;

	numBodies = 0;

	inputString >> word >> word; // take off "{" and read next word

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

void Reaction::ReadData(std::ifstream& inFile)
{
	tw::Int i,num;
	PrimitiveReaction::ReadData(inFile);
	inFile.read((char *)&catalyst,sizeof(catalyst));
	inFile.read((char *)&numBodies,sizeof(numBodies));
	inFile.read((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
	{
		sub.push_back(new SubReaction);
		sub.back()->ReadData(inFile);
	}
}

void Reaction::WriteData(std::ofstream& outFile)
{
	tw::Int i,num;
	PrimitiveReaction::WriteData(outFile);
	outFile.write((char *)&catalyst,sizeof(catalyst));
	outFile.write((char *)&numBodies,sizeof(numBodies));
	num = sub.size();
	outFile.write((char *)&num,sizeof(tw::Int));
	for (i=0;i<num;i++)
		sub[i]->WriteData(outFile);
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

void Excitation::ReadData(std::ifstream& inFile)
{
	PrimitiveReaction::ReadData(inFile);
	inFile.read((char *)&h1,sizeof(h1));
	inFile.read((char *)&h2,sizeof(h2));
	inFile.read((char *)&e1,sizeof(e1));
	inFile.read((char *)&e2,sizeof(e2));
	inFile.read((char *)&m1,sizeof(m1));
	inFile.read((char *)&m2,sizeof(m2));
	inFile.read((char *)&level,sizeof(tw::Float));
}

void Excitation::WriteData(std::ofstream& outFile)
{
	PrimitiveReaction::WriteData(outFile);
	outFile.write((char *)&h1,sizeof(h1));
	outFile.write((char *)&h2,sizeof(h2));
	outFile.write((char *)&e1,sizeof(e1));
	outFile.write((char *)&e2,sizeof(e2));
	outFile.write((char *)&m1,sizeof(m1));
	outFile.write((char *)&m2,sizeof(m2));
	outFile.write((char *)&level,sizeof(tw::Float));
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

void Collision::ReadData(std::ifstream& inFile)
{
	inFile.read((char *)&type,sizeof(sparc::collisionType));
	inFile.read((char *)&h1,sizeof(h1));
	inFile.read((char *)&h2,sizeof(h2));
	inFile.read((char *)&e1,sizeof(e1));
	inFile.read((char *)&e2,sizeof(e2));
	inFile.read((char *)&m1,sizeof(m1));
	inFile.read((char *)&m2,sizeof(m2));
	inFile.read((char *)&crossSection,sizeof(tw::Float));
	inFile.read((char *)&ks,sizeof(tw::Float));
	inFile.read((char *)&T_ref,sizeof(tw::Float));
	inFile.read((char *)&n_ref,sizeof(tw::Float));
}

void Collision::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&type,sizeof(sparc::collisionType));
	outFile.write((char *)&h1,sizeof(h1));
	outFile.write((char *)&h2,sizeof(h2));
	outFile.write((char *)&e1,sizeof(e1));
	outFile.write((char *)&e2,sizeof(e2));
	outFile.write((char *)&m1,sizeof(m1));
	outFile.write((char *)&m2,sizeof(m2));
	outFile.write((char *)&crossSection,sizeof(tw::Float));
	outFile.write((char *)&ks,sizeof(tw::Float));
	outFile.write((char *)&T_ref,sizeof(tw::Float));
	outFile.write((char *)&n_ref,sizeof(tw::Float));
}
