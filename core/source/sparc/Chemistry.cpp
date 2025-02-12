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

/// @brief read either arrhenius or janev rate
/// @param curs starting on rate node
/// @param src source document
/// @param numBodies number of reactants
/// @param native units
void PrimitiveReaction::ReadRate(TSTreeCursor *curs,const std::string& src,tw::Int numBodies,const tw::UnitConverter& native)
{
	tw::UnitConverter cgs(tw::units::cgs,native);
	unit_T_eV = 1.0*tw::dims::temperature >> native >> cgs;
	unit_rate_cgs = std::pow(1.0*tw::dims::density>>native>>cgs,tw::Float(1-numBodies)) / (1.0*tw::dims::time>>native>>cgs);
	std::string word;

	if (tw::input::node_kind(curs) != "rate")
		throw tw::FatalError("expected rate node");
	ts_tree_cursor_goto_first_child(curs);
	if (tw::input::node_kind(curs) == "arrhenius") {
		c1 = stod(tw::input::next_named_node_text(curs,src),NULL);
		c2 = stod(tw::input::next_named_node_text(curs,src),NULL);
		c3 = stod(tw::input::next_named_node_text(curs,src),NULL);
		c1 *= pow(unit_T_eV,c2);
		c1 /= unit_rate_cgs;
		c3 = c3*tw::dims::temperature >> cgs >> native;
	} else if (tw::input::node_kind(curs) == "janev") {
		c1 = 0.0;
		for (tw::Int i=0;i<9;i++) {
			b[i] = stod(tw::input::next_named_node_text(curs,src),NULL);
		}
		// janev coefficients are left in eV-cgs
	} else {
		throw tw::FatalError("unexpected node in reaction rate");
	}
}

/// @brief read half of subreaction
/// @param curs on chems node
/// @param src source document
/// @param names vector to hold the names
/// @return the heat of reaction
tw::Float ReadChemList(TSTreeCursor *curs,const std::string& src,std::vector<std::string>& names) {
	tw::Float sign = 1.0;
	tw::Float heat = 0.0;
	ts_tree_cursor_goto_first_child(curs); // go in chems
	do {
		if (tw::input::trim(tw::input::node_text(curs,src)) == "+")
			sign = 1.0;
		else if (tw::input::trim(tw::input::node_text(curs,src)) == "-")
			sign = -1.0;
		else if (tw::input::node_kind(curs) == "identifier") {
			std::string word = tw::input::node_text(curs,src);
			tw::input::StripQuotes(word);
			names.push_back(word);
		} else if (tw::input::node_kind(curs) == "decimal") {
			heat = sign*std::stod(tw::input::node_text(curs,src));
		}
	} while (ts_tree_cursor_goto_next_sibling(curs));
	ts_tree_cursor_goto_parent(curs);
	return heat;
}


Reaction::~Reaction()
{
	for (auto s : sub)
		delete s;
}

void Reaction::ReadInputFile(TSTreeCursor *curs,const std::string& src,const tw::UnitConverter& native)
{
	std::string word;
	tw::Float sign,energy;
	bool rhs;
	numBodies = 0;
	tw::UnitConverter cgs(tw::units::cgs,native);

	ts_tree_cursor_goto_first_child(curs);
	tw::input::next_named_node(curs,false);
	ts_tree_cursor_goto_first_child(curs); // go in full_formula node
	// Read in reactants and products
	while (tw::input::next_named_node(curs,false))
	{
		ts_tree_cursor_goto_first_child(curs); // go in subformula
		sub.push_back(new SubReaction);
		sub.back()->heat = ReadChemList(curs,src,sub.back()->reactant_names);
		numBodies += sub.back()->reactant_names.size();
		tw::input::next_named_node(curs,false);
		sub.back()->heat += ReadChemList(curs,src,sub.back()->product_names);
		sub.back()->vheat = 0.0;
		ts_tree_cursor_goto_parent(curs);
	}

	ts_tree_cursor_goto_parent(curs);
	tw::input::next_named_node(curs,false);
	ReadRate(curs,src,numBodies,native);

	// Get temperature range and catalyst
	catalyst_name = tw::input::PythonRange(curs,src,&T0,&T1);
	T0 = T0*tw::dims::temperature >> cgs >> native;
	T1 = T1*tw::dims::temperature >> cgs >> native;
}

void Excitation::ReadInputFile(TSTreeCursor *curs,const std::string& src,const tw::UnitConverter& native)
{
	ts_tree_cursor_goto_first_child(curs);
	name1 = tw::input::next_named_node_text(curs,src);
	name2 = tw::input::next_named_node_text(curs,src);
	tw::input::StripQuotes(name1);
	tw::input::StripQuotes(name2);
	level = std::stod(tw::input::next_named_node_text(curs,src));

	T0 = 0.0;
	T1 = tw::big_pos;
	ts_tree_cursor_goto_next_sibling(curs);
	ReadRate(curs,src,2,native);
}

void Collision::ReadInputFile(TSTreeCursor *curs,const std::string& src,const tw::UnitConverter& native)
{
	// hard sphere example
	// new collision = e <-> N , cross section = 1e-15

	// coulomb collision example
	// new collision = e <-> N[+] , coulomb

	// metallic collision example
	// new collision = e <-> Cu[+] , metallic , ks = 2.4 , fermi_energy_ev = 7.0 , ref_density = 1e23

	ts_tree_cursor_goto_first_child(curs);
	name1 = tw::input::next_named_node_text(curs,src);
	name2 = tw::input::next_named_node_text(curs,src);
	tw::input::StripQuotes(name1);
	tw::input::StripQuotes(name2);

	std::string word,species;
	tw::UnitConverter cgs(tw::units::cgs,native);

	ts_tree_cursor_goto_next_sibling(curs);
	if (tw::input::node_text(curs,src) == "coulomb") {
		type = sparc::coulomb;
	} else if (tw::input::node_text(curs,src) == "metallic") {
		type = sparc::metallic;
		ks = std::stod(tw::input::next_named_node_text(curs,src));
		T_ref = std::stod(tw::input::next_named_node_text(curs,src));
		n_ref = std::stod(tw::input::next_named_node_text(curs,src));
		T_ref = T_ref*tw::dims::temperature >> cgs >> native;
		n_ref = n_ref*tw::dims::density >> cgs >> native;
	} else if (tw::input::node_text(curs,src) == "cross") {
		type = sparc::hard_sphere;
		crossSection = std::stod(tw::input::next_named_node_text(curs,src));
		crossSection = crossSection * tw::dims::cross_section >> cgs >> native;
	} else {
		throw tw::FatalError("unexpected token in collision");
	}
}
