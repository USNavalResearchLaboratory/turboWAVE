module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"

export module chemistry;
import input;
import compute_tool;
import physics;
import logger;

export namespace sparc
{
	enum collisionType { hard_sphere, coulomb, metallic };
}

// Quasitools for chemical reactions and collisions

tw::Float ReadChemList(TSTreeCursor *curs,const std::string& src,std::vector<std::string>& names);

struct SubReaction
{
	std::vector<std::string> reactant_names,product_names;
	std::vector<sparc::hydro_set> reactants,products;
	std::vector<sparc::material> mat_r,mat_p;
	tw::Float heat,vheat;
};

struct PrimitiveReaction
{
	tw::Float T0,T1; // temperature range
	tw::Float c1,c2,c3; // arrhenius form
	tw::Float b[9]; // janev coefficients
	tw::Float unit_T_eV,unit_rate_cgs; // normalization help for janev

	tw::Float PrimitiveRate(tw::Float T);
	void ReadRate(TSTreeCursor *curs,const std::string& src,tw::Int numBodies,const tw::UnitConverter& uc);
};

export struct Reaction : PrimitiveReaction
{
	std::vector<SubReaction*> sub;
	std::string catalyst_name;
	sparc::eos_set catalyst;
	tw::Int numBodies;

	virtual ~Reaction();
	virtual void ReadInputFile(TSTreeCursor *curs,const std::string& src,const tw::UnitConverter& uc);
};

export struct Excitation : PrimitiveReaction
{
	std::string name1,name2;// 1 excites 2
	sparc::hydro_set h1,h2;
	sparc::eos_set e1,e2;
	sparc::material m1,m2;
	tw::Float level;

	virtual ~Excitation() {}
	virtual void ReadInputFile(TSTreeCursor *curs,const std::string& src,const tw::UnitConverter& uc);
};

export struct Collision
{
	std::string name1,name2;
	sparc::hydro_set h1,h2;
	sparc::eos_set e1,e2;
	sparc::material m1,m2;
	sparc::collisionType type;
	tw::Float crossSection;
	tw::Float ks,T_ref,n_ref;

	virtual ~Collision() {}
	virtual void ReadInputFile(TSTreeCursor *curs,const std::string& src,const tw::UnitConverter& uc);
};

tw::Float PrimitiveReaction::PrimitiveRate(tw::Float T)
{
	tw::Float rate;

	if (T<T0 || T>T1)
		return 0.0;

	if (c1==0.0) // use janev form
	{
		const tw::Float logT_eV = std::log(tw::small_pos + std::fabs(T)*unit_T_eV);
		rate = 0.0;
		for (tw::Int s=0;s<9;s++)
			rate += b[s]*pow(logT_eV,tw::Float(s));
		rate = std::exp(rate) / unit_rate_cgs;
	}
	else // use arrhenius form
		rate = c1 * pow(T,c2) * std::exp(-c3/T);

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
		tw::input::ThrowParsingError(curs,src,"expected rate node");
	ts_tree_cursor_goto_first_child(curs); // go in rate
	if (tw::input::node_kind(curs) == "arrhenius") {
		ts_tree_cursor_goto_first_child(curs); // go in arrhenius
		c1 = stod(tw::input::next_named_node_text(curs,src),NULL);
		c2 = stod(tw::input::next_named_node_text(curs,src),NULL);
		c3 = stod(tw::input::next_named_node_text(curs,src),NULL);
		c1 *= pow(unit_T_eV,c2);
		c1 /= unit_rate_cgs;
		c3 = c3*tw::dims::temperature >> cgs >> native;
		ts_tree_cursor_goto_parent(curs);
	} else if (tw::input::node_kind(curs) == "janev") {
		ts_tree_cursor_goto_first_child(curs); // go in janev
		c1 = 0.0;
		for (tw::Int i=0;i<9;i++) {
			b[i] = stod(tw::input::next_named_node_text(curs,src),NULL);
		}
		ts_tree_cursor_goto_parent(curs);
		// janev coefficients are left in eV-cgs
	} else {
		tw::input::ThrowParsingError(curs,src,"unexpected node in reaction rate");
	}
	ts_tree_cursor_goto_parent(curs);
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
	numBodies = 0;
	tw::UnitConverter cgs(tw::units::cgs,native);

	ts_tree_cursor_goto_first_child(curs);
	tw::input::next_named_node(curs,false);
	
	logger::TRACE("parse the formula");
	ts_tree_cursor_goto_first_child(curs); // go in full_formula node
	while (tw::input::next_named_node(curs,false))
	{
		logger::TRACE(std::format("parse sub-formula {}",tw::input::node_text(curs,src)));
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

	logger::TRACE("parse the rate");
	tw::input::next_named_node(curs,false);
	ReadRate(curs,src,numBodies,native);

	logger::TRACE("parse the catalyst");
	tw::input::next_named_node(curs,false);
	catalyst_name = tw::input::PythonRange(curs,src,&T0,&T1);
	T0 = T0*tw::dims::temperature >> cgs >> native;
	T1 = T1*tw::dims::temperature >> cgs >> native;
}

void Excitation::ReadInputFile(TSTreeCursor *curs,const std::string& src,const tw::UnitConverter& native)
{
	ts_tree_cursor_goto_first_child(curs);
	name1 = tw::input::next_named_node_text(curs,src);
	name2 = tw::input::next_named_node_text(curs,src);
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

	std::string word,species;
	tw::UnitConverter cgs(tw::units::cgs,native);

	ts_tree_cursor_goto_next_sibling(curs);
	if (tw::input::node_text(curs,src,true) == "coulomb") {
		type = sparc::coulomb;
	} else if (tw::input::node_text(curs,src,true) == "metallic") {
		type = sparc::metallic;
		ks = std::stod(tw::input::next_named_node_text(curs,src));
		T_ref = std::stod(tw::input::next_named_node_text(curs,src));
		n_ref = std::stod(tw::input::next_named_node_text(curs,src));
		T_ref = T_ref*tw::dims::temperature >> cgs >> native;
		n_ref = n_ref*tw::dims::density >> cgs >> native;
	} else if (tw::input::node_text(curs,src,true) == "cross") {
		type = sparc::hard_sphere;
		crossSection = std::stod(tw::input::next_named_node_text(curs,src));
		crossSection = crossSection * tw::dims::cross_section >> cgs >> native;
	} else {
		tw::input::ThrowParsingError(curs, src, std::format("unexpected token `{}`",tw::input::node_text(curs,src)));
	}
}
