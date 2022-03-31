namespace sparc
{
	enum collisionType { hard_sphere, coulomb, metallic };
}

// Quasitools for chemical reactions and collisions

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
	void ReadRate(std::stringstream& inputString,tw::Int numBodies,const tw::UnitConverter& uc);
};

struct Reaction : PrimitiveReaction
{
	std::vector<SubReaction*> sub;
	std::string catalyst_name;
	sparc::eos_set catalyst;
	tw::Int numBodies;

	~Reaction();
	virtual void ReadInputFile(std::stringstream& inputString,const tw::UnitConverter& uc);
};

struct Excitation : PrimitiveReaction
{
	std::string name1,name2;// 1 excites 2
	sparc::hydro_set h1,h2;
	sparc::eos_set e1,e2;
	sparc::material m1,m2;
	tw::Float level;

	virtual void ReadInputFile(std::stringstream& inputString,const tw::UnitConverter& uc);
};

struct Collision
{
	std::string name1,name2;
	sparc::hydro_set h1,h2;
	sparc::eos_set e1,e2;
	sparc::material m1,m2;
	sparc::collisionType type;
	tw::Float crossSection;
	tw::Float ks,T_ref,n_ref;

	virtual void ReadInputFile(std::stringstream& inputString,const tw::UnitConverter& uc);
};
