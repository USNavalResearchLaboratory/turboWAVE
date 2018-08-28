enum tw_ionization_model {noIonization,ADKTunneling,pptIonization,mpiSimple};

enum tw_dimensions
{
	time_dim,
	length_dim,
	number_dim,
	mass_dim,
	energy_dim,
	momentum_dim,
	angular_momentum_dim,
	density_dim,
	power_dim,
	fluence_dim,
	intensity_dim,
	energy_density_dim,
	power_density_dim,
	charge_dim,
	current_dim,
	current_density_dim,
	charge_density_dim,
	electric_field_dim,
	magnetic_field_dim,
	scalar_potential_dim,
	conductivity_dim,
	rate_coefficient_2_dim,
	rate_coefficient_3_dim,
	mobility_dim,
	temperature_dim,
	cross_section_dim,
	susceptibility_dim,
	susceptibility_2_dim,
	susceptibility_3_dim
};

struct UnitConverter
{
	tw::Float n1,wp;
	tw::Float c,qe,me,eps0,kB,hbar;

	UnitConverter(tw::Float unitDensityCGS);
	tw::Float MKSValue(tw_dimensions dim) const;
	tw::Float CGSValue(tw_dimensions dim) const;

	tw::Float SimToMKS(tw_dimensions dim,tw::Float val) const
	{
		return val*MKSValue(dim);
	}
	tw::Float MKSToSim(tw_dimensions dim,tw::Float val) const
	{
		return val/MKSValue(dim);
	}
	tw::Float SimToCGS(tw_dimensions dim,tw::Float val) const
	{
		return val*CGSValue(dim);
	}
	tw::Float CGSToSim(tw_dimensions dim,tw::Float val) const
	{
		return val/CGSValue(dim);
	}
	tw::Float sim_to_eV(tw::Float val) const
	{
		return val*MKSValue(energy_dim)/qe;
	}
	tw::Float eV_to_sim(tw::Float val) const
	{
		return val*qe/MKSValue(energy_dim);
	}
};

struct IonizationData
{
	tw::Float ionizationPotential;
	tw::Float electrons,protons;

	tw_ionization_model ionizationModel;
	tw::Float adkMultiplier,pptMultiplier;
	tw::Float C_PPT,C_ADK,C_ADK_AVG,nstar,lstar;
	tw::Float t_atomic,E_atomic,f_atomic_to_sim,E_sim_to_atomic;
	tw::Float photons,w0,E_MPI,max_rate;
	tw::Int terms;

	bool ionizeFromGas;
	tw::Int ionSpecies,electronSpecies;

	IonizationData();
	void Initialize(tw::Float unitDensity,tw::Float* carrierFrequency);
	tw::Float wfunc(tw::Float x);
	tw::Float Rate(tw::Float instant,tw::Float peak);
	void ReadInputFileTerm(std::stringstream& inputString,std::string& command);
};

// current eos variables (for reference)
// -------------------------------------
// t    : temperature
// tv   : vibrational temperature
// p    : pressure
// k    : thermometricConductivity * nmcv
// visc : kinematicViscosity * nm

// ASHER_MOD
// This is the master class of all EOS computation related objects
struct EOSDataTool:ComputeTool
{
	EOSDataTool(MetricSpace *m,Task *tsk,bool shared);

	// The following are modifications of functions copied over from the EquilibriumGroup class
	// definition in 'fluid.h'. They seem to be an integral part of the calculation of
	// various EOS-related quantities, and so may be advantageous to have them included inside
	// the EOSDataTool object. Some arguments were added because EOSDataTool cannot directly
	// access the EquilibriumGroup's internal variables.

	// std::valarray<tw::Float> mass,charge,cvm,excitationEnergy,thermo_cond_cvm,k_visc_m;

	tw::Float DensitySum(Element e,const Field& f,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=e.low;s<=e.high;s++)
			ans += f(cell,s);
		return ans;
	}
	tw::Float DensityWeightedSum(Element e,const Field& f,std::valarray<tw::Float>& qty,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=e.low;s<=e.high;s++)
			ans += f(cell,s)*qty[s-e.low];
		return ans;
	}
	tw::Float ConditionalDensitySum(Element e,const Field& f,std::valarray<tw::Float>& qty,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=e.low;s<=e.high;s++)
			ans += qty[s-e.low] > 0.0 ? f(cell,s) : 0.0;
		return ans;
	}

	tw::vec3 Velocity(Element e,tw::Int npx,tw::Int npy,tw::Int npz, const Field& f,std::valarray<tw::Float> mass, const CellIterator& cell)
	{
		tw::Float nm = DensityWeightedSum(e,f,mass,cell);
		tw::vec3 chemVelocity = tw::vec3(f(cell,npx),f(cell,npy),f(cell,npz));
		return tw::vec3(f(cell,npx),f(cell,npy),f(cell,npz))/(tw::small_pos + nm);
	}

	// The following functions have many arguments because EOSDataTool cannot directly access
	// the many EOS related internal variables of EquilibriumGroup, and therefore must be
	// passed on as function arguments. In the long term we hope to migrate these terms into
	// the EOSDataTool object, and if they need to referenced in the EquilibriumGroup level,
	// it may be done so with polymorphism
	// (for example 'group->eosData->T' instead of 'group->T' for the temperature index)
	// The result of which would be that the following functions will be vastly simplified.
	virtual void ApplyEOS(
								std::valarray<tw::Float> mass,std::valarray<tw::Float> charge, std::valarray<tw::Float> cvm,
								std::valarray<tw::Float> excitationEnergy, std::valarray<tw::Float> thermo_cond_cvm,
								std::valarray<tw::Float> k_visc_m, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos) {;}

	// mass, charge, and cvm are simply float values for a single Chemical object (including electrons)
    // notice two new arguments, 'ie' and 'nu_e', which are actually elements of Chemistry
	virtual void ApplyEOS(
								tw::Float mass,tw::Float charge,tw::Float cvm,
								tw::Float excitationEnergy,tw::Float thermometricConductivity,
								tw::Float kinematicViscosity, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos) {;}

	static tw_tool ReadInputFileTerm_GetToolType(std::stringstream& inputString,std::string& command);
	virtual void ReadInputFileBlock(std::stringstream& inputString) {;}
};

// EOS Calculation for a single Chemical object with ideal gas properties,
//   embedded in an EqiulibriumGroup object that uses EOSMixture
//
// The following is different from the EOSIdealGasMixture in that the calculations are
// made at the eosData (or Chemical) level.
struct EOSIdealGas:EOSDataTool
{
	EOSIdealGas(MetricSpace *m,Task *tsk,bool shared);

	// mass, charge, and cvm are simply float values for electrons
    // notice two new arguments, 'ie' and 'nu_e', which are actually elements of Chemistry
	virtual void ApplyEOS(
								tw::Float mass,tw::Float charge,tw::Float cvm,
								tw::Float excitationEnergy,tw::Float thermometricConductivity,
								tw::Float kinematicViscosity, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos);
};

// EOS calculations for the electrons, treating them like a type of ideal gas
//
// This contains the same calculations that were originally embedded in Turbowave
// Since they are EOS calculations of a particular Chemical type, I place it here
// as its own EOSDataTool object. Currently this is the only electron EOS model
// there is, but in the future we may include other electron models
// (like Thomas-Fermi, etc.)
struct EOSHotElectrons:EOSDataTool
{
	EOSHotElectrons(MetricSpace *m,Task *tsk,bool shared);

	virtual void ApplyEOS(
								tw::Float mass,tw::Float charge,tw::Float cvm,
								tw::Float excitationEnergy,tw::Float thermometricConductivity,
								tw::Float kinematicViscosity, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos);
};

// EOS calculations for a EquilibriumGroup object, which may contain one or multiple Chemical objects
//
// This contains a vector of EOSDataTools, each EOSData corresponding to the EOS of an EqiulibriumGroup
// component chemical. The EOS quantities that are singular to the group (like T, Tv) are calculated
// at the EOSMixture level, and the EOS quantities that have different model contributions for each
// Chemical (like P, visc), are calculated at the eosData level and added/weighted at the
// EOSMixture level accordingly.
struct EOSMixture:EOSDataTool
{
	std::vector<EOSDataTool*> eosComponents;
	Field eos_tmp, eos_tmp2;

	EOSMixture(MetricSpace *m,Task *tsk,bool shared);

	virtual void ApplyEOS(
								std::valarray<tw::Float> mass,std::valarray<tw::Float> charge, std::valarray<tw::Float> cvm,
								std::valarray<tw::Float> excitationEnergy, std::valarray<tw::Float> thermo_cond_cvm,
								std::valarray<tw::Float> k_visc_m, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos);
};

// EOS for a mix of gases, assuming all components are ideal gases
//
// this is the code that is identical to the original EOS calculation in turbowave
// you can use this to compare with results that would be produced with the original calculations
// you would merely have to initialize this datatool instead of eosMixture in EquilibriumGroup
// all component chemicals are treated as ideal gases, regardless of the eosDataTool they contain
struct EOSIdealGasMix:EOSMixture
{
	std::vector<EOSDataTool*> eosComponents;

	EOSIdealGasMix(MetricSpace *m,Task *tsk,bool shared);

	virtual void ApplyEOS(
								std::valarray<tw::Float> mass,std::valarray<tw::Float> charge, std::valarray<tw::Float> cvm,
								std::valarray<tw::Float> excitationEnergy, std::valarray<tw::Float> thermo_cond_cvm,
								std::valarray<tw::Float> k_visc_m, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos);
};

// This is the most basic implementation of the MieGruneisen EOS
//
// It assumes that GRUN = GRUN0 = const. at all times
// This is typically not used in literature concerning MieGruneisen EOSs,
// as the results are rarely physicaly accurate
// If you produce sound waves with this model (for example with Cu),
// you'll notice the sound speed is off. Qualitatively, it gives a broad picture.
struct EOSMieGruneisen:EOSDataTool
{
	tw::Float GRUN; // Gruneisen coefficient

	EOSMieGruneisen(MetricSpace *m,Task *tsk,bool shared);

	virtual void ApplyEOS(
								tw::Float mass,tw::Float charge,tw::Float cvm,
								tw::Float excitationEnergy,tw::Float thermometricConductivity,
								tw::Float kinematicViscosity, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos);

	virtual void ReadInputFileBlock(std::stringstream& inputString);
};

// This is a MieGruneisen EOS that assumes \rho * GRUN = const., and a linear Hugoniot fit
//
// This is what is more typically what is found in literature regarding MieGruneisen models
// You'll get a more accurate sound speed and shock speeds, assuming that the simulation is
// within the range of relevant Hugoniot data and model assumptions are properly met
struct EOSMieGruneisen2:EOSDataTool
{
	tw::Float GRUN; // Gruneisen coefficient
	tw::Float n0;   // Reference density
	tw::Float c0;   // y - intercept of Hugoniot fit (usually appriximately speed of sound)
	tw::Float S1;   // coefficient of linear fit of Hugoniot data

	EOSMieGruneisen2(MetricSpace *m,Task *tsk,bool shared);

	virtual void ApplyEOS(
								tw::Float mass,tw::Float charge,tw::Float cvm,
								tw::Float excitationEnergy,tw::Float thermometricConductivity,
								tw::Float kinematicViscosity, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos);

	virtual void ReadInputFileBlock(std::stringstream& inputString);
};
