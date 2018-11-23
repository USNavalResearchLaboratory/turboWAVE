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

namespace sparc
{
	struct hydro_set
	{
		// num = number of constituent chemicals
		// first = index of first chemical density (assuming contiguous sequence)
		// last = index of last hydro quantity in the set (assuming contiguous sequence)
		// ni = index to a particular chemical density (optionally used)
		// npx,npy,npz = index to momentum density
		// u = index to energy density
		// x = index to vibrational density
		tw::Int num,first,last,ni,npx,npy,npz,u,x;
		tw::Int Load(tw::Int i,tw::Int N)
		{
			num=N;
			first=i;
			npx=i+N;
			npy=i+N+1;
			npz=i+N+2;
			u=i+N+3;
			x=i+N+4;
			last=x;
			return last+1;
		}
		static const tw::Int count = 5; // gives count of state components not counting densities
	};

	struct eos_set
	{
		// T = index to translational and rotational temperature
		// Tv = vibrational temperature (if used)
		// P = pressure
		// K = heat conductivity (check)
		// visc = dynamic viscosity (check)
		// Cv = heat capacity at constant volume
		tw::Int T,Tv,P,K,visc,Cv;
		tw::Int Load(tw::Int i) { T=i;Tv=i+1;P=i+2;K=i+3;visc=i+4;Cv=i+5;return i+6; }
		static const tw::Int count = 6;
	};

	struct material
	{
		// Characteristics of a chemical
		tw::Float mass,charge,cvm,excitationEnergy,thermometricConductivity,kinematicViscosity,eps_r,eps_i;
		void ReadInputFileDirective(std::stringstream& inputString,const std::string& command,const UnitConverter& uc);
	};

	struct material_set
	{
		// package multiple materials in a way that could promote optimization
		std::valarray<tw::Float> mass,charge,cvm,excitationEnergy,thermo_cond_cvm,k_visc_m,eps_r,eps_i;
		void Allocate(tw::Int num)
		{
			mass.resize(num);
			charge.resize(num);
			cvm.resize(num);
			excitationEnergy.resize(num);
			thermo_cond_cvm.resize(num);
			k_visc_m.resize(num);
			eps_r.resize(num);
			eps_i.resize(num);
		}
		void AddMaterial(const material& mat,const tw::Int& i)
		{
			mass[i] = mat.mass;
			charge[i] = mat.charge;
			cvm[i] = mat.cvm;
			excitationEnergy[i] = mat.excitationEnergy;
			thermo_cond_cvm[i] = mat.thermometricConductivity * mat.cvm;
			k_visc_m[i] = mat.kinematicViscosity * mat.mass;
			eps_r[i] = mat.eps_r;
			eps_i[i] = mat.eps_i;
		}
	};

	tw::Float CoulombCrossSection(const UnitConverter& uc,tw::Float q1,tw::Float q2,tw::Float m12,tw::Float v12,tw::Float N1,tw::Float N2,tw::Float T1,tw::Float T2);
	tw::Float ElectronPhononRateCoeff(const UnitConverter& uc,tw::Float Ti,tw::Float EFermi,tw::Float ks,tw::Float nref);
}

namespace tw
{
	enum class ionization_model {none,ADK,PPT,MPI};
}

struct IonizationData
{
	// this quasitool requires the owner to manage indexing of ionized species
	tw::Float ionizationPotential;
	tw::Float electrons,protons;

	tw::ionization_model ionizationModel;
	tw::Float adkMultiplier,pptMultiplier;
	tw::Float C_PPT,C_ADK,C_ADK_AVG,nstar,lstar;
	tw::Float t_atomic,E_atomic,f_atomic_to_sim,E_sim_to_atomic;
	tw::Float photons,w0,E_MPI,max_rate;
	tw::Int terms;

	// members determining species involved
	bool ionizeFromGas;
	std::string ion_name,electron_name;
	tw::Int ionSpecies,electronSpecies; // for PIC use the module index
	sparc::hydro_set hi,he,hgas; // for hydro use the field indexing

	IonizationData();
	void Initialize(tw::Float unitDensity,tw::Float* carrierFrequency);
	tw::Float wfunc(tw::Float x);
	tw::Float Rate(tw::Float instant,tw::Float peak);
	void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	void ReadData(std::ifstream& inFile);
	void WriteData(std::ofstream& outFile);
};

// This is the master class of all EOS computation related objects.
// This now retains its own indexing and material information in 3 lightweight structs:
// sparc::hydro_set and sparc::eos_set contain indices into a Field object
// sparc::material contains floats characterizing a chemical
struct EOSDataTool:ComputeTool
{
	sparc::hydro_set hidx;
	sparc::eos_set eidx;
	sparc::material mat;

	EOSDataTool(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize(const sparc::hydro_set&,const sparc::eos_set&,const sparc::material&);

	// The velocity was only used to get the internal energy, so make a function to do that instead
	// Move summation functions to the mixture class
	tw::Float InternalEnergy(const tw::Float& nm,const Field& f,const CellIterator& cell)
	{
		const tw::vec3 np = tw::vec3(f(cell,hidx.npx),f(cell,hidx.npy),f(cell,hidx.npz));
		const tw::vec3 vel = np/(tw::small_pos + nm);
		const tw::Float KE = 0.5*nm*Norm(vel);
		const tw::Float primitive = f(cell,hidx.u) - KE;
		const tw::Float failsafe = 1e-6*KE + tw::small_pos;
		const tw::Float sel = tw::Float(primitive>0.0);
		return sel*primitive + (1.0-sel)*failsafe;
	}

	// Used to have two versions of ApplyEOS for vector and scalar materials parameters.
	// Here there is a different function for component contributions (see other comments below)
	virtual void ComponentContribution(Field& hydro,Field& eos) {;}
	// Note ie is not used anywhere, but leave it for the present.
	// Conceptulaly, nu_e is something that should be computed by EOS tools rather than provided externally, but how?
	virtual void ApplyEOS(tw::Int ie, ScalarField& nu_e,Field& hydro, Field& eos) {;}
};

// EOS Calculation for a single Chemical object with ideal gas properties,
//   embedded in an EquilibriumGroup object that uses EOSMixture
// The distinction between mixture and non-mixture classes is maybe not needed anymore?
// Instead we are distinguishing two functions: ComponentContribution puts in something additive from each chemical (e.g., pressure)
// The ApplyEOS figures out the whole EOS Field for whatever is in the group.
// Now the information about the group is in the 3 lightweight classes sparc::hydro_set, sparc::eos_set, sparc::material
struct EOSIdealGas:EOSDataTool
{
	EOSIdealGas(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void ComponentContribution(Field& hydro,Field& eos);
	virtual void ApplyEOS(tw::Int ie, ScalarField& nu_e,Field& hydro, Field& eos);
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
	EOSHotElectrons(const std::string& name,MetricSpace *m,Task *tsk);
	// Missing component function implies we expect hot electrons will always be the lone member in a group.
	// Perhaps this is why we need special single-component versions of ApplyEOS?
	virtual void ApplyEOS(tw::Int ie, ScalarField& nu_e,Field& hydro, Field& eos);
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
	tw::Float DensitySum(const Field& f,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += f(cell,s+hidx.first);
		return ans;
	}
	tw::Float MassDensity(const Field& f,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += f(cell,s+hidx.first)*eosComponents[s]->mat.mass;
		return ans;
	}
	tw::Float MixHeatCapacities(const Field& f,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += f(cell,s+hidx.first)*eosComponents[s]->mat.cvm;
		return ans;
	}
	tw::Float MixVibrationalEnergy(const Field& f,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += f(cell,s+hidx.first)*eosComponents[s]->mat.excitationEnergy;
		return ans;
	}
	tw::Float MixVibrationalStates(const Field& f,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += eosComponents[s]->mat.excitationEnergy > 0.0 ? f(cell,s+hidx.first) : 0.0;
		return ans;
	}

	EOSMixture(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void ApplyEOS(tw::Int ie, ScalarField& nu_e,Field& hydro, Field& eos);
};

// DFG - refined input file control eliminates need for special ideal gas mix class
// So now we keep it here but just as a copy of the base mix class.

struct EOSIdealGasMix:EOSMixture
{
	// it seems nothing is different, unless we really want to override the constituents.
	EOSIdealGasMix(const std::string& name,MetricSpace *m,Task *tsk);
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

	EOSMieGruneisen(const std::string& name,MetricSpace *m,Task *tsk);

	virtual void ComponentContribution(Field& hydro, Field& eos);

	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
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

	EOSMieGruneisen2(const std::string& name,MetricSpace *m,Task *tsk);

	virtual void ComponentContribution(Field& hydro, Field& eos);

	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};
