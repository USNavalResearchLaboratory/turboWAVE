namespace sparc
{
	struct hydro_set
	{
		// num = number of constituent chemicals
		// first = index of first hydro quantity in the set, must be the first chemical density
		// last = index of last hydro quantity in the set
		// ni = index to a particular chemical density (optionally used)
		// npx,npy,npz = index to momentum density
		// u = index to energy density
		// x = index to vibrational density
		tw::Int num,first,last,ni,npx,npy,npz,u,x;
		tw::Int Load(tw::Int i,tw::Int N)
		{
			num = N;
			ni = i;
			first = i;
			npx = i+N;
			npy = i+N+1;
			npz = i+N+2;
			u = i+N+3;
			x = i+N+4;
			last = x;
			return last+1;
		}
		tw::Float DensitySum(const Field& f,const tw::cell& cell) const
		{
			tw::Float ans = 0.0;
			for (tw::Int i=0;i<num;i++)
				ans += f(cell,first+i);
			return ans;
		}
		static const tw::Int count = 5; // gives count of state components not counting densities
	};

	struct eos_set
	{
		// T = index to translational and rotational temperature
		// Tv = vibrational temperature (if used)
		// P = pressure
		// K = heat conductivity
		// visc = dynamic viscosity
		// nmcv = mass density weighted heat capacity at constant volume
		tw::Int T,Tv,P,K,visc,nmcv;
		tw::Int Load(tw::Int i) { T=i;Tv=i+1;P=i+2;K=i+3;visc=i+4;nmcv=i+5;return i+6; }
		static const tw::Int count = 6;
	};

	struct material
	{
		// Characteristics of a chemical
		tw::Float mass,charge,cvm,excitationEnergy,thermometricConductivity,kinematicViscosity,eps[2];
		void AddDirectives(tw::input::DirectiveReader& directives);
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
			eps_r[i] = mat.eps[0];
			eps_i[i] = mat.eps[1];
		}
	};

	tw::Float CoulombCrossSection(const UnitConverter& uc,tw::Float q1,tw::Float q2,tw::Float m12,tw::Float v12,tw::Float N1,tw::Float N2,tw::Float T1,tw::Float T2);
	tw::Float ElectronPhononRateCoeff(const UnitConverter& uc,tw::Float Ti,tw::Float EFermi,tw::Float ks,tw::Float nref);
}

struct Ionizer : ComputeTool
{
	// this tool requires the owner to manage indexing of ionized species
	tw::Float ionizationPotential;
	tw::Float electrons,protons;
	tw::Float multiplier,max_rate;

	// members determining species involved
	std::string ion_name,electron_name;
	sparc::hydro_set hi,he,hgas; // for hydro save the field indexing

	// members that are assigned in Initialize or in constructor
	tw::Float Z,Uion,nstar,lstar,I1,I2,I3,A1,A2,A3;

	Ionizer(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float InstantRate(tw::Float w0,tw::Float E) { return 0.0; }
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E) { return 0.0; }
	tw::Float ThresholdEstimate() { return space->units->ConvertToNative(A3,tw::dimensions::electric_field,tw::units::atomic); }
};

struct Multiphoton : Ionizer
{
	tw::Float E_MPI;
	Multiphoton(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E);
};

struct ADK : Ionizer
{
	ADK(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float InstantRate(tw::Float w0,tw::Float E);
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E);
};

struct PPT : Ionizer
{
	tw::Int terms;
	PPT(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E);
	tw::Float wfunc(tw::Float x);
};

// DFG - redesigned to:
// (i) take advantage of encapsulated data structures
// (ii) conform to new ComputeTool spec (strong preference that containment tree be reserved for modules)
// In this process I split components and mixtures more strongly.
// Higher levels must orchestrate the following sequence (see Chemistry::ApplyEOS):
// 1. Components load the nmcv array
// 2. Mixture computes T, and returns IE and nm for components to use
// 3. Components load the P, K, and visc arrays

// EOS classes now retain their own indexing and material information in 3 lightweight structs:
// sparc::hydro_set and sparc::eos_set contain indices into a Field object
// sparc::material contains floats characterizing a chemical
struct EOSComponent:ComputeTool
{
	sparc::hydro_set hidx;
	sparc::eos_set eidx;
	sparc::material mat;

	EOSComponent(const std::string& name,MetricSpace *m,Task *tsk);
	void SetupIndexing(tw::Int component_index,const sparc::hydro_set& h,const sparc::eos_set& e,const sparc::material& m)
	{
		hidx = h;
		hidx.ni = component_index;
		eidx = e;
		mat = m;
	}

	virtual void SetHeatCapacity(ScalarField& nm,Field& eos);
	virtual void AddHeatCapacity(Field& hydro,Field& eos);
	virtual void AddPKV(ScalarField& IE,ScalarField& nm,ScalarField& nu_e,Field& hydro,Field& eos);
};

struct EOSIdealGas:EOSComponent
{
	EOSIdealGas(const std::string& name,MetricSpace *m,Task *tsk);
};

// Braginskii model for plasma electrons, hot enough to ignore quantum effects
struct EOSHotElectrons:EOSComponent
{
	EOSHotElectrons(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void AddPKV(ScalarField& IE,ScalarField& nm,ScalarField& nu_e,Field& hydro,Field& eos);
};

// This is the most basic implementation of the MieGruneisen EOS
//
// It assumes that GRUN = GRUN0 = const. at all times
// This is typically not used in literature concerning MieGruneisen EOSs,
// as the results are rarely physicaly accurate
// If you produce sound waves with this model (for example with Cu),
// you'll notice the sound speed is off. Qualitatively, it gives a broad picture.
struct EOSSimpleMieGruneisen:EOSComponent
{
	tw::Float GRUN; // Gruneisen coefficient

	EOSSimpleMieGruneisen(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void AddPKV(ScalarField& IE,ScalarField& nm,ScalarField& nu_e,Field& hydro,Field& eos);
};

// This is a MieGruneisen EOS that assumes \rho * GRUN = const., and a linear Hugoniot fit
//
// This is what is more typically what is found in literature regarding MieGruneisen models
// You'll get a more accurate sound speed and shock speeds, assuming that the simulation is
// within the range of relevant Hugoniot data and model assumptions are properly met
struct EOSLinearMieGruneisen:EOSComponent
{
	tw::Float GRUN; // Gruneisen coefficient
	tw::Float n0;   // Reference density
	tw::Float c0;   // y - intercept of Hugoniot fit (usually appriximately speed of sound)
	tw::Float S1;   // coefficient of linear fit of Hugoniot data

	EOSLinearMieGruneisen(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void AddPKV(ScalarField& IE,ScalarField& nm,ScalarField& nu_e,Field& hydro,Field& eos);
};

// Tillotson EOS for modeling vaporization, cavitation, and shocks.
// Coefficients for water can be found at [A.L. Brundage, Procedia Engineering (2013)]
struct EOSTillotson:EOSComponent
{
	tw::Float n0;   // Reference density

	tw::Float a;   // Tillotson Coefficient
	tw::Float b;   // Tillotson Coefficient
	tw::Float A;   // Bulk Modulus [sim]
	tw::Float B;   // Tillotson Parameter [sim]
	tw::Float alpha;   // Tillotson Coefficient
	tw::Float beta;   // Tillotson Coefficient

	tw::Float nIV;   // Vaporization Pressure
	tw::Float E0;   // Reference energy 
	tw::Float EIV;   // Vaporization Energy
	tw::Float ECV;   // Cavitation Energy

	EOSTillotson(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void AddPKV(ScalarField& IE,ScalarField& nm,ScalarField& nu_e,Field& hydro,Field& eos);
};

// DFG - the role of the mixture is still the same, it computes non-additive things like temperature.
// The components are then called at higher levels to add in additive things like pressure.
struct EOSMixture:ComputeTool
{
	sparc::hydro_set hidx;
	sparc::eos_set eidx;
	sparc::material_set matset;
	void SetupIndexing(const sparc::hydro_set& h,const sparc::eos_set& e,const sparc::material_set m)
	{
		hidx = h;
		eidx = e;
		matset = m;
	}
	tw::Float DensitySum(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += f(cell,s+hidx.first);
		return ans;
	}
	tw::Float MassDensity(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += f(cell,s+hidx.first)*matset.mass[s];
		return ans;
	}
	tw::Float MixVibrationalEnergy(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += f(cell,s+hidx.first)*matset.excitationEnergy[s];
		return ans;
	}
	tw::Float MixVibrationalStates(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += matset.excitationEnergy[s] > 0.0 ? f(cell,s+hidx.first) : 0.0;
		return ans;
	}
	tw::Float InternalEnergy(const tw::Float& nm,const Field& f,const tw::cell& cell)
	{
		const tw::vec3 np = tw::vec3(f(cell,hidx.npx),f(cell,hidx.npy),f(cell,hidx.npz));
		const tw::vec3 vel = np/(tw::small_pos + nm);
		const tw::Float KE = 0.5*nm*Norm(vel);
		const tw::Float primitive = f(cell,hidx.u) - KE;
		const tw::Float failsafe = 1e-6*KE + tw::small_pos;
		const tw::Float sel = tw::Float(primitive>0.0);
		return sel*primitive + (1.0-sel)*failsafe;
	}

	EOSMixture(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void ComputeTemperature(ScalarField& IE,ScalarField& nm,Field& hydro,Field& eos);
	virtual void ComputeTemperature(ScalarField& IE,ScalarField& nm,Field& hydroRef,Field& hydro,Field& eosRef,Field& eos);
	virtual void UpdateEnergy(ScalarField& nm,ScalarField& T0,Field& hydro,Field& eos);
};

struct EOSIdealGasMix:EOSMixture
{
	// This uses the polytropic ideal gas caloric EOS in its original form.
	// However, see comments in UpdateEnergy.
	EOSIdealGasMix(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void ComputeTemperature(ScalarField& IE,ScalarField& nm,Field& hydroRef,Field& hydro,Field& eosRef,Field& eos);
};
