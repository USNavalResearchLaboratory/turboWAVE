module;

#include "tw_includes.h"

export module physics;

import input;
import compute_tool;
import fields;
import functions;

export namespace sparc
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

	tw::Float CoulombCrossSectionCGS(tw::Float q1,tw::Float q2,tw::Float m12,tw::Float v12,tw::Float N1,tw::Float N2,tw::Float T1,tw::Float T2);
	tw::Float ElectronPhononFrequencyCGS(tw::Float Ti,tw::Float EFermi,tw::Float ks);
}

export struct Ionizer : ComputeTool
{
	// this tool requires the owner to manage indexing of ionized species
	tw::Float ionizationPotential;
	tw::Float electrons,protons;
	tw::Float multiplier,max_rate,cutoff_field;

	// members determining species involved
	std::string ion_name,electron_name;
	sparc::hydro_set hi,he,hgas; // for hydro save the field indexing

	// members that are assigned in Initialize or in constructor
	tw::Float Z,Uion,nstar,lstar,l,m,I1,I2,I3,A1,A2,A3;

	Ionizer(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	void DeduceAveragedCoeff();
	void DeduceStaticCoeff();
	virtual tw::Float InstantRate(tw::Float w0,tw::Float E) { return 0.0; }
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E) { return 0.0; }
	//tw::Float ThresholdEstimate() { return A3*tw::dims::electric_field >> atomic >> native; }
};

export struct Multiphoton : Ionizer
{
	tw::Float E_MPI;
	Multiphoton(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E);
};

export struct Tunneling : Ionizer
{
	// abstract class
	Tunneling(const std::string& name,MetricSpace *m,Task *tsk) : Ionizer(name,m,tsk) {}
	virtual tw::Float InstantRate(tw::Float w0,tw::Float E);
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E);
};

export struct KYH : Tunneling
{
	KYH(const std::string& name,MetricSpace *m,Task *tsk) : Tunneling(name,m,tsk) {}
	virtual void Initialize();
};

export struct ADK : Tunneling
{
	ADK(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual bool Test(tw::Int& id);
};

export struct PPT_Tunneling : ADK
{
	PPT_Tunneling(const std::string& name,MetricSpace *m,Task *tsk) : ADK(name,m,tsk) {}
	virtual void Initialize();
	virtual bool Test(tw::Int& id);
};

export struct PPT : Ionizer
{
	tw::Int terms;
	PPT(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E);
	tw::Float FourierSum(tw::Float gam,tw::Float nu);
};

export struct PMPB : PPT
{
	PMPB(const std::string& name,MetricSpace *m,Task *tsk) : PPT(name,m,tsk) {}
	virtual void Initialize();
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E);
};

// notes on EOS:
// Higher levels must orchestrate the following sequence (see Chemistry::ApplyEOS):
// 1. Components load the nmcv array
// 2. Mixture computes T, and returns IE and nm for components to use
// 3. Components load the P, K, and visc arrays

// EOS classes retain their own indexing and material information in 3 lightweight structs:
// sparc::hydro_set and sparc::eos_set contain indices into a Field object
// sparc::material contains floats characterizing a chemical
export struct EOSComponent:ComputeTool
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

export struct EOSIdealGas:EOSComponent
{
	EOSIdealGas(const std::string& name,MetricSpace *m,Task *tsk);
};

// Braginskii model for plasma electrons, hot enough to ignore quantum effects
export struct EOSHotElectrons:EOSComponent
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
export struct EOSSimpleMieGruneisen:EOSComponent
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
export struct EOSLinearMieGruneisen:EOSComponent
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
export struct EOSTillotson:EOSComponent
{
	tw::Float rho0;   // Reference density

	tw::Float a;   // Tillotson Coefficient
	tw::Float b;   // Tillotson Coefficient
	tw::Float A;   // Bulk Modulus [pressure]
	tw::Float B;   // Tillotson Parameter [pressure]
	tw::Float alpha;   // Tillotson Coefficient
	tw::Float beta;   // Tillotson Coefficient

	tw::Float rhoIV;   // Incipient vaporization density
	tw::Float E0;   // Reference specific energy
	tw::Float EIV;   // Incipient vaporization specific energy
	tw::Float ECV;   // Complete vaporization specific energy

	EOSTillotson(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void AddPKV(ScalarField& IE,ScalarField& nm,ScalarField& nu_e,Field& hydro,Field& eos);
};

// DFG - the role of the mixture is still the same, it computes non-additive things like temperature.
// The components are then called at higher levels to add in additive things like pressure.
export struct EOSMixture:ComputeTool
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

export struct EOSIdealGasMix:EOSMixture
{
	// This uses the polytropic ideal gas caloric EOS in its original form.
	// However, see comments in UpdateEnergy.
	EOSIdealGasMix(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void ComputeTemperature(ScalarField& IE,ScalarField& nm,Field& hydroRef,Field& hydro,Field& eosRef,Field& eos);
};

// DFG - Populating the input hash is partly delegated to this low level material object.
// We are now requiring the user to specify units, otherwise normalized units are assumed.
// *This breaks old input files*, but is necessary for the sake of consistency.
void sparc::material::AddDirectives(tw::input::DirectiveReader& directives)
{
	directives.Add("mass",new tw::input::Float(&mass));
	directives.Add("charge",new tw::input::Float(&charge));
	directives.Add("cv",new tw::input::Float(&cvm));
	directives.Add("vibrational energy",new tw::input::Float(&excitationEnergy),false);
	directives.Add("thermometric conductivity",new tw::input::Float(&thermometricConductivity),false);
	directives.Add("kinematic viscosity",new tw::input::Float(&kinematicViscosity),false);
	directives.Add("permittivity",new tw::input::Numbers<tw::Float>(&eps[0],2),false);
}

tw::Float sparc::CoulombCrossSectionCGS(tw::Float q1,tw::Float q2,tw::Float m12,tw::Float v12,tw::Float N1,tw::Float N2,tw::Float T1,tw::Float T2)
{
	// temperatures in ergs
	tw::Float rmin,rmax,rmin_alt,coulombLog;
	rmin = std::fabs(q1*q2)/(tw::small_pos + m12*v12*v12);
	rmin_alt = cgs::hbar/(tw::small_pos + 2*m12*v12);
	rmin = (rmin*rmin_alt)/(rmin+rmin_alt); // efficiently estimate the smaller of the two
	rmax = 1/std::sqrt(tw::small_pos + 4*pi*N1*q1*q1/T1 + 4*pi*N2*q2*q2/T2);
	coulombLog = std::log((3*rmin+rmax)/rmin); // well behaved approximation of std::log(rmax/rmin)
	return (32/pi)*pow(v12,-4)*sqr(q1*q2/m12)*coulombLog;
}

tw::Float sparc::ElectronPhononFrequencyCGS(tw::Float Ti,tw::Float EFermi,tw::Float ks)
{
	// temperatures in ergs
	const tw::Float vF = std::sqrt(2*EFermi/cgs::me);
	// collision frequency in cgs units
	return 2*ks*sqr(cgs::qe)*Ti/(sqr(cgs::hbar)*vF);
}

//////////////////////////
//                      //
//   Photoionization    //
//                      //
//////////////////////////


Ionizer::Ionizer(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk)
{
	ionizationPotential = 1e-5; // in native units
	cutoff_field = 1e-3; // in atomic units
	electrons = 0;
	protons = 0;
	multiplier = 1.0;
	I1 = I2 = I3 = A1 = A2 = A3 = 0.0;
	nstar = 1.0;
	lstar = l = m = 0.0;
	max_rate = tw::big_pos;
	// read ionspecies and electronspecies indices in Species::Initialize
	// setup hydro indexing during Chemical::Initialize
	directives.Add("ionization potential",new tw::input::Float(&ionizationPotential));
	directives.Add("protons",new tw::input::Float(&protons));
	directives.Add("electrons",new tw::input::Float(&electrons));
	directives.Add("multiplier",new tw::input::Float(&multiplier),false);
	directives.Add("saturated rate",new tw::input::Float(&max_rate),false);
	directives.Add("ion species",new tw::input::String(&ion_name));
	directives.Add("electron species",new tw::input::String(&electron_name));
}

void Ionizer::Initialize()
{
	ComputeTool::Initialize();
	Z = protons - electrons + 1;
	Uion = ionizationPotential * tw::dims::energy >> native >> atomic;
	nstar = Z / std::sqrt(2*Uion);
}

void Ionizer::DeduceStaticCoeff()
{
	I1 = A1 * pow(2*Uion,0.75) / std::sqrt(3/pi);
	I2 = A2 - 0.5;
	I3 = A3;
}

void Ionizer::DeduceAveragedCoeff()
{
	A1 = I1 * std::sqrt(3/pi) / pow(2*Uion,0.75);
	A2 = I2 + 0.5;
	A3 = I3;
}

Multiphoton::Multiphoton(const std::string& name,MetricSpace *m,Task *tsk) : Ionizer(name,m,tsk)
{
	directives.Add("reference field",new tw::input::Float(&E_MPI));
}

void Multiphoton::Initialize()
{
	A1 = multiplier*two*pi;
}

tw::Float Multiphoton::AverageRate(tw::Float w0,tw::Float E)
{
	const tw::Float wa = w0 * tw::dims::angular_frequency >> native >> atomic;
	const tw::Float photons = MyFloor(Uion/wa + 1);
	return A1*w0*pow(std::fabs(E)/E_MPI,two*photons) / Factorial(photons-1);
}

tw::Float Tunneling::InstantRate(tw::Float w0,tw::Float E)
{
	const tw::Float Ea = (std::fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
	return I1*pow(Ea,I2)*std::exp(I3/Ea);
}

tw::Float Tunneling::AverageRate(tw::Float w0,tw::Float E)
{
	const tw::Float Ea = (std::fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
	return A1*pow(Ea,A2)*std::exp(A3/Ea);
}

void KYH::Initialize()
{
	Ionizer::Initialize();
	const tw::Float alpha = 0.0072973525693;
	const tw::Float Ua2 = Uion*alpha*alpha;
	const tw::Float Ea = pow(2*Uion,1.5);
	A1 = multiplier * pow(2.0,3.0-4*Ua2) * std::sqrt(3/pi) * (1.0-7*Ua2/72) * std::exp(4*Ua2) * pow(2*Uion,1.75-3.0*Ua2);
	A1 /= std::tgamma(3.0-2*Ua2);
	A2 = 2.0*Ua2 - 0.5;
	A3 = -(2.0/3.0)*Ea*(1.0-Ua2/12);
	DeduceStaticCoeff();
	I1 = I1*tw::dims::angular_frequency >> atomic >> native;
	A1 = A1*tw::dims::angular_frequency >> atomic >> native;
}

ADK::ADK(const std::string& name,MetricSpace *ms,Task *tsk) : Tunneling(name,ms,tsk)
{
	directives.Add("orbital number",new tw::input::Float(&l),false);
	directives.Add("orbital projection",new tw::input::Float(&m),false);
	directives.Add("effective orbital number",new tw::input::Float(&lstar),false);
}

void ADK::Initialize()
{
	Ionizer::Initialize();
	m = std::fabs(m);
	const tw::Float e = std::exp(1.0);
	const tw::Float sn2l2 = std::sqrt(nstar*nstar-lstar*lstar);
	A1 = std::sqrt(3/cub(pi)) * (2*l+1) * std::tgamma(l+m+1) / std::tgamma(m+1) / std::tgamma(l-m+1);
	A1 *= pow(e/sn2l2,m+1.5); // exponent is of poor print quality in JETP
	A1 *= pow((nstar+lstar)/(nstar-lstar),lstar+0.5);
	A1 *= (Z*Z/cub(nstar));
	A1 *= pow(4*e*cub(Z/nstar)/sn2l2,2*nstar-m-1.5);
	A2 = m+1.5-2*nstar;
	A3 = -2*cub(Z/nstar)/3; // ADK 1986 has nstar**4 in Eq. 21, must be a typo?
	DeduceStaticCoeff();
	I1 = I1*tw::dims::angular_frequency >> atomic >> native;
	A1 = A1*tw::dims::angular_frequency >> atomic >> native;
}

void PPT_Tunneling::Initialize()
{
	Ionizer::Initialize();
	m = std::fabs(m);
	const tw::Float F0 = cub(Z/nstar);
	// First without the Coulomb factor
	A1 = Uion;
	A1 *= pow(2,2*nstar) / (nstar*std::tgamma(nstar+lstar+1)*std::tgamma(nstar-lstar)); // |C|^2
	A1 *= std::sqrt(6/pi);
	A1 *= (2*l+1) * std::tgamma(l+m+1) / pow(2,m) / std::tgamma(m+1) / std::tgamma(l-m+1);
	A1 *= pow(0.5/F0,m+1.5);
	A2 = m+1.5;
	A3 = -2*F0/3;
	// Account for Coulomb correction
	A1 *= pow(2*F0,2*nstar);
	A2 -= 2*nstar;
	DeduceStaticCoeff();
	I1 = I1*tw::dims::angular_frequency >> atomic >> native;
	A1 = A1*tw::dims::angular_frequency >> atomic >> native;
}

PPT::PPT(const std::string& name,MetricSpace *m,Task *tsk) : Ionizer(name,m,tsk)
{
	terms = 1;
	directives.Add("terms",new tw::input::Int(&terms));
	directives.Add("orbital number",new tw::input::Float(&l),false);
	directives.Add("effective orbital number",new tw::input::Float(&lstar),false);
	// the projection must be zero, do not accept input for it
}

void PPT::Initialize()
{
	Ionizer::Initialize();
	// use A1 to hold C_nl^2
	A1 = pow(2,2*nstar);
	A1 /= nstar*std::tgamma(nstar+lstar+1)*std::tgamma(nstar-lstar);
	// use A3 to hold F0
	A3 = pow(2*Uion,tw::Float(1.5));
}

tw::Float PPT::FourierSum(tw::Float gam,tw::Float nu)
{
	tw::Float A0 = 0.0;
	const tw::Float alpha = 2*(std::asinh(gam) - gam/std::sqrt(1 + gam*gam));
	const tw::Float beta = 2*gam/std::sqrt(1 + gam*gam);
	const tw::Float dnu = MyCeil(nu)-nu;
	for (tw::Int n=0;n<terms;n++)
		A0 += std::exp(-alpha*(n+dnu))*tw::dawsoni(std::sqrt(beta*(n+dnu)));
	return A0;
}

tw::Float PPT::AverageRate(tw::Float w0,tw::Float E)
{
	tw::Float ans;
	const tw::Float wa = (w0*tw::dims::angular_frequency >> native >> atomic);
	const tw::Float Ea = (std::fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
	const tw::Float F = Ea/A3;
	const tw::Float gam = std::sqrt(2*Uion)*wa/Ea;
	const tw::Float g = (3/(2*gam))*((1 + 1/(2*gam*gam))*std::asinh(gam) - std::sqrt(1 + gam*gam)/(2*gam));
	const tw::Float nu = (Uion/wa) * (1 + 1/(2*gam*gam));
	ans = Uion*A1*std::sqrt(6/pi)*(2*l+1);
	ans *= pow(F*std::sqrt(1 + gam*gam)/2,tw::Float(1.5));
	ans *= (4/std::sqrt(3*pi)) * (gam*gam/(1 + gam*gam)) * FourierSum(gam,nu);
	ans *= std::exp(-2*g/(3*F));
	ans *= pow(2/F,2*nstar); // coulomb correction
	return ans*tw::dims::angular_frequency >> atomic >> native;
}

void PMPB::Initialize()
{
	Ionizer::Initialize();
	// use A1 to hold C_nl^2, n.b. PMPB has a different convention from PPT (factor of 4)
	A1 = pow(2,2*nstar-2);
	A1 /= nstar*std::tgamma(nstar+lstar+1)*std::tgamma(nstar-lstar);
	// use A3 to hold F0
	A3 = pow(2*Uion,tw::Float(1.5));
}

tw::Float PMPB::AverageRate(tw::Float w0,tw::Float E)
{
	// Mappings from our PPT notation to PMPB notation:
	// beta -> beta , gam -> gam , Uion -> I , Uion/w -> K0
	// E/F0 -> F , alpha -> 2*c1 , nu -> nth , g -> g , dawson_integral -> script-F
	tw::Float ans;
	const tw::Float wa = (w0*tw::dims::angular_frequency >> native >> atomic);
	const tw::Float Ea = (std::fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
	const tw::Float F = Ea/A3;
	const tw::Float gam = std::sqrt(2*Uion)*wa/Ea;
	const tw::Float g = (3/(2*gam))*((1 + 1/(2*gam*gam))*std::asinh(gam) - std::sqrt(1 + gam*gam)/(2*gam));
	const tw::Float nu = (Uion/wa) * (1 + 1/(2*gam*gam));
	ans = (2/pi)*Uion*A1*(2*l+1);
	ans *= pow(Uion/wa,-1.5);
	ans *= std::sqrt(2*gam/std::sqrt(1+gam*gam)) * FourierSum(gam,nu);
	ans *= std::exp(-2*g/(3*F));
	// Following is the improved Coulomb correction factor
	ans *= pow(2/F,2*nstar) * pow(1+2*gam/std::exp(1),-2*nstar);
	return ans*tw::dims::angular_frequency >> atomic >> native;
}

// ASHER_MOD

//////////////////////////////////
//                              //
// Equation of State Base Class //
//                              //
//////////////////////////////////


EOSComponent::EOSComponent(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	// DFG - must not set names anymore here, let the base class constructor handle it
	// Note enumerated types are now strongly typed and namespaced.
}

void EOSComponent::SetHeatCapacity(ScalarField& nm,Field& eos)
{
	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
			eos(cell,eidx.nmcv) = nm(cell) * mat.cvm / mat.mass;
	}
}

void EOSComponent::AddHeatCapacity(Field& hydro,Field& eos)
{
	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
			eos(cell,eidx.nmcv) += hydro(cell,hidx.ni) * mat.cvm;
	}
}

void EOSComponent::AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
{
	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
		{
			const tw::Float ngas = hydro(cell,hidx.ni);
			eos(cell,eidx.P) += ngas*eos(cell,eidx.T);
			eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * ngas;
			eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * ngas;
		}
	}
}

////////////////////////////////
//                            //
// Ideal Gas Law for Chemical //
//                            //
////////////////////////////////

// nothing to do, ideal gas is default

EOSIdealGas::EOSIdealGas(const std::string& name,MetricSpace *m, Task *tsk) : EOSComponent(name,m,tsk)
{
}

/////////////////////////////////
//                             //
//    EOS for hot electrons    //
//                             //
/////////////////////////////////

EOSHotElectrons::EOSHotElectrons(const std::string& name,MetricSpace *m, Task *tsk) : EOSComponent(name,m,tsk)
{
}

void EOSHotElectrons::AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
{
	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
		{
			const tw::Float ne = hydro(cell,hidx.ni);
			eos(cell,eidx.P) += ne*eos(cell,eidx.T);
			eos(cell,eidx.K) += 3.2*ne*eos(cell,eidx.T)/(mat.mass*nu_e(cell));
			//eos(cell,eidx.visc) += 0.0; // don't touch, may help caching.
			// Braginskii has for e-viscosity 0.73*ne*eos(cell,eidx.T)/nu_e(cell)
			// However, we are forcing electrons to move with ions and so should not diffuse velocity field
		}
	}
}

//////////////////////////////////
//                              //
// EOS Simplified Mie Gruneisen //
//                              //
//////////////////////////////////

EOSSimpleMieGruneisen::EOSSimpleMieGruneisen(const std::string& name,MetricSpace *m, Task *tsk) : EOSComponent(name,m,tsk)
{
	GRUN = 2.0; // value for Cu on p. 257 of "Shock Wave Physics and Equation of State Modeling"
	// GRUN = 0.1; // value for water in the above book.
	directives.Add("gruneisen parameter",new tw::input::Float(&GRUN));
}

void EOSSimpleMieGruneisen::AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
{
	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
		{
			const tw::Float nion = hydro(cell,hidx.ni);
			const tw::Float partial_IE = IE(cell) * nion * mat.mass / nm(cell);
			eos(cell,eidx.P) += GRUN*partial_IE;
			eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * nion;
			eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * nion;
		}
	}
}


/////////////////////////
//                     //
//  EOS Mie Gruneisen  //
//                     //
/////////////////////////
// better MieGruneisen EOS model that uses a linear Hugoniot fit

EOSLinearMieGruneisen::EOSLinearMieGruneisen(const std::string& name,MetricSpace *m, Task *tsk) : EOSComponent(name,m,tsk)
{
	// Hugoniot data fit for Cu
	n0 = 3.3e3;
	c0 = 1.3248e-5;
	S1 = 1.5;

	// Hugoniot data fit for H20
	// n0 = 1334.0;
	// c0 = 5.197e-6;
	// S1 = 1.8153;

	GRUN = 2.0; // value for Cu on p. 257 of "Shock Wave Physics and Equation of State Modeling"
	// GRUN = 0.1; // value for water in the above book.

	directives.Add("gruneisen parameter",new tw::input::Float(&GRUN));
	directives.Add("reference density",new tw::input::Float(&n0));
	directives.Add("hugoniot intercept",new tw::input::Float(&c0));
	directives.Add("hugoniot slope", new tw::input::Float(&S1));
}

void EOSLinearMieGruneisen::AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
{
	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
		{
			const tw::Float nion = hydro(cell,hidx.ni);
			const tw::Float partial_IE = IE(cell) * nion * mat.mass / nm(cell);
			const tw::Float mu = nion/n0 - 1;
			const tw::Float sel = tw::Float(mu>=0.0);
			// Temperature is not treated as additive, worked out by parent object
			// Pressure from <http://bluevistasw.com/2016/02/16/mie-gruneisen-eos-implementation/>
			eos(cell,eidx.P) += (1-sel)*(mat.mass*n0*c0*c0*mu + GRUN*(mu+1)*partial_IE);
			eos(cell,eidx.P) += sel*(mat.mass*n0*c0*c0*mu*(1 + (1 - GRUN*(mu+1)/2)*mu)/(1 - (S1-1)*mu) + GRUN*(mu+1)*partial_IE);
			eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * nion;
			eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * nion;
		}
	}
}

//////////////////////
//                  //
//  EOS Tillotson   //
//                  //
//////////////////////
// EOS for shock and cavitation
// for materials such as water

EOSTillotson::EOSTillotson(const std::string& name,MetricSpace *m, Task *tsk) : EOSComponent(name,m,tsk)
{
	// Tillotson parameters for H20
	rho0 = tw::dnum("0.998 [g/cm3]") >> native;

	a = 0.7;   // Tillotson Coefficient
	b = 0.15;   // Tillotson Coefficient
	A = tw::dnum("21.8e3 [bar]") >> native;   // Tillotson Coefficient, pressure
	B = tw::dnum("132.5e3 [bar]") >> native;   // Tillotson Coefficient, pressure
	alpha = 10.0;   // Tillotson Coefficient
	beta = 5.0;   // Tillotson Coefficient

	rhoIV = tw::dnum("0.958 [g/cm3]") >> native;   // Incipient vaporization density
	E0 = tw::dnum("0.07e12 [ergs/g]") >> native;   // Reference energy
	EIV = tw::dnum("0.00419e12 [ergs/g]") >> native;   // Incipient vaporization specific energy
	ECV = tw::dnum("0.025e12 [ergs/g]") >> native;   // Complete vaporization specific energy

	directives.Add("reference mass density",new tw::input::Float(&rho0));

	directives.Add("parameter a",new tw::input::Float(&a));
	directives.Add("parameter b",new tw::input::Float(&b));
	directives.Add("pressure A",new tw::input::Float(&A));
	directives.Add("pressure B",new tw::input::Float(&B));
	directives.Add("parameter alpha",new tw::input::Float(&alpha));
	directives.Add("parameter beta",new tw::input::Float(&beta));

	directives.Add("incipient vaporization mass density",new tw::input::Float(&rhoIV));
	directives.Add("reference specific energy",new tw::input::Float(&E0));
	directives.Add("incipient vaporization specific energy",new tw::input::Float(&EIV));
	directives.Add("complete vaporization specific energy",new tw::input::Float(&ECV));
}

// There are four (arguably 5) regions in the revised Tillotson EOS
// [A.L. Brundage, Procedia Engineering (2013)]
//
// (1) Compressed States           : (\rho > \rho_0 & E > 0)
// (2) Cold Expanded States        : (\rho_0 > \rho > \rho_IV & E < E_{IV})
// (3) Hot Expanded States         : (\rho_0 > \rho & E >= E_{CV})
// (4) Low Energy Expansion States : (\rho < \rho_IV & E < E_{CV})
// (5) Mixed Region                : ( \rho_0 > \rho > \rho_IV & E_{CV} > E > E_{IV} )
void EOSTillotson::AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
{
	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
		{
			const tw::Float ndens = hydro(cell,hidx.ni);
			const tw::Float u = hydro(cell,hidx.u);
			const tw::Float rho = mat.mass * ndens;
			const tw::Float u0 = rho*E0 + tw::small_pos;

			const tw::Float eta = rho/rho0; // compression
			const tw::Float mew = eta - 1.0; // strain

			// Determine Region
			tw::Int region = 0; // 1, 2, 3, 4, or 5 : 0 is for error detection
			if ( rho >= rho0 ) region = 1;
			if ( ((rho >= rhoIV) && (rho < rho0)) && (u <= rho*EIV) ) region = 2;
			if ((rho < rho0) && (u >= rho*ECV)) region = 3;
			if ((rho < rhoIV) && (u < rho*ECV)) region = 4;
			if (((rho > rhoIV) && (rho < rho0)) && ((u > rho*EIV) && (u < rho*ECV))) region = 5;

			if (region == 0) {
				std::stringstream err_mess;
				err_mess << "Unrecognized Region Detected in Tillotson EOS." << std::endl;
				err_mess << "rho = " << (rho*tw::dims::mass_density>>native>>cgs) << " [g/cm3]" << std::endl;
				err_mess << "E = " << ((u/rho)*tw::dims::specific_energy>>native>>cgs) << " [ergs/g]" << std::endl;
				throw tw::FatalError(err_mess.str());
			}

			// Pressure Calculation
			tw::Float denom = (u/(u0*sqr(eta))) + 1.0; // this quantity is repeated in expressions
			tw::Float expo, P2, P3; // might be used depending on state
			switch (region)
			{
				case 1:
					eos(cell,eidx.P) += (a + b/denom)*u + A*mew + B*sqr(mew);
					break;
				case 2:
					eos(cell,eidx.P) += (a + b/denom)*u + A*mew + B*sqr(mew);
					break;
				case 3:
					expo = (rho0/rho) - 1.0;
					eos(cell,eidx.P) += a*u + ((b*u/denom) + A*mew*std::exp(-beta*expo))*std::exp(-alpha*sqr(expo));
					break;
				case 4:
					eos(cell,eidx.P) += (a + b/denom)*u + A*mew;
					break;
				case 5: // this is an interpolation of region 2 and 3
					P2 = (a + b/denom)*u + A*mew + B*sqr(mew);
					expo = (rho0/rho) - 1.0;
					P3 = a*u + ((b*u/denom) + A*mew*std::exp(-beta*expo))*std::exp(-alpha*sqr(expo));
					eos(cell,eidx.P) += ((u - rho*EIV)*P3 + (rho*ECV - u)*P2)/(rho*(ECV-EIV));
					break;
			}

			// keep rest of EOS the same
			eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * ndens;
			eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * ndens;
		}
	}

}

///////////////////////////////////////
//                                   //
// EOS Mixture for EquilibriumGroup  //
//                                   //
///////////////////////////////////////

EOSMixture::EOSMixture(const std::string& name,MetricSpace *m, Task *tsk) : ComputeTool(name,m,tsk)
{
}

void EOSMixture::ComputeTemperature(ScalarField& IE, ScalarField& nm, Field& hydro, Field& eos)
{
	// DFG - caloric EOS with zero-reference (use if only one time level is available)
	// Compute the temperature assuming nmcv has been loaded.
	// Pass IE and nm back out for use by component EOS classes

	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
		{
			nm(cell) = MassDensity(hydro,cell);
			IE(cell) = InternalEnergy(nm(cell),hydro,cell);
			const tw::Float epsvn = MixVibrationalEnergy(hydro,cell);
			const tw::Float nv = MixVibrationalStates(hydro,cell);

			eos(cell,eidx.T) = IE(cell)/(tw::small_pos + eos(cell,eidx.nmcv));
			eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/std::log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
		}
	}
}

void EOSMixture::ComputeTemperature(ScalarField& IE, ScalarField& nm, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos)
{
	// DFG - caloric EOS accounting for nearby reference state (usually previous time level)
	// Compute the temperature assuming nmcv has been loaded.
	// Pass IE and nm back out for use by component EOS classes

	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
		{
			const tw::Float nm1 = tw::small_pos + MassDensity(hydro,cell);
			const tw::Float IE1 = InternalEnergy(nm1,hydro,cell);
			const tw::Float nm0 = tw::small_pos + MassDensity(hydroRef,cell);
			const tw::Float IE0 = InternalEnergy(nm0,hydroRef,cell);
			const tw::Float epsvn = MixVibrationalEnergy(hydro,cell);
			const tw::Float nv = MixVibrationalStates(hydro,cell);
			const tw::Float nmcv_sum = tw::small_pos + eosRef(cell,eidx.nmcv) + eos(cell,eidx.nmcv);

			//eos(cell,eidx.T) = eosRef(cell,eidx.T) + 2.0*(IE1 - IE0)/nmcv_sum; // wrong
			eos(cell,eidx.T) = eosRef(cell,eidx.T) + (IE1*(1+nm0/nm1) - IE0*(1+nm1/nm0))/nmcv_sum;
			eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/std::log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
			nm(cell) = nm1;
			IE(cell) = IE1;
		}
	}
}

void EOSMixture::UpdateEnergy(ScalarField& nm,ScalarField& T0,Field& hydro,Field& eos)
{
	// Add energy corresponding to a change in temperature only.
	// Not centered, because cv is not updated.
	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
			hydro(cell,hidx.u) += eos(cell,eidx.nmcv) * (eos(cell,eidx.T) - T0(cell));
	}
}

///////////////////////////
//                       //
//     Ideal Gas Mix     //
//                       //
///////////////////////////

EOSIdealGasMix::EOSIdealGasMix(const std::string& name,MetricSpace *m, Task *tsk) : EOSMixture(name,m,tsk)
{
}

void EOSIdealGasMix::ComputeTemperature(ScalarField& IE, ScalarField& nm, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos)
{
	// DFG - polytropic ideal gas caloric EOS in the original SPARC mode of calculation.  Ignores reference states.
	// (does not force pressure to be ideal gas)
	// Compute the temperature assuming nmcv has been loaded.
	// Pass IE and nm back out for use by component EOS classes

	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
		{
			const tw::Float nm1 = MassDensity(hydro,cell);
			const tw::Float IE1 = InternalEnergy(nm1,hydro,cell);
			const tw::Float epsvn = MixVibrationalEnergy(hydro,cell);
			const tw::Float nv = MixVibrationalStates(hydro,cell);

			eos(cell,eidx.T) = IE1/(tw::small_pos + eos(cell,eidx.nmcv));
			eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/std::log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
			nm(cell) = nm1;
			IE(cell) = IE1;
		}
	}
}
