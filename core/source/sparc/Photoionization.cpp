module;

#include "tw_includes.h"
#include "tw_test.h"

export module photoionization;

import input;
import driver;
import fields;
import functions;
import hydro_primitives;

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

	Ionizer(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk) {
		ionizationPotential = 1e-5; // in native units
		cutoff_field = 1e-3; // in atomic units
		electrons = 0;
		protons = 0;
		multiplier = 1.0;
		I1 = I2 = I3 = A1 = A2 = A3 = 0.0;
		nstar = 1.0;
		lstar = l = m = 0.0;
		max_rate = tw::max_pos;
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

	virtual void Initialize() {
		ComputeTool::Initialize();
		Z = protons - electrons + 1;
		Uion = ionizationPotential * tw::dims::energy >> native >> atomic;
		nstar = Z / std::sqrt(2*Uion);
	}

	void DeduceAveragedCoeff() {
		A1 = I1 * std::sqrt(3/pi) / std::pow(2*Uion,0.75);
		A2 = I2 + 0.5;
		A3 = I3;
	}

	void DeduceStaticCoeff() {
		I1 = A1 * std::pow(2*Uion,0.75) / std::sqrt(3/pi);
		I2 = A2 - 0.5;
		I3 = A3;
	}

	virtual tw::Float InstantRate(tw::Float w0,tw::Float E) { return 0.0; }
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E) { return 0.0; }
	//tw::Float ThresholdEstimate() { return A3*tw::dims::electric_field >> atomic >> native; }
};

export struct Multiphoton : Ionizer
{
	tw::Float E_MPI;
	Multiphoton(const std::string& name,MetricSpace *m,Task *tsk)  : Ionizer(name,m,tsk) {
		directives.Add("reference field",new tw::input::Float(&E_MPI));
	}
	virtual void Initialize() {
		A1 = multiplier*two*pi;
	}
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E) {
		const tw::Float wa = w0 * tw::dims::angular_frequency >> native >> atomic;
		const tw::Float photons = MyFloor(Uion/wa + 1);
		return A1*w0*std::pow(std::fabs(E)/E_MPI,two*photons) / Factorial(photons-1);
	}
};

export struct Tunneling : Ionizer
{
	// abstract class
	Tunneling(const std::string& name,MetricSpace *m,Task *tsk) : Ionizer(name,m,tsk) {}
	virtual tw::Float InstantRate(tw::Float w0,tw::Float E) {
		const tw::Float Ea = (std::fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
		return I1*std::pow(Ea,I2)*std::exp(I3/Ea);
	}
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E) {
		const tw::Float Ea = (std::fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
		return A1*std::pow(Ea,A2)*std::exp(A3/Ea);
	}
};

export struct KYH : Tunneling
{
	KYH(const std::string& name,MetricSpace *m,Task *tsk) : Tunneling(name,m,tsk) {}
	virtual void Initialize() {
		Ionizer::Initialize();
		const tw::Float alpha = 0.0072973525693;
		const tw::Float Ua2 = Uion*alpha*alpha;
		const tw::Float Ea = std::pow(2*Uion,1.5);
		A1 = multiplier * std::pow(2.0,3.0-4*Ua2) * std::sqrt(3/pi) * (1.0-7*Ua2/72) * std::exp(4*Ua2) * std::pow(2*Uion,1.75-3.0*Ua2);
		A1 /= std::tgamma(3.0-2*Ua2);
		A2 = 2.0*Ua2 - 0.5;
		A3 = -(2.0/3.0)*Ea*(1.0-Ua2/12);
		DeduceStaticCoeff();
		I1 = I1*tw::dims::angular_frequency >> atomic >> native;
		A1 = A1*tw::dims::angular_frequency >> atomic >> native;
	}
};

export struct ADK : Tunneling
{
	ADK(const std::string& name,MetricSpace *ms,Task *tsk) : Tunneling(name,ms,tsk) {
		directives.Add("orbital number",new tw::input::Float(&l),false);
		directives.Add("orbital projection",new tw::input::Float(&m),false);
		directives.Add("effective orbital number",new tw::input::Float(&lstar),false);
	}
	virtual void Initialize() {
		Ionizer::Initialize();
		m = std::fabs(m);
		const tw::Float e = std::exp(1.0);
		const tw::Float sn2l2 = std::sqrt(nstar*nstar-lstar*lstar);
		A1 = std::sqrt(3/cub(pi)) * (2*l+1) * std::tgamma(l+m+1) / std::tgamma(m+1) / std::tgamma(l-m+1);
		A1 *= std::pow(e/sn2l2,m+1.5); // exponent is of poor print quality in JETP
		A1 *= std::pow((nstar+lstar)/(nstar-lstar),lstar+0.5);
		A1 *= (Z*Z/cub(nstar));
		A1 *= std::pow(4*e*cub(Z/nstar)/sn2l2,2*nstar-m-1.5);
		A2 = m+1.5-2*nstar;
		A3 = -2*cub(Z/nstar)/3; // ADK 1986 has nstar**4 in Eq. 21, must be a typo?
		DeduceStaticCoeff();
		I1 = I1*tw::dims::angular_frequency >> atomic >> native;
		A1 = A1*tw::dims::angular_frequency >> atomic >> native;
	}
	virtual void RegisterTests() {
		REGISTER(ADK,HeTest);
	}
	void HeTest();
};

export struct PPT_Tunneling : ADK
{
	PPT_Tunneling(const std::string& name,MetricSpace *m,Task *tsk) : ADK(name,m,tsk) {}
	virtual void Initialize() {
		Ionizer::Initialize();
		m = std::fabs(m);
		const tw::Float F0 = cub(Z/nstar);
		// First without the Coulomb factor
		A1 = Uion;
		A1 *= std::pow(2,2*nstar) / (nstar*std::tgamma(nstar+lstar+1)*std::tgamma(nstar-lstar)); // |C|^2
		A1 *= std::sqrt(6/pi);
		A1 *= (2*l+1) * std::tgamma(l+m+1) / std::pow(2,m) / std::tgamma(m+1) / std::tgamma(l-m+1);
		A1 *= std::pow(0.5/F0,m+1.5);
		A2 = m+1.5;
		A3 = -2*F0/3;
		// Account for Coulomb correction
		A1 *= std::pow(2*F0,2*nstar);
		A2 -= 2*nstar;
		DeduceStaticCoeff();
		I1 = I1*tw::dims::angular_frequency >> atomic >> native;
		A1 = A1*tw::dims::angular_frequency >> atomic >> native;
	}
	virtual void RegisterTests() {
		REGISTER(PPT_Tunneling,HeTest);
	}
	void HeTest();
};

export struct PPT : Ionizer
{
	tw::Int terms;
	PPT(const std::string& name,MetricSpace *m,Task *tsk) : Ionizer(name,m,tsk) {
		terms = 1;
		directives.Add("terms",new tw::input::Int(&terms));
		directives.Add("orbital number",new tw::input::Float(&l),false);
		directives.Add("effective orbital number",new tw::input::Float(&lstar),false);
		// the projection must be zero, do not accept input for it
	}
	virtual void Initialize() {
		Ionizer::Initialize();
		// use A1 to hold C_nl^2
		A1 = std::pow(2,2*nstar);
		A1 /= nstar*std::tgamma(nstar+lstar+1)*std::tgamma(nstar-lstar);
		// use A3 to hold F0
		A3 = std::pow(2*Uion,tw::Float(1.5));
	}
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E) {
		tw::Float ans;
		const tw::Float wa = (w0*tw::dims::angular_frequency >> native >> atomic);
		const tw::Float Ea = (std::fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
		const tw::Float F = Ea/A3;
		const tw::Float gam = std::sqrt(2*Uion)*wa/Ea;
		const tw::Float g = (3/(2*gam))*((1 + 1/(2*gam*gam))*std::asinh(gam) - std::sqrt(1 + gam*gam)/(2*gam));
		const tw::Float nu = (Uion/wa) * (1 + 1/(2*gam*gam));
		ans = Uion*A1*std::sqrt(6/pi)*(2*l+1);
		ans *= std::pow(F*std::sqrt(1 + gam*gam)/2,tw::Float(1.5));
		ans *= (4/std::sqrt(3*pi)) * (gam*gam/(1 + gam*gam)) * FourierSum(gam,nu);
		ans *= std::exp(-2*g/(3*F));
		ans *= std::pow(2/F,2*nstar); // coulomb correction
		return ans*tw::dims::angular_frequency >> atomic >> native;
	}
	tw::Float FourierSum(tw::Float gam,tw::Float nu) {
		tw::Float A0 = 0.0;
		const tw::Float alpha = 2*(std::asinh(gam) - gam/std::sqrt(1 + gam*gam));
		const tw::Float beta = 2*gam/std::sqrt(1 + gam*gam);
		const tw::Float dnu = MyCeil(nu)-nu;
		for (tw::Int n=0;n<terms;n++)
			A0 += std::exp(-alpha*(n+dnu))*tw::dawsoni(std::sqrt(beta*(n+dnu)));
		return A0;
	}
};

export struct PMPB : PPT
{
	PMPB(const std::string& name,MetricSpace *m,Task *tsk) : PPT(name,m,tsk) {}
	virtual void Initialize()
	{
		Ionizer::Initialize();
		// use A1 to hold C_nl^2, n.b. PMPB has a different convention from PPT (factor of 4)
		A1 = std::pow(2,2*nstar-2);
		A1 /= nstar*std::tgamma(nstar+lstar+1)*std::tgamma(nstar-lstar);
		// use A3 to hold F0
		A3 = std::pow(2*Uion,tw::Float(1.5));
	}
	virtual tw::Float AverageRate(tw::Float w0,tw::Float E)
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
		ans *= std::pow(Uion/wa,-1.5);
		ans *= std::sqrt(2*gam/std::sqrt(1+gam*gam)) * FourierSum(gam,nu);
		ans *= std::exp(-2*g/(3*F));
		// Following is the improved Coulomb correction factor
		ans *= std::pow(2/F,2*nstar) * std::pow(1+2*gam/std::exp(1),-2*nstar);
		return ans*tw::dims::angular_frequency >> atomic >> native;
	}
};