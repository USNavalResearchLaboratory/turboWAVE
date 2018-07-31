#include "definitions.h"
#include "ctools.h"
#include "3dmath.h"
#include "functions.h"
#include "physics.h"

/////////////////////
//                 //
// IONIZATION DATA //
//                 //
/////////////////////


IonizationData::IonizationData()
{
	ionizationPotential = 1.0; // normalized to hydrogen
	electrons = 0;
	protons = 0;
	ionizationModel = noIonization;
	adkMultiplier = 1.0;
	pptMultiplier = 1.0;
	C_ADK = C_ADK_AVG = C_PPT = nstar = 0.0;
	terms = 3;
	t_atomic = 1.0/4.13404e16;
	E_atomic = 5.14203e11; // V/m
	photons = 1;
	w0 = 1.0;
	E_MPI = 1.0;
	max_rate = tw::big_pos;
	ionizeFromGas = false;
	ionSpecies = 0;
	electronSpecies = 0;
}

void IonizationData::Initialize(tw::Float unitDensity,tw::Float* carrierFrequency)
{
	tw::vec3 pos;
	UnitConverter uc(unitDensity);

	f_atomic_to_sim = uc.MKSValue(time_dim)/t_atomic;
	E_sim_to_atomic = uc.MKSValue(electric_field_dim)/E_atomic;
	tw::Float Z = protons - electrons + 1.0;
	tw::Float Uion = 0.5*ionizationPotential; // put in a.u. for consistency of notation
	nstar = Z / sqrt(2.0*Uion);
	lstar = nstar - 1.0;

	// Number of photons for simple MPI model

	if (carrierFrequency)
	{
		w0 = (*carrierFrequency) / f_atomic_to_sim; // laser freq. in a.u.
		photons = MyFloor(Uion/w0 + 1.0);
	}

	// ADK tunneling constants

	C_ADK = pow(two*two*exp(one)*pow(Z,tw::Float(3.0))/pow(nstar,tw::Float(4.0)),two*nstar)/(tw::Float(8.0)*pi*Z);
	C_ADK *= adkMultiplier;

	C_ADK_AVG = C_ADK * sqrt(3.0*cub(nstar)/(pi*cub(Z)));

	// PPT ionization constants

	tw::Float C2 = pow(two,two*nstar) / (nstar*tgamma(nstar+lstar+one)*tgamma(nstar-lstar));
	C_PPT = Uion * C2 * sqrt(6.0/pi) * pptMultiplier;
}

tw::Float IonizationData::wfunc(tw::Float x)
{
	// the argument of this function varies from 0 up to sqrt(2*terms)
	// where terms is the number of terms kept in the ppt expansion
	return 0.5 * exp(-x*x) * erfi(x) * sqrt(pi);
}

tw::Float IonizationData::Rate(tw::Float instant,tw::Float peak)
{
	tw::Int i;
	tw::Float ans,C_EXP;
	tw::Float Uion = 0.5*ionizationPotential; // put in a.u.
	instant = fabs(instant)*E_sim_to_atomic;
	peak = fabs(peak)*E_sim_to_atomic;

	if (ionizationModel==mpiSimple)
	{
		ans = two*pi*w0*pow(peak/E_MPI,two*photons) / Factorial(photons-1);
		if (ans > max_rate) ans = max_rate;
		return ans*f_atomic_to_sim;
	}
	if (ionizationModel==ADKTunneling)
	{
		ans = 0.0;
		C_EXP = two*pow(two*Uion,one+half)/tw::Float(3.0);
		instant += 0.01*C_EXP; // avoid underflow and divide by zero
		peak += 0.01*C_EXP;
		ans += C_ADK*pow(instant,one-two*nstar)*exp(-C_EXP/instant);
		ans += C_ADK_AVG*pow(peak,one+half-two*nstar)*exp(-C_EXP/peak);
		if (ans > max_rate) ans = max_rate;
		return ans*f_atomic_to_sim;
	}
	if (ionizationModel==pptIonization)
	{
		ans = 0.0;
		if (peak>tw::small_pos)
		{
			tw::Float F0 = pow(two*Uion,one+half);
			C_EXP = two*F0/tw::Float(3.0);
			peak += 0.01*C_EXP; // avoid underflow and divide by zero
			tw::Float gam = sqrt(2*Uion)*w0/peak;
			tw::Float gam2 = gam*gam;
			tw::Float alpha = 2.0*(asinh(gam) - gam/sqrt(1.0 + gam2));
			tw::Float beta = 2.0*gam/sqrt(1.0 + gam2);
			tw::Float g = (3.0/(2.0*gam))*((1.0 + 1.0/(2.0*gam2))*asinh(gam) - sqrt(1.0 + gam2)/(2.0*gam));
			tw::Float nu = (Uion/w0)*(1.0 + 1.0/(2.0*gam2));
			for (i=tw::Int(MyCeil(nu));i<tw::Int(MyCeil(nu))+terms;i++)
				ans += exp(-alpha*(tw::Float(i)-nu))*wfunc(sqrt(beta*(tw::Float(i)-nu)));
			ans *= (4.0/sqrt(3.0*pi)) * (gam2/(1.0 + gam2));
			ans *= pow(peak*sqrt(one + gam2)/(two*F0),one+half);
			ans *= pow(two*F0/peak,two*nstar); // coulomb correction
			ans *= exp(-C_EXP*g/peak);
			ans *= C_PPT;
		}
		if (ans > max_rate) ans = max_rate;
		return ans*f_atomic_to_sim;
	}
	return 0.0;
}

void IonizationData::ReadInputFileTerm(std::stringstream& inputString,std::string& command)
{
	// read ionspecies and electronspecies indices in Species method

	std::string word;

	if (command=="ionization") // eg, ionization potential = 1.0
	{
		inputString >> word;
		if (word=="potential")
			inputString >> word >> ionizationPotential;
		else
		{
			inputString >> word >> word;
			if (word=="none")
				ionizationModel = noIonization;
			if (word=="adk")
				ionizationModel = ADKTunneling;
			if (word=="ppt")
				ionizationModel = pptIonization;
			if (word=="mpi")
				ionizationModel = mpiSimple;
		}
	}
	if (command=="mpi") // eg, mpi reference field = 1.0
		inputString >> word >> word >> word >> E_MPI;
	if (command=="terms")
		inputString >> word >> terms;
	if (command=="protons")
		inputString >> word >> protons;
	if (command=="electrons")
		inputString >> word >> electrons;
	if (command=="ionize") // eg, ionize from gas = yes
	{
		inputString >> word >> word >> word >> word;
		ionizeFromGas = (word=="yes" || word=="true" || word=="on");
	}
	if (command=="adk")
		inputString >> word >> word >> adkMultiplier;
	if (command=="ppt")
		inputString >> word >> word >> pptMultiplier;
	if (command=="saturated") // eg, saturated rate = 0.01
		inputString >> word >> word >> max_rate;
}

// Normalization Functions

UnitConverter::UnitConverter(tw::Float unitDensityCGS)
{
	n1 = unitDensityCGS*1e6;
	wp = sqrt(n1*sqr(mks::qe)/(mks::eps0*mks::me));
	c = mks::c;
	qe = mks::qe;
	me = mks::me;
	eps0 = mks::eps0;
	kB = mks::kB;
	hbar = mks::hbar;
}

tw::Float UnitConverter::MKSValue(tw_dimensions dim) const
{
	switch (dim)
	{
		case time_dim:
			return 1.0/wp;
		case length_dim:
			return c/wp;
		case number_dim:
			return n1*cub(c/wp);
		case mass_dim:
			return me;
		case energy_dim:
			return me*c*c;
		case momentum_dim:
			return me*c;
		case angular_momentum_dim:
			return me*c*c/wp;
		case density_dim:
			return n1;
		case power_dim:
			return me*c*c*wp;
		case fluence_dim:
			return sqr(me*c*wp/qe)*c*eps0/wp; // E^2/eta0/wp
		case intensity_dim:
			return sqr(me*c*wp/qe)*c*eps0; // E^2/eta0
		case energy_density_dim:
			return me*c*c*n1;
		case power_density_dim:
			return me*c*c*n1*wp;
		case charge_dim:
			return qe;
		case current_dim:
			return n1*qe*c*sqr(c/wp);
		case current_density_dim:
			return n1*qe*c;
		case charge_density_dim:
			return n1*qe;
		case electric_field_dim:
			return me*c*wp/qe;
		case magnetic_field_dim:
			return me*wp/qe;
		case scalar_potential_dim:
			return me*c*c/qe;
		case conductivity_dim:
			return (n1*qe*c)/(me*c*wp/qe); // j/E
		case rate_coefficient_2_dim:
			return wp/n1;
		case rate_coefficient_3_dim:
			return wp/(n1*n1);
		case mobility_dim:
			return c/(me*c*wp/qe); // c/E
		case temperature_dim:
			return me*c*c/kB;
		case cross_section_dim:
			return wp/(n1*c);
		case susceptibility_dim:
			return 1.0;
		case susceptibility_2_dim:
			return 1.0/(me*c*wp/qe); // 1/E
		case susceptibility_3_dim:
			return 1.0/sqr(me*c*wp/qe); // 1/E^2
		default:
			return 0.0;
	}
}

tw::Float UnitConverter::CGSValue(tw_dimensions dim) const
{
	switch (dim)
	{
		case time_dim:
			return 1.0/wp;
		case length_dim:
			return 1e2*c/wp;
		case number_dim:
			return n1*cub(c/wp);
		case mass_dim:
			return 1e3*me;
		case energy_dim:
			return 1e7*me*c*c;
		case momentum_dim:
			return 1e5*me*c;
		case angular_momentum_dim:
			return 1e7*me*c*c/wp;
		case density_dim:
			return 1e-6*n1;
		case power_dim:
			return 1e7*me*c*c*wp;
		case fluence_dim:
			return 1e3*sqr(me*c*wp/qe)*c*eps0/wp; // E^2/eta0/wp
		case intensity_dim:
			return 1e3*sqr(me*c*wp/qe)*c*eps0; // E^2/eta0
		case energy_density_dim:
			return 1e1*me*c*c*n1;
		case power_density_dim:
			return 1e1*me*c*c*n1*wp;
		case charge_dim:
			return 3e9*qe;
		case current_dim:
			return 3e9*n1*qe*c*sqr(c/wp);
		case current_density_dim:
			return 3e5*n1*qe*c;
		case charge_density_dim:
			return 3e3*n1*qe;
		case electric_field_dim:
			return 0.333333e-4*me*c*wp/qe;
		case magnetic_field_dim:
			return 1e4*me*wp/qe;
		case scalar_potential_dim:
			return 0.333333e-2*me*c*c/qe;
		case conductivity_dim:
			return 9e9*(n1*qe*c)/(me*c*wp/qe); // j/E
		case rate_coefficient_2_dim:
			return 1e6*wp/n1;
		case rate_coefficient_3_dim:
			return 1e12*wp/(n1*n1);
		case mobility_dim:
			return 1e2*c/(0.333333e-4*me*c*wp/qe); // c/E
		case temperature_dim:
			return me*c*c/kB;
		case cross_section_dim:
			return 1e4*wp/(n1*c);
		case susceptibility_dim:
			return 1.0/(4.0*pi);
		case susceptibility_2_dim:
			return 1.0/(4.0*pi*0.333333e-4*me*c*wp/qe); // 1/4*pi*E
		case susceptibility_3_dim:
			return 1.0/sqr(4.0*pi*0.333333e-4*me*c*wp/qe); // 1/4*pi*E^2
		default:
			return 0.0;
	}
}
