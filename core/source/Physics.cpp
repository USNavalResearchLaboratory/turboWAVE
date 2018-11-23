#include "definitions.h"
#include "ctools.h"
#include "3dmath.h"
#include "functions.h"

// ASHER_MOD
//#include "sim.h"
//#include "fluid.h"
#include "tasks.h"
#include "metricSpace.h"
#include "computeTool.h"
#include "3dfields.h"

#include "physics.h"

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

// DFG - encapsulated material parameters can now read their own input file directives
// The CGS units in the input file should really be handled with conversion macros; that would break input files however.
// To keep doing it the usual way we need the unit converter as an argument.
void sparc::material::ReadInputFileDirective(std::stringstream& inputString,const std::string& command,const UnitConverter& uc)
{
	std::string word;
	if (command=="mass") // eg, mass = 1800.0
		inputString >> word >> mass;
	if (command=="charge") // eg, charge = 1.0
		inputString >> word >> charge;
	if (command=="rotational") // eg, rotational degrees of freedom = 2
	{
		inputString >> word >> word >> word >> word >> cvm;
		cvm = 0.5*(3.0 + cvm);
	}
	if (command=="cv") // eg, cv = 2.5 (this is really m*cv, or m*cv/kB)
	{
		inputString >> word >> cvm;
	}
	if (command=="vibrational" || command=="vibration" || command=="excitation") // eg, vibrational energy = 0.3
	{
		inputString >> word >> word >> excitationEnergy;
		excitationEnergy = uc.eV_to_sim(excitationEnergy);
	}
	if (command=="thermometric") // eg, thermometric conductivity = 1.0 // cm^2/s
	{
		inputString >> word >> word >> thermometricConductivity;
		thermometricConductivity *= uc.CGSValue(time_dim)/sqr(uc.CGSValue(length_dim));
	}
	if (command=="kinematic") // eg, kinematic viscosity = 1.0 // cm^2/s
	{
		inputString >> word >> word >> kinematicViscosity;
		kinematicViscosity *= uc.CGSValue(time_dim)/sqr(uc.CGSValue(length_dim));
	}
	if (command=="permittivity") // eg, permittivity = 5.0 , 0.0
		inputString >> word >> eps_r >> eps_i;
}

tw::Float sparc::CoulombCrossSection(const UnitConverter& uc,tw::Float q1,tw::Float q2,tw::Float m12,tw::Float v12,tw::Float N1,tw::Float N2,tw::Float T1,tw::Float T2)
{
	tw::Float rmin,rmax,rmin_alt,coulombLog,hbar;
	hbar = uc.hbar / (uc.MKSValue(time_dim) * uc.MKSValue(energy_dim));
	rmin = fabs(q1*q2)/(tw::small_pos + 4.0*pi*m12*v12*v12*uc.MKSValue(number_dim));
	rmin_alt = hbar/(tw::small_pos + 2.0*m12*v12);
	if (rmin_alt<rmin) rmin = rmin_alt;
	rmax = 1.0/sqrt(tw::small_pos + N1*q1*q1/T1 + N2*q2*q2/T2);
	coulombLog = log(rmax/rmin);
	if (coulombLog<1.0) coulombLog = 1.0;
	return (32.0/pi)*pow(v12,tw::Float(-4.0))*sqr(q1*q2/(4*pi*m12))*coulombLog/uc.MKSValue(number_dim);
}

tw::Float sparc::ElectronPhononRateCoeff(const UnitConverter& uc,tw::Float Ti,tw::Float EFermi,tw::Float ks,tw::Float nref)
{
	// here, the rate coefficient is collision frequency divided by a reference density
	// this allows us to multiply away the collisions in regions where the metallic density is at "background level"
	// get quantities into cgs
	tw::Float vF = uc.c*100.0*sqrt(2.0*EFermi);
	tw::Float qe = uc.SimToCGS(charge_dim,1.0);
	tw::Float kB_Ti = uc.SimToCGS(energy_dim,Ti);
	tw::Float hbar = 1e7 * uc.hbar;
	// collision frequency in real units
	tw::Float nu = 2.0*ks*qe*qe*kB_Ti/(sqr(hbar)*vF);
	// return normalized rate coefficient
	return (uc.CGSValue(time_dim)/nref) * nu;
}

//////////////////////////
//                      //
// Ionization Quasitool //
//                      //
//////////////////////////


IonizationData::IonizationData()
{
	ionizationPotential = 1.0; // normalized to hydrogen
	electrons = 0;
	protons = 0;
	ionizationModel = tw::ionization_model::none;
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

	if (ionizationModel==tw::ionization_model::MPI)
	{
		ans = two*pi*w0*pow(peak/E_MPI,two*photons) / Factorial(photons-1);
		if (ans > max_rate) ans = max_rate;
		return ans*f_atomic_to_sim;
	}
	if (ionizationModel==tw::ionization_model::ADK)
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
	if (ionizationModel==tw::ionization_model::PPT)
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

void IonizationData::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	// read ionspecies and electronspecies indices in Species::ReadInputFileDirective
	// setup hydro indexing during Chemical::Initialize

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
				ionizationModel = tw::ionization_model::none;
			if (word=="adk")
				ionizationModel = tw::ionization_model::ADK;
			if (word=="ppt")
				ionizationModel = tw::ionization_model::PPT;
			if (word=="mpi")
				ionizationModel = tw::ionization_model::MPI;
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
	if (command=="adk")
		inputString >> word >> word >> adkMultiplier;
	if (command=="ppt")
		inputString >> word >> word >> pptMultiplier;
	if (command=="saturated") // eg, saturated rate = 0.01
		inputString >> word >> word >> max_rate;
	if (command=="ion") // eg, ion species = N3
		inputString >> word >> word >> ion_name;
	if (command=="electron") // eg, electron species = electrons
		inputString >> word >> word >> electron_name;
}

void IonizationData::ReadData(std::ifstream& inFile)
{
	inFile.read((char *)&ionizationModel,sizeof(ionizationModel));
	inFile.read((char *)&ionizationPotential,sizeof(ionizationPotential));
	inFile.read((char *)&electrons,sizeof(electrons));
	inFile.read((char *)&protons,sizeof(protons));
	inFile.read((char *)&adkMultiplier,sizeof(adkMultiplier));
	inFile.read((char *)&pptMultiplier,sizeof(pptMultiplier));
	inFile.read((char *)&photons,sizeof(photons));
	inFile.read((char *)&w0,sizeof(w0));
	inFile.read((char *)&E_MPI,sizeof(E_MPI));
	inFile.read((char *)&max_rate,sizeof(max_rate));
	inFile.read((char *)&terms,sizeof(terms));
	inFile.read((char *)&ionSpecies,sizeof(ionSpecies));
	inFile.read((char *)&electronSpecies,sizeof(electronSpecies));
	inFile.read((char *)&hi,sizeof(hi));
	inFile.read((char *)&he,sizeof(he));
	inFile.read((char *)&hgas,sizeof(hgas));
}

void IonizationData::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)&ionizationModel,sizeof(ionizationModel));
	outFile.write((char *)&ionizationPotential,sizeof(ionizationPotential));
	outFile.write((char *)&electrons,sizeof(electrons));
	outFile.write((char *)&protons,sizeof(protons));
	outFile.write((char *)&adkMultiplier,sizeof(adkMultiplier));
	outFile.write((char *)&pptMultiplier,sizeof(pptMultiplier));
	outFile.write((char *)&photons,sizeof(photons));
	outFile.write((char *)&w0,sizeof(w0));
	outFile.write((char *)&E_MPI,sizeof(E_MPI));
	outFile.write((char *)&max_rate,sizeof(max_rate));
	outFile.write((char *)&terms,sizeof(terms));
	outFile.write((char *)&ionSpecies,sizeof(ionSpecies));
	outFile.write((char *)&electronSpecies,sizeof(electronSpecies));
	outFile.write((char *)&hi,sizeof(hi));
	outFile.write((char *)&he,sizeof(he));
	outFile.write((char *)&hgas,sizeof(hgas));
}

// ASHER_MOD

//////////////////////////////////
//                              //
// Equation of State Base Class //
//                              //
//////////////////////////////////


EOSDataTool::EOSDataTool(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	// DFG - must not set names anymore here, let the base class constructor handle it
	// Note enumerated types are now strongly typed and namespaced.
	typeCode = tw::tool_type::eosData;
}

////////////////////////////////
//                            //
// Ideal Gas Law for Chemical //
//                            //
////////////////////////////////

EOSIdealGas::EOSIdealGas(const std::string& name,MetricSpace *m, Task *tsk) : EOSDataTool(name,m,tsk)
{
	typeCode = tw::tool_type::eosIdealGas;
}

void EOSIdealGas::ComponentContribution(Field& hydro, Field& eos)
{
	#pragma omp parallel firstprivate(hidx,eidx,mat)
	{
		for (CellIterator cell(*space,false);cell<cell.end();++cell)
		{
			const tw::Float ngas = hydro(cell,hidx.ni);
			// DFG - temperature is not treated as additive, worked out by parent object (same as previous approach)
			eos(cell,eidx.P) += ngas*eos(cell,eidx.T);
			eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * ngas;
			eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * ngas;
		}
	}
}

void EOSIdealGas::ApplyEOS(tw::Int ie, ScalarField& nu_e,Field& hydro, Field& eos)
{
	#pragma omp parallel firstprivate(ie,hidx,eidx,mat)
	{
		for (CellIterator cell(*space,false);cell<cell.end();++cell)
		{
			const tw::Float ngas = hydro(cell,hidx.ni);
			const tw::Float nmcv = ngas*mat.cvm;
			eos(cell,eidx.T) = hydro(cell,hidx.u)/(tw::small_pos + nmcv);
			eos(cell,eidx.Tv) = 0.0;
			eos(cell,eidx.P) += ngas*eos(cell,eidx.T);
			eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * ngas;
			eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * ngas;
		}
	}
}


/////////////////////////////////
//                             //
//    EOS for hot electrons    //
//                             //
/////////////////////////////////

EOSHotElectrons::EOSHotElectrons(const std::string& name,MetricSpace *m, Task *tsk) : EOSDataTool(name,m,tsk)
{
	typeCode = tw::tool_type::eosHotElectrons;
}

void EOSHotElectrons::ApplyEOS(tw::Int ie, ScalarField& nu_e,Field& hydro, Field& eos)
{
	#pragma omp parallel firstprivate(ie,hidx,eidx,mat)
	{
		for (CellIterator cell(*space,false);cell<cell.end();++cell)
		{
			const tw::Float ne = hydro(cell,ie);
			// EOS for classical electrons
			const tw::Float nmcv = 1.5*ne;
			eos(cell,eidx.T) = hydro(cell,hidx.u)/(tw::small_pos + nmcv);
			eos(cell,eidx.Tv) = 0.0;
			eos(cell,eidx.P) = ne*eos(cell,eidx.T);
			eos(cell,eidx.K) = 3.2*ne*eos(cell,eidx.T)/(mat.mass*nu_e(cell));
			eos(cell,eidx.visc) = 0.0;
			// Braginskii has for e-viscosity 0.73*ne*eos(cell,eidx.T)/nu_e(cell)
			// However, we are forcing electrons to move with ions and so should not diffuse velocity field
		}
	}
}


///////////////////////////////////////
//                                   //
// EOS Mixture for EquilibriumGroup  //
//                                   //
///////////////////////////////////////

EOSMixture::EOSMixture(const std::string& name,MetricSpace *m, Task *tsk) : EOSDataTool(name,m,tsk)
{
	typeCode = tw::tool_type::eosMixture;
	// DFG - it seems we don't need these anymore...
	// eos_tmp.Initialize(5,*space,task);
	// eos_tmp2.Initialize(5,*space,task);
}

void EOSMixture::ApplyEOS(tw::Int ie, ScalarField& nu_e,Field& hydro, Field& eos)
{
	// DFG - Modification from previous version:
	// First set temperature and reset other EOS quantities.  IE is encapsulated as an inline function.
	// Next loop over components and call a distinct component function to add in P,K,visc.

	for (CellIterator cell(*space,false);cell<cell.end();++cell)
	{
		const tw::Float nm = MassDensity(hydro,cell);
		const tw::Float nmcv = MixHeatCapacities(hydro,cell);
		const tw::Float epsvn = MixVibrationalEnergy(hydro,cell);
		const tw::Float nv = MixVibrationalStates(hydro,cell);
		const tw::Float IE = InternalEnergy(nm,hydro,cell);

		eos(cell,eidx.T) = IE/(tw::small_pos + nmcv);
		eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
		eos(cell,eidx.P) = 0.0;
		eos(cell,eidx.K) = 0.0;
		eos(cell,eidx.visc) = 0.0;
		eos(cell,eidx.Cv) = 0.0; // How to handle this...
	}

	// component contribution Loop - load P,K,visc additively
	for (tw::Int s=0;s<eosComponents.size();s++)
			eosComponents[s]->ComponentContribution(hydro,eos);
}

///////////////////////////////////////
//                                   //
//    Ideal Gas Mix - Eliminate?     //
//                                   //
///////////////////////////////////////

EOSIdealGasMix::EOSIdealGasMix(const std::string& name,MetricSpace *m, Task *tsk) : EOSMixture(name,m,tsk)
{
	typeCode = tw::tool_type::eosIdealGasMix;
}

///////////////////////
//                   //
// EOS Mie Gruneisen //
//                   //
///////////////////////

EOSMieGruneisen::EOSMieGruneisen(const std::string& name,MetricSpace *m, Task *tsk) : EOSDataTool(name,m,tsk)
{
	typeCode = tw::tool_type::eosMieGruneisen;
	GRUN = 2.0; // value for Cu on p. 257 of "Shock Wave Physics and Equation of State Modeling"
	// GRUN = 0.1; // value for water in the above book.
}

void EOSMieGruneisen::ComponentContribution(Field& hydro, Field& eos)
{
	#pragma omp parallel firstprivate(hidx,eidx,mat)
	{
		for (CellIterator cell(*space,false);cell<cell.end();++cell)
		{
			const tw::Float nion = hydro(cell,hidx.ni);
			const tw::Float IE = InternalEnergy(nion*mat.mass,hydro,cell);
			// Temperature is not treated as additive, worked out by parent object
			eos(cell,eidx.P) += GRUN*IE;
			eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * nion;
			eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * nion;
		}
	}
}

void EOSMieGruneisen::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	// DFG - this is how the tool gets data from the input file.
	// This supports both named tools and creating on the fly as before.
	std::string word;
	ComputeTool::ReadInputFileDirective(inputString,command);
	// expected form is : gruneisen parameter = 2.0
	if (command=="gruneisen")
		inputString >> word >> word >> GRUN;
}

void EOSMieGruneisen::ReadData(std::ifstream& inFile)
{
	// DFG - now tools can have ReadData/WriteData to support their own restarts
	EOSDataTool::ReadData(inFile);
	inFile.read((char*)&GRUN,sizeof(GRUN));
}

void EOSMieGruneisen::WriteData(std::ofstream& outFile)
{
	EOSDataTool::WriteData(outFile);
	outFile.write((char*)&GRUN,sizeof(GRUN));
}

/////////////////////////
//                     //
// EOS Mie Gruneisen 2 //
//                     //
/////////////////////////
// better MieGruneisen EOS model that uses a linear Hugoniot fit

EOSMieGruneisen2::EOSMieGruneisen2(const std::string& name,MetricSpace *m, Task *tsk) : EOSDataTool(name,m,tsk)
{
	typeCode = tw::tool_type::eosMieGruneisen2;
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
}

void EOSMieGruneisen2::ComponentContribution(Field& hydro, Field& eos)
{
	#pragma omp parallel firstprivate(GRUN,n0,c0,S1,hidx,eidx,mat)
	{
		for (CellIterator cell(*space,false);cell<cell.end();++cell)
		{
			const tw::Float nion = hydro(cell,hidx.ni);
			const tw::Float IE = InternalEnergy(nion*mat.mass,hydro,cell);
			const tw::Float mu = nion/n0 - 1.0;
			const tw::Float sel = tw::Float(mu>=0.0);
			// Temperature is not treated as additive, worked out by parent object
			// Pressure from <http://bluevistasw.com/2016/02/16/mie-gruneisen-eos-implementation/>
			eos(cell,eidx.P) += (1.0-sel)*(mat.mass*n0*c0*c0*mu + GRUN*(mu+1.0)*IE);
			eos(cell,eidx.P) += sel*(mat.mass*n0*c0*c0*mu*(1.0 + (1.0 - GRUN*(mu+1.0)/2.0)*mu)/(1.0 - (S1-1.0)*mu) + GRUN*(mu+1.0)*IE);
			eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * nion;
			eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * nion;
		}
	}
}

void EOSMieGruneisen2::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;
	ComputeTool::ReadInputFileDirective(inputString,command);
	if (command=="gruneisen") // eg, gruneisen parameter = 2.0
		inputString >> word >> word >> GRUN;
	if (command=="reference") // eg, reference density = 1.0
		inputString >> word >> word >> n0;
	if (command=="hugoniot")
	{
		inputString >> word;
		if (word=="intercept")
			inputString >> word >> c0;
		if (word=="slope")
			inputString >> word >> S1;
	}
}

void EOSMieGruneisen2::ReadData(std::ifstream& inFile)
{
	EOSDataTool::ReadData(inFile);
	inFile.read((char*)&GRUN,sizeof(GRUN));
	inFile.read((char*)&n0,sizeof(n0));
	inFile.read((char*)&c0,sizeof(c0));
	inFile.read((char*)&S1,sizeof(S1));
}

void EOSMieGruneisen2::WriteData(std::ofstream& outFile)
{
	EOSDataTool::WriteData(outFile);
	outFile.write((char*)&GRUN,sizeof(GRUN));
	outFile.write((char*)&n0,sizeof(n0));
	outFile.write((char*)&c0,sizeof(c0));
	outFile.write((char*)&S1,sizeof(S1));
}
