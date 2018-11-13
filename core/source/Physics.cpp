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

void IonizationData::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
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

// ASHER_MOD

//////////////////////////////////
//                              //
// Equation of State Base Class //
//                              //
//////////////////////////////////


EOSDataTool::EOSDataTool(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
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
			const tw::Float ngas = hydro(cell,hidx.n);
			// Temperature is not treated as additive, worked out by parent object (same as previous approach)
			eos(cell,eidx.P) += ngas*eos(cell,eidx.T);
			eos(cell,eidx.K) += mat.thermo_cond_cvm*ngas;//thermometricConductivity*nmcv;
			eos(cell,eidx.visc) += mat.k_visc_m*ngas;//kinematicViscosity*nm;
		}
	}
}

void EOSIdealGas::ApplyEOS(tw::Int ie, ScalarField& nu_e,Field& hydro, Field& eos)
{
	#pragma omp parallel firstprivate(ie,hidx,eidx,mat)
	{
		for (CellIterator cell(*space,false);cell<cell.end();++cell)
		{
			const tw::Float ngas = hydro(cell,hidx.n);
			const tw::Float nmcv = ngas*mat.cvm;
			eos(cell,eidx.T) = hydro(cell,hidx.u)/(tw::small_pos + nmcv);
			eos(cell,eidx.Tv) = 0.0;
			eos(cell,eidx.P) += ngas*eos(cell,eidx.T);
			eos(cell,eidx.K) += mat.thermo_cond_cvm*ngas;//thermometricConductivity*nmcv;
			eos(cell,eidx.visc) += mat.k_visc_m*ngas;//kinematicViscosity*nm;
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
	// it seems we don't need these anymore...
	// eos_tmp.Initialize(5,*space,task);
	// eos_tmp2.Initialize(5,*space,task);
}

void EOSMixture::ApplyEOS(tw::Int ie, ScalarField& nu_e,Field& hydro, Field& eos)
{
	// Modification from previous version:
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
			const tw::Float n = hydro(cell,hidx.n);
			const tw::Float IE = InternalEnergy(n*mat.mass,hydro,cell);
			// Temperature is not treated as additive, worked out by parent object
			eos(cell,eidx.P) += GRUN*IE;
			eos(cell,eidx.K) += mat.thermo_cond_cvm*n;//thermometricConductivity*nmcv;
			eos(cell,eidx.visc) += mat.k_visc_m*n;//kinematicViscosity*nm;
		}
	}
}

void EOSMieGruneisen::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;
	ComputeTool::ReadInputFileDirective(inputString,command);
	// expected form is : gruneisen parameter = 2.0
	if (command=="gruneisen")
		inputString >> word >> word >> GRUN;
}

void EOSMieGruneisen::ReadData(std::ifstream& inFile)
{
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
			const tw::Float n = hydro(cell,hidx.n);
			const tw::Float IE = InternalEnergy(n*mat.mass,hydro,cell);
			const tw::Float mu = n/n0 - 1.0;
			const tw::Float sel = tw::Float(mu>=0.0);
			// Temperature is not treated as additive, worked out by parent object
			// Pressure from <http://bluevistasw.com/2016/02/16/mie-gruneisen-eos-implementation/>
			eos(cell,eidx.P) += (1.0-sel)*(mat.mass*n0*c0*c0*mu + GRUN*(mu+1.0)*IE);
			eos(cell,eidx.P) += sel*(mat.mass*n0*c0*c0*mu*(1.0 + (1.0 - GRUN*(mu+1.0)/2.0)*mu)/(1.0 - (S1-1.0)*mu) + GRUN*(mu+1.0)*IE);
			eos(cell,eidx.K) += mat.thermo_cond_cvm*n;//thermometricConductivity*nmcv;
			eos(cell,eidx.visc) += mat.k_visc_m*n;//kinematicViscosity*nm;
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
