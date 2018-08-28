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

// ASHER_MOD

//////////////////////////////////
//                              //
// Equation of State Base Class //
//                              //
//////////////////////////////////


EOSDataTool::EOSDataTool(MetricSpace *m,Task *tsk, bool shared) : ComputeTool(m, tsk,shared)
{
	name = "eos_data";
	typeCode = eosDataTool;
}


EOSDataTool::~EOSDataTool()
{
// In subclasses which use table arrays, they must be destroyed in an overloaded version of this function.
}

void EOSDataTool::Initialize()
{
	ComputeTool::Initialize();
 // When this is overloaded by an appripriate subclass, table arrays may perhaps be initialized here.
}

void EOSDataTool::ApplyEOS(
								std::valarray<tw::Float> mass,std::valarray<tw::Float> charge, std::valarray<tw::Float> cvm,
								std::valarray<tw::Float> excitationEnergy, std::valarray<tw::Float> thermo_cond_cvm,
								std::valarray<tw::Float> k_visc_m, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos)
{
	// include statements here
}

void EOSDataTool::ApplyEOS(
								tw::Float mass,tw::Float charge,tw::Float cvm,
								tw::Float excitationEnergy,tw::Float thermometricConductivity,
								tw::Float kinematicViscosity, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos)
{
	// include statements here
}


////////////////////////////////
//                            //
// Ideal Gas Law for Chemical //
//                            //
////////////////////////////////

EOSIdealGas::EOSIdealGas(MetricSpace *m, Task *tsk,bool shared) : EOSDataTool(m,tsk,shared)
{
	name = "eos_ideal_gas";
	typeCode = eosIdealGas;
}

EOSIdealGas::~EOSIdealGas()
{
}

void EOSIdealGas::Initialize()
{
	EOSDataTool::Initialize();
}


void EOSIdealGas::ApplyEOS(
								tw::Float mass,tw::Float charge,tw::Float cvm,
								tw::Float excitationEnergy,tw::Float thermometricConductivity,
								tw::Float kinematicViscosity, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos)
{
	tw::Float ne;

	// loop over cells
	for (CellIterator cell(*space,false);cell<cell.end();++cell)
	{
		ne = hydro(cell,ie);
		// T and Tv are calculated at the EOSMixture level
		//
		// eos(cell,T) = IE/(tw::small_pos + nmcv); // do temperature in EOSMixture loop?
		// eos(cell,Tv) = (epsvn/(nv+tw::small_pos))/log(1.0001 + epsvn/(hydro(cell,Xi)+tw::small_pos)); // in EOSMix loop
		eos(cell,P) = ne*eos(cell,T);
		eos(cell,K) = thermometricConductivity*ne;//thermometricConductivity*nmcv;
		eos(cell,visc) = kinematicViscosity*ne;//kinematicViscosity*nm;
	}

}


/////////////////////////////////
//                             //
// Ideal Gas Law for electrons //
//                             //
/////////////////////////////////
EOSIdealGasElectrons::EOSIdealGasElectrons(MetricSpace *m, Task *tsk,bool shared) : EOSDataTool(m,tsk,shared)
{
	name = "eos_ideal_gas_electrons";
	typeCode = eosIdealGasElectrons;
}

EOSIdealGasElectrons::~EOSIdealGasElectrons()
{
}

void EOSIdealGasElectrons::Initialize()
{
	EOSDataTool::Initialize();
}


void EOSIdealGasElectrons::ApplyEOS(
								tw::Float mass,tw::Float charge,tw::Float cvm,
								tw::Float excitationEnergy,tw::Float thermometricConductivity,
								tw::Float kinematicViscosity, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos)
{
	tw::Float ne, nmcv;

	// loop over cells
	for (CellIterator cell(*space,false);cell<cell.end();++cell)
	{
		ne = hydro(cell,ie);

		// EOS for electrons
		nmcv = 1.5*ne;
		eos(cell,T) = hydro(cell,U)/(tw::small_pos + nmcv);
		eos(cell,Tv) = 0.0;
		eos(cell,P) = ne*eos(cell,T);
		eos(cell,K) = 3.2*ne*eos(cell,T)/(mass*nu_e(cell));
		eos(cell,visc) = 0.0;
		// Braginskii has for e-viscosity 0.73*ne*eos(cell,electrons->group->T)/nu_e(cell)
		// However, we are forcing electrons to move with ions and so should not diffuse velocity field
	}

}


///////////////////////////////////////
//                                   //
// EOS Mixture for EquilibriumGroup  //
//                                   //
///////////////////////////////////////

EOSMixture::EOSMixture(MetricSpace *m, Task *tsk,bool shared) : EOSDataTool(m,tsk,shared)
{
	name = "eos_mixture";
	typeCode = eosMixture;

	// initialize temporary eos field for averaging calculations
	// I just need enough temporary space to evaluate the eos values for each eos variable (5 of them)
	eos_tmp.Initialize(5,*space,task);
	eos_tmp2.Initialize(5,*space,task);
}

EOSMixture::~EOSMixture()
{
}

void EOSMixture::Initialize()
{
	EOSDataTool::Initialize();
}

void EOSMixture::ApplyEOS(
								std::valarray<tw::Float> mass,std::valarray<tw::Float> charge, std::valarray<tw::Float> cvm,
								std::valarray<tw::Float> excitationEnergy, std::valarray<tw::Float> thermo_cond_cvm,
								std::valarray<tw::Float> k_visc_m, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos)
{
	tw::Float n_ratio, nm_ratio, nmcv_ratio, ntot, nm, nmcv, epsvn, nv, nm_chem, nmcv_chem, IE, KE;
	tw::Int n_idx;
	tw::vec3 chemVelocity;

	// tw::Float SC1, SC2, SC3, SC4, SC5; // for debugging purposes

	// loop over cells and calculate temperature
	// only the subsequent EOSs are calculated and updated accordingly
	// an EqiulibriumGroup is defined to have a single temperature, after all
	for (CellIterator cell(*space,false);cell<cell.end();++cell)
		{
			// total density at this cell
			ntot = this->DensitySum(e,hydro,cell);
			// total mass density at this cell
			nm = this->DensityWeightedSum(e,hydro,mass,cell);
			nmcv = this->DensityWeightedSum(e,hydro,cvm,cell);
			epsvn = this->DensityWeightedSum(e,hydro,excitationEnergy,cell);
			nv = this->ConditionalDensitySum(e,hydro,excitationEnergy,cell);

			//chemVelocity = vec3(hydro(cell,npx),hydro(cell,npy),hydro(cell,npz))/(tw::small_pos + nm);
			KE = 0.5*nm*Norm(this->Velocity(e,npx,npy,npz,hydro,mass,cell));
			IE = hydro(cell,U) - KE;
			if (IE<=0.0)
			{
				hydro(cell,U) = KE*1.00001 + tw::small_pos;
				IE = hydro(cell,U) - KE;
			}

			// calculate temperature
			eos_tmp(cell,0) = IE/(tw::small_pos + nmcv);
			eos_tmp(cell,1) = IE; // pass on internal energy to eosComponents

			eos_tmp2(cell,0) = IE/(tw::small_pos + nmcv);
			// I imagine the vibrational temperature is also at equilibrium by definition
			eos_tmp2(cell,1) = (epsvn/(nv+tw::small_pos))/log(1.0001 + epsvn/(hydro(cell,Xi)+tw::small_pos));

	}


	// loop over elements (Chemical objects) of the group
	n_idx = e.low;
	for (tw::Int s=0;s<eosComponents.size();s++) {

		// for (CellIterator cell(*space,false);cell<cell.end();++cell)
		// {
		// 		eos_tmp(cell,0) = 0.0;
		// 		eos_tmp(cell,1) = 0.0;
		// 	 	eos_tmp(cell,2) = 0.0;
		// 	 	eos_tmp(cell,3) = 0.0;
		// 	 	eos_tmp(cell,4) = 0.0;
		// }

		// apply component EOS, deposit into eos_tml Field
		eosComponents[s]->ApplyEOS( mass[s],charge[s],cvm[s],
									 excitationEnergy[s],thermo_cond_cvm[s],k_visc_m[s],
									 e,0,1,2,3,
									 4,5,npx,npy,npz,
									 U,Xi,n_idx,nu_e,hydro,eos_tmp);

		// loop over cells
		for (CellIterator cell(*space,false);cell<cell.end();++cell)
		{
			// if first element then wipe eos_tmp2
			if (s == 0) {
//				eos_tmp2(cell,0) = 0.0;
//				eos_tmp2(cell,1) = 0.0;
			 	eos_tmp2(cell,2) = 0.0;
			 	eos_tmp2(cell,3) = 0.0;
			 	eos_tmp2(cell,4) = 0.0;
			}

			// SC1 = eos_tmp(cell,0);
			// SC2 = eos_tmp(cell,1);
			// SC3 = eos_tmp(cell,2);
			// SC4 = eos_tmp(cell,3);
			// SC5 = eos_tmp(cell,4);

			// if (std::isnan(SC1)||std::isnan(SC2)||std::isnan(SC3)||std::isnan(SC4)||std::isnan(SC5)) {
			// 	std::cout << "eos_tmp Loop\n";
			// 	std::cout << "SC1: " << SC1 << "\n";
			// 	std::cout << "SC2: " << SC2 << "\n";
			// 	std::cout << "SC3: " << SC3 << "\n";
			// 	std::cout << "SC4: " << SC4 << "\n";
			// 	std::cout << "SC5: " << SC5 << "\n";
			// exit(0);
			// }

			eos_tmp2(cell,2) += eos_tmp(cell,2); // P
			eos_tmp2(cell,3) += eos_tmp(cell,3); // K
			eos_tmp2(cell,4) += eos_tmp(cell,4); // visc

		}

		// increment density index
		n_idx++;
	}


	// loop over cells -- copy over the density weighted values in eos_tmp2 to the output eos field
	for (CellIterator cell(*space,false);cell<cell.end();++cell)
	{
		// SC1 = eos_tmp2(cell,0);
		// SC2 = eos_tmp2(cell,1);
		// SC3 = eos_tmp2(cell,2);
		// SC4 = eos_tmp2(cell,3);
		// SC5 = eos_tmp2(cell,4);

		// if (std::isnan(SC1)||std::isnan(SC2)||std::isnan(SC3)||std::isnan(SC4)||std::isnan(SC5)) {
		// 	std::cout << "Final EOS Copy Loop\n";
		// 	std::cout << "SC1: " << SC1 << "\n";
		// 	std::cout << "SC2: " << SC2 << "\n";
		// 	std::cout << "SC3: " << SC3 << "\n";
		// 	std::cout << "SC4: " << SC4 << "\n";
		// 	std::cout << "SC5: " << SC5 << "\n";
		// 	exit(0);
		// }

		eos(cell,T)    = eos_tmp2(cell,0);
		eos(cell,Tv)   = eos_tmp2(cell,1);
		eos(cell,P)    = eos_tmp2(cell,2);
		eos(cell,K)    = eos_tmp2(cell,3);
		eos(cell,visc) = eos_tmp2(cell,4);
	}

}


/////////////////////////////////////////////////
//                                             //
// EOS Ideal Gas Mixture for EquilibriumGroup  //
//                                             //
/////////////////////////////////////////////////

EOSIdealGasMix::EOSIdealGasMix(MetricSpace *m, Task *tsk,bool shared) : EOSMixture(m,tsk,shared)
{
	name = "eos_ideal_gas_mix";
	typeCode = eosIdealGasMix;
}

EOSIdealGasMix::~EOSIdealGasMix()
{
}

void EOSIdealGasMix::Initialize()
{
	EOSDataTool::Initialize();
}

void EOSIdealGasMix::ApplyEOS(
								std::valarray<tw::Float> mass,std::valarray<tw::Float> charge, std::valarray<tw::Float> cvm,
								std::valarray<tw::Float> excitationEnergy, std::valarray<tw::Float> thermo_cond_cvm,
								std::valarray<tw::Float> k_visc_m, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos)
{
	tw::Float ntot,ne,ionChargeDensity,nm,nmcv,epsvn,nv,KE,IE;

	// include statements here

	// loop over cells
	for (CellIterator cell(*space,false);cell<cell.end();++cell)
	{
		ntot = this->DensitySum(e,hydro,cell);
		nm = this->DensityWeightedSum(e,hydro,mass,cell);
		nmcv = this->DensityWeightedSum(e,hydro,cvm,cell);
		epsvn = this->DensityWeightedSum(e,hydro,excitationEnergy,cell);
		nv = this->ConditionalDensitySum(e,hydro,excitationEnergy,cell);
		KE = 0.5*nm*Norm(this->Velocity(e,npx,npy,npz,hydro,mass,cell));
		IE = hydro(cell,U) - KE;
		if (IE<=0.0)
		{
			hydro(cell,U) = KE*1.00001 + tw::small_pos;
			IE = hydro(cell,U) - KE;
		}

		eos(cell,T) = IE/(tw::small_pos + nmcv);
		eos(cell,Tv) = (epsvn/(nv+tw::small_pos))/log(1.0001 + epsvn/(hydro(cell,Xi)+tw::small_pos));
		eos(cell,P) = ntot * eos(cell,T);
		eos(cell,K) = this->DensityWeightedSum(e,hydro,thermo_cond_cvm,cell); // thermometricConductivity * nmcv;
		eos(cell,visc) = this->DensityWeightedSum(e,hydro,k_visc_m,cell); // kinematicViscosity * nm;
	}


}

///////////////////////
//                   //
// EOS Mie Gruneisen //
//                   //
///////////////////////

EOSMieGruneisen::EOSMieGruneisen(MetricSpace *m, Task *tsk,bool shared) : EOSDataTool(m,tsk,shared)
{
	name = "eos_mie_gruneisen";
	typeCode = eosMieGruneisen;
}

EOSMieGruneisen::~EOSMieGruneisen()
{
}

void EOSMieGruneisen::Initialize()
{
	EOSDataTool::Initialize();
}


void EOSMieGruneisen::ApplyEOS(
								tw::Float mass,tw::Float charge,tw::Float cvm,
								tw::Float excitationEnergy,tw::Float thermometricConductivity,
								tw::Float kinematicViscosity, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos)
{
	tw::Float ne, n0; //, nm, nmcv, nv, epsvn;

	tw::Float IE;

	// GRUN = 2.0; // value for Cu on p. 257 of "Shock Wave Physics and Equation of State Modeling"
	// GRUN = 0.1; // value for water in the above book.

	// loop over cells
	for (CellIterator cell(*space,false);cell<cell.end();++cell)
	{
		ne = hydro(cell,ie);
		IE = eos(cell,Tv);

		eos(cell,P) = GRUN*eos(cell,Tv); // 'eos(cell,Tv)' is actuall IE at this point
		eos(cell,K) = thermometricConductivity*ne;//thermometricConductivity*nmcv;
		eos(cell,visc) = kinematicViscosity*ne;//kinematicViscosity*nm;

	}

}


// a different, better MieGruneisen EOS model that uses a linear Hugoniot fit
/////////////////////////
//                     //
// EOS Mie Gruneisen 2 //
//                     //
/////////////////////////

EOSMieGruneisen2::EOSMieGruneisen2(MetricSpace *m, Task *tsk,bool shared) : EOSDataTool(m,tsk,shared)
{
	name = "eos_mie_gruneisen2";
	typeCode = eosMieGruneisen2;
}

EOSMieGruneisen2::~EOSMieGruneisen2()
{
}


void EOSMieGruneisen2::Initialize()
{
	EOSDataTool::Initialize();
}


void EOSMieGruneisen2::ApplyEOS(
								tw::Float mass,tw::Float charge,tw::Float cvm,
								tw::Float excitationEnergy,tw::Float thermometricConductivity,
								tw::Float kinematicViscosity, Element& e,
								tw::Int T,tw::Int Tv,tw::Int P,tw::Int K,tw::Int visc,
								tw::Int Cv,
								tw::Int npx,tw::Int npy,tw::Int npz,tw::Int U,tw::Int Xi,
								tw::Int ie, ScalarField& nu_e,
								Field& hydro, Field& eos)
{
	tw::Float ne; //, nm, nmcv, nv, epsvn;
	tw::Float mu;

	tw::Float IE, IE0, IE1, P0, P1, DP;

	// Hugoniot data fit for Cu
	// n0 = 3.3e3;
	// c0 = 1.3248e-5;
	// S1 = 1.5;

	// Hugoniot data fit for H20
	// n0 = 1334.0;
	// c0 = 5.197e-6;
	// S1 = 1.8153;

	// GRUN = 2.0; // value for Cu on p. 257 of "Shock Wave Physics and Equation of State Modeling"
	// GRUN = 0.1; // value for water in the above book.

	// loop over cells
	for (CellIterator cell(*space,false);cell<cell.end();++cell)
	{
		ne = hydro(cell,ie);
		IE = eos(cell,Tv); // 'eos(cell,Tv)' is actuall IE at this point

		mu = ne/n0 - 1.0;

		// from <http://bluevistasw.com/2016/02/16/mie-gruneisen-eos-implementation/>
		if (mu >= 0){
			eos(cell,P) = mass*n0*c0*c0*mu*(1.0 + (1.0 - GRUN*(mu+1.0)/2.0)*mu)/(1.0 - (S1-1.0)*mu) + GRUN*(mu+1.0)*IE;
		} else {
			eos(cell,P) = mass*n0*c0*c0*mu + GRUN*(mu+1.0)*IE;
		}

		eos(cell,K) = thermometricConductivity*ne;//thermometricConductivity*nmcv;
		eos(cell,visc) = kinematicViscosity*ne;//kinematicViscosity*nm;

	}

}
