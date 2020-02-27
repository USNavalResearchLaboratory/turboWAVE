#include "meta_base.h"
#include "computeTool.h"
#include "physics.h"

// DFG - Populating the input hash is partly delegated to this low level material object.
// We are now requiring the user to specify units, otherwise normalized units are assumed.
// *This breaks old input files*, but is necessary for the sake of consistency.
void sparc::material::AddDirectives(tw::input::DirectiveReader& directives)
{
	directives.Add("mass",new tw::input::Float(&mass));
	directives.Add("charge",new tw::input::Float(&charge));
	directives.Add("cv",new tw::input::Float(&cvm));
	directives.Add("vibrational energy",new tw::input::Float(&excitationEnergy));
	directives.Add("thermometric conductivity",new tw::input::Float(&thermometricConductivity));
	directives.Add("kinematic viscosity",new tw::input::Float(&kinematicViscosity));
	directives.Add("permittivity",new tw::input::Numbers<tw::Float>(&eps[0],2));
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
	return 0.5 * exp(-x*x) * tw::erfi(x) * sqrt(pi);
}

tw::Float IonizationData::Echar()
{
	// return characteristic field in simulation units
	const tw::Float Uion = 0.5*ionizationPotential; // put in a.u.
	if (ionizationModel==tw::ionization_model::MPI)
	{
		return E_MPI/E_sim_to_atomic;
	}
	if (ionizationModel==tw::ionization_model::ADK)
	{
		return two*pow(two*Uion,one+half)/tw::Float(3.0)/E_sim_to_atomic;
	}
	if (ionizationModel==tw::ionization_model::PPT)
	{
		const tw::Float F0 = pow(two*Uion,one+half);
		return two*F0/tw::Float(3.0)/E_sim_to_atomic;
	}
	return 0.0;
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

void IonizationData::AddDirectives(tw::input::DirectiveReader& directives)
{
	// read ionspecies and electronspecies indices in Species::Initialize
	// setup hydro indexing during Chemical::Initialize
	directives.Add("ionization potential",new tw::input::Float(&ionizationPotential));
	std::map<std::string,tw::ionization_model> imod = {{"none",tw::ionization_model::none},
		{"adk",tw::ionization_model::ADK},{"ppt",tw::ionization_model::PPT},{"mpi",tw::ionization_model::MPI}};
	directives.Add("ionization model",new tw::input::Enums<tw::ionization_model>(imod,&ionizationModel));
	directives.Add("mpi reference field",new tw::input::Float(&E_MPI));
	directives.Add("terms",new tw::input::Int(&terms));
	directives.Add("protons",new tw::input::Float(&protons));
	directives.Add("electrons",new tw::input::Float(&electrons));
	directives.Add("adk multiplier",new tw::input::Float(&adkMultiplier));
	directives.Add("ppt multiplier",new tw::input::Float(&pptMultiplier));
	directives.Add("saturated rate",new tw::input::Float(&max_rate));
	directives.Add("ion species",new tw::input::String(&ion_name));
	directives.Add("electron species",new tw::input::String(&electron_name));
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


EOSComponent::EOSComponent(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	// DFG - must not set names anymore here, let the base class constructor handle it
	// Note enumerated types are now strongly typed and namespaced.
	typeCode = tw::tool_type::eosData;
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
	typeCode = tw::tool_type::eosIdealGas;
}

/////////////////////////////////
//                             //
//    EOS for hot electrons    //
//                             //
/////////////////////////////////

EOSHotElectrons::EOSHotElectrons(const std::string& name,MetricSpace *m, Task *tsk) : EOSComponent(name,m,tsk)
{
	typeCode = tw::tool_type::eosHotElectrons;
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
	typeCode = tw::tool_type::eosSimpleMieGruneisen;
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

void EOSSimpleMieGruneisen::ReadData(std::ifstream& inFile)
{
	// DFG - now tools can have ReadData/WriteData to support their own restarts
	EOSComponent::ReadData(inFile);
	inFile.read((char*)&GRUN,sizeof(GRUN));
}

void EOSSimpleMieGruneisen::WriteData(std::ofstream& outFile)
{
	EOSComponent::WriteData(outFile);
	outFile.write((char*)&GRUN,sizeof(GRUN));
}

/////////////////////////
//                     //
//  EOS Mie Gruneisen  //
//                     //
/////////////////////////
// better MieGruneisen EOS model that uses a linear Hugoniot fit

EOSLinearMieGruneisen::EOSLinearMieGruneisen(const std::string& name,MetricSpace *m, Task *tsk) : EOSComponent(name,m,tsk)
{
	typeCode = tw::tool_type::eosLinearMieGruneisen;
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
			const tw::Float mu = nion/n0 - 1.0;
			const tw::Float sel = tw::Float(mu>=0.0);
			// Temperature is not treated as additive, worked out by parent object
			// Pressure from <http://bluevistasw.com/2016/02/16/mie-gruneisen-eos-implementation/>
			eos(cell,eidx.P) += (1.0-sel)*(mat.mass*n0*c0*c0*mu + GRUN*(mu+1.0)*partial_IE);
			eos(cell,eidx.P) += sel*(mat.mass*n0*c0*c0*mu*(1.0 + (1.0 - GRUN*(mu+1.0)/2.0)*mu)/(1.0 - (S1-1.0)*mu) + GRUN*(mu+1.0)*partial_IE);
			eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * nion;
			eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * nion;
		}
	}
}

void EOSLinearMieGruneisen::ReadData(std::ifstream& inFile)
{
	EOSComponent::ReadData(inFile);
	inFile.read((char*)&GRUN,sizeof(GRUN));
	inFile.read((char*)&n0,sizeof(n0));
	inFile.read((char*)&c0,sizeof(c0));
	inFile.read((char*)&S1,sizeof(S1));
}

void EOSLinearMieGruneisen::WriteData(std::ofstream& outFile)
{
	EOSComponent::WriteData(outFile);
	outFile.write((char*)&GRUN,sizeof(GRUN));
	outFile.write((char*)&n0,sizeof(n0));
	outFile.write((char*)&c0,sizeof(c0));
	outFile.write((char*)&S1,sizeof(S1));
}

///////////////////////////////////////
//                                   //
// EOS Mixture for EquilibriumGroup  //
//                                   //
///////////////////////////////////////

EOSMixture::EOSMixture(const std::string& name,MetricSpace *m, Task *tsk) : ComputeTool(name,m,tsk)
{
	typeCode = tw::tool_type::eosMixture;
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
			eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
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
			eos(cell,eidx.T) = eosRef(cell,eidx.T) + (IE1*(1.0+nm0/nm1) - IE0*(1.0+nm1/nm0))/nmcv_sum;
			eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
			nm(cell) = nm1;
			IE(cell) = IE1;
		}
	}
}

void EOSMixture::UpdateEnergy(ScalarField& nm,ScalarField& T0,Field& hydro,Field& eos)
{
	// Add energy corresponding to a change in temperature only.
	// Not centered, because cv is not updated.
	// Must not assume hydro field has complete data, because of Chemical::GenerateFluid
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
	typeCode = tw::tool_type::eosIdealGasMix;
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
			eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
			nm(cell) = nm1;
			IE(cell) = IE1;
		}
	}
}

void EOSIdealGasMix::UpdateEnergy(ScalarField& nm,ScalarField& T0,Field& hydro,Field& eos)
{
	// DFG - polytropic ideal gas caloric EOS in the original SPARC mode of calculation.  Ignores reference states.
	// CANNOT USE: incompatible with generalized Chemical::GenerateFluid
	// Therefore call the inherited function and comment out the rest.
	EOSMixture::UpdateEnergy(nm,T0,hydro,eos);
	// #pragma omp parallel
	// {
	// 	for (auto cell : EntireCellRange(*space))
	// 	{
	// 		hydro(cell,hidx.u) = eos(cell,eidx.nmcv) * eos(cell,eidx.T);
	// 		hydro(cell,hidx.u) += 0.5*sqr(hydro(cell,hidx.npx))/(tw::small_pos + nm(cell));
	// 		hydro(cell,hidx.u) += 0.5*sqr(hydro(cell,hidx.npy))/(tw::small_pos + nm(cell));
	// 		hydro(cell,hidx.u) += 0.5*sqr(hydro(cell,hidx.npz))/(tw::small_pos + nm(cell));
	// 	}
	// }
}
