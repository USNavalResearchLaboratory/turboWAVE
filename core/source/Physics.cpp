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
	directives.Add("vibrational energy",new tw::input::Float(&excitationEnergy),false);
	directives.Add("thermometric conductivity",new tw::input::Float(&thermometricConductivity),false);
	directives.Add("kinematic viscosity",new tw::input::Float(&kinematicViscosity),false);
	directives.Add("permittivity",new tw::input::Numbers<tw::Float>(&eps[0],2),false);
}

tw::Float sparc::CoulombCrossSection(const UnitConverter& uc,tw::Float q1,tw::Float q2,tw::Float m12,tw::Float v12,tw::Float N1,tw::Float N2,tw::Float T1,tw::Float T2)
{
	tw::Float rmin,rmax,rmin_alt,coulombLog,hbar;
	hbar = uc.hbar / (uc.MKSValue(time_dim) * uc.MKSValue(energy_dim));
	rmin = fabs(q1*q2)/(tw::small_pos + 4*pi*m12*v12*v12*uc.MKSValue(number_dim));
	rmin_alt = hbar/(tw::small_pos + 2*m12*v12);
	if (rmin_alt<rmin) rmin = rmin_alt;
	rmax = 1/sqrt(tw::small_pos + N1*q1*q1/T1 + N2*q2*q2/T2);
	coulombLog = log(rmax/rmin);
	if (coulombLog<1.0) coulombLog = 1.0;
	return (32/pi)*pow(v12,-4)*sqr(q1*q2/(4*pi*m12))*coulombLog/uc.MKSValue(number_dim);
}

tw::Float sparc::ElectronPhononRateCoeff(const UnitConverter& uc,tw::Float Ti,tw::Float EFermi,tw::Float ks,tw::Float nref)
{
	// here, the rate coefficient is collision frequency divided by a reference density
	// this allows us to multiply away the collisions in regions where the metallic density is at "background level"
	// get quantities into cgs
	tw::Float vF = uc.c*100*sqrt(2*EFermi);
	tw::Float qe = uc.SimToCGS(charge_dim,1);
	tw::Float kB_Ti = uc.SimToCGS(energy_dim,Ti);
	tw::Float hbar = 1e7 * uc.hbar;
	// collision frequency in real units
	tw::Float nu = 2*ks*qe*qe*kB_Ti/(sqr(hbar)*vF);
	// return normalized rate coefficient
	return (uc.CGSValue(time_dim)/nref) * nu;
}

//////////////////////////
//                      //
//   Photoionization    //
//                      //
//////////////////////////


Ionizer::Ionizer(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	ionizationPotential = 1e-5; // in simulation units
	electrons = 0;
	protons = 0;
	multiplier = 1.0;
	I1 = I2 = I3 = A1 = A2 = A3 = 0.0;
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
	if (space->units==NULL)
		throw tw::FatalError("Ionizer tool requires units (set unit density).");
	Z = protons - electrons + 1;
	Uion = space->units->SimToAtomic(energy_dim,ionizationPotential);
	nstar = Z / sqrt(2*Uion);
	lstar = nstar - 1;
}

Multiphoton::Multiphoton(const std::string& name,MetricSpace *m,Task *tsk) : Ionizer(name,m,tsk)
{
	directives.Add("reference field",new tw::input::Float(&E_MPI));
}

void Multiphoton::Initialize()
{
	A3 = space->units->SimToAtomic(electric_field_dim,E_MPI);
}

tw::Float Multiphoton::AverageRate(tw::Float w0,tw::Float E)
{
	const tw::Float wa = space->units->SimToAtomic(angular_frequency_dim,w0); // laser freq. in a.u.
	const tw::Float photons = MyFloor(Uion/wa + 1);
	return multiplier*two*pi*w0*pow(fabs(E)/E_MPI,two*photons) / Factorial(photons-1);
}

ADK::ADK(const std::string& name,MetricSpace *m,Task *tsk) : Ionizer(name,m,tsk)
{
}

void ADK::Initialize()
{
	Ionizer::Initialize();
	I1 = multiplier*pow(4*exp(1)*pow(Z,3)/pow(nstar,4),2*nstar)/(8*pi*Z);
	I2 = 1 - 2*nstar;
	I3 = 2*pow(2*Uion,tw::Float(1.5))/3;
	A1 = I1 * sqrt(3*cub(nstar)/(pi*cub(Z)));
	A2 = I2 + tw::Float(0.5);
	A3 = I3;
	I1 = space->units->AtomicToSim(angular_frequency_dim,I1);
	A1 = space->units->AtomicToSim(angular_frequency_dim,A1);
}

tw::Float ADK::InstantRate(tw::Float w0,tw::Float E)
{
	const tw::Float Ea = space->units->SimToAtomic(electric_field_dim,fabs(E)) + 0.01;
	return I1*pow(Ea,I2)*exp(-I3/Ea);
}

tw::Float ADK::AverageRate(tw::Float w0,tw::Float E)
{
	const tw::Float Ea = space->units->SimToAtomic(electric_field_dim,fabs(E)) + 0.01;
	return A1*pow(Ea,A2)*exp(-A3/Ea);
}

PPT::PPT(const std::string& name,MetricSpace *m,Task *tsk) : Ionizer(name,m,tsk)
{
	terms = 1;
	directives.Add("terms",new tw::input::Int(&terms));
}

tw::Float PPT::wfunc(tw::Float x)
{
	// the argument of this function varies from 0 up to sqrt(2*terms)
	// where terms is the number of terms kept in the ppt expansion
	return tw::Float(0.5) * exp(-x*x) * tw::erfi(x) * sqrt(pi);
}

void PPT::Initialize()
{
	Ionizer::Initialize();
	const tw::Float F0 = pow(2*Uion,tw::Float(1.5));
	A1 =  multiplier * Uion * (4/sqrt(3*pi)) * sqrt(6/pi) * pow(2,2*nstar) / (nstar*tgamma(nstar+lstar+1)*tgamma(nstar-lstar));
	A3 = 2*F0/3;
	A1 = space->units->AtomicToSim(angular_frequency_dim,A1);
}

tw::Float PPT::AverageRate(tw::Float w0,tw::Float E)
{
	tw::Float ans = 0.0;
	const tw::Float wa = space->units->SimToAtomic(angular_frequency_dim,w0);
	const tw::Float Ea = space->units->SimToAtomic(electric_field_dim,fabs(E)) + 0.01;
	const tw::Float F0 = pow(2*Uion,tw::Float(1.5));
	const tw::Float gam = sqrt(2*Uion)*wa/Ea;
	const tw::Float gam2 = gam*gam;
	const tw::Float alpha = 2*(asinh(gam) - gam/sqrt(1 + gam2));
	const tw::Float beta = 2*gam/sqrt(1 + gam2);
	const tw::Float g = (3/(2*gam))*((1 + 1/(2*gam2))*asinh(gam) - sqrt(1 + gam2)/(2*gam));
	const tw::Float nu = (Uion/wa)*(1 + 1/(2*gam2));
	const tw::Float dnu = MyCeil(nu)-nu;
	for (tw::Int n=0;n<terms;n++)
		ans += exp(-alpha*(n+dnu))*wfunc(sqrt(beta*(n+dnu)));
	ans *= (gam2/(1 + gam2));
	ans *= pow(Ea*sqrt(1 + gam2)/(2*F0),tw::Float(1.5));
	ans *= pow(2*F0/Ea,2*nstar); // coulomb correction
	ans *= exp(-A3*g/Ea);
	ans *= A1;
	return ans;
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
	// arbitrary normalization value I used when testing this EOS
	//unitDensityCGS = 2.8e19; 

	// Tillotson parameters for H20
	n0 = 1193.0;

	a = 0.7;   // Tillotson Coefficient
	b = 0.15;   // Tillotson Coefficient
	A = 0.000950973;   // Tillotson Coefficient [sim - pressure unit]
	B = 0.00577999;   // Tillotson Coefficient [sim - pressure unit]
	alpha = 10.0;   // Tillotson Coefficient
	beta = 5.0;   // Tillotson Coefficient

	nIV = 1144.06;   // vaporization density 
	E0 = 2.55698e-06;   // Reference energy  
	EIV = 1.53054e-07;   // Vaporization Energy 
	ECV = 9.13208e-07;   // Cavitation Energy 

	// std::cout << "Start Creating Till tool" << std::endl;

	//directives.Add("unit density",new tw::input::Float(&unitDensityCGS));
	directives.Add("reference density",new tw::input::Float(&n0));

	directives.Add("parameter a",new tw::input::Float(&a));
	directives.Add("parameter b",new tw::input::Float(&b));
	directives.Add("parameter A",new tw::input::Float(&A));
	directives.Add("parameter B",new tw::input::Float(&B));
	directives.Add("parameter alpha",new tw::input::Float(&alpha));
	directives.Add("parameter beta",new tw::input::Float(&beta));

	directives.Add("vaporization density",new tw::input::Float(&nIV));
	directives.Add("reference energy",new tw::input::Float(&E0));
	directives.Add("vaporization energy",new tw::input::Float(&EIV));
	directives.Add("cavitation energy",new tw::input::Float(&ECV));

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

	UnitConverter *uc = space->units;

	// threshold parameters in cgs units
	tw::Float RhoIV_cgs = uc->SimToCGS(mass_dim,mat.mass)*uc->SimToCGS(density_dim,nIV); // Vaporization Pressure [g/cm3]
	tw::Float E0_cgs = uc->SimToCGS(energy_dim,E0)/uc->SimToCGS(mass_dim,mat.mass); // Reference energy density [erg/g]
	tw::Float EIV_cgs = uc->SimToCGS(energy_dim,EIV)/uc->SimToCGS(mass_dim,mat.mass); // Vaporization Energy [erg/g]
	tw::Float ECV_cgs = uc->SimToCGS(energy_dim,ECV)/uc->SimToCGS(mass_dim,mat.mass); // Cavitation Energy [erg/g]

	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
		{
			const tw::Float ndens = hydro(cell,hidx.ni);
			const tw::Float udens = hydro(cell,hidx.u);
			const tw::Float rho_cgs = uc->SimToCGS(mass_dim,mat.mass)*uc->SimToCGS(density_dim,ndens);
			const tw::Float u_cgs = uc->SimToCGS(energy_density_dim,udens);

			const tw::Float rho0_cgs = uc->SimToCGS(mass_dim,mat.mass)*uc->SimToCGS(density_dim,n0);
			const tw::Float u0_cgs = rho_cgs*E0_cgs + tw::small_pos; // U [etg/cm3] = rho [g/cm3] * E [erg/g]

			const tw::Float eta = ndens/n0; // compression
			const tw::Float mew = eta - 1.0; // strain

			// Determine Region
			tw::Int region = 0; // 1, 2, 3, 4, or 5 : 0 is for error detection
			if ( rho_cgs >= rho0_cgs ) region = 1;
			if ( ((rho_cgs >= RhoIV_cgs) and (rho_cgs < rho0_cgs)) and (u_cgs <= rho_cgs*EIV_cgs) ) region = 2;
			if ((rho_cgs < rho0_cgs) and (u_cgs >= rho_cgs*ECV_cgs)) region = 3;
			if ((rho_cgs < RhoIV_cgs) and (u_cgs < rho_cgs*ECV_cgs)) region = 4;
			if (((rho_cgs > RhoIV_cgs) and (rho_cgs < rho0_cgs)) and ((u_cgs > rho_cgs*EIV_cgs) and (u_cgs < rho_cgs*ECV_cgs))) {
				region = 5;
			}
			if (region == 0) {
				std::cout << "ERROR: Unrecognized Region Detected in Tillotson EOS." << std::endl;
				std::cout << "rho_cgs = " << rho_cgs << std::endl;
				std::cout << "E_cgs = " << (u_cgs/rho_cgs) << std::endl;
				exit(1); 
			}

			// Pressure Calculation
			// eos(cell,eidx.P) += ndens*eos(cell,eidx.T);
			tw::Float denom = (u_cgs/(u0_cgs*sqr(eta))) + 1.0; // this quantity is repeated in expressions
			tw::Float expo, P3, P2, P_sim; // might be used depending on state
			switch (region)
			{
				case 1:
					eos(cell,eidx.P) += (a + b/denom)*udens + A*mew + B*sqr(mew);
					break;
				case 2:
					eos(cell,eidx.P) += (a + b/denom)*udens + A*mew + B*sqr(mew);
					break;
				case 3:
					expo = (rho0_cgs/rho_cgs) - 1.0;
					eos(cell,eidx.P) += a*udens + ((b*udens/denom) + A*mew*exp(-beta*expo))*exp(-alpha*sqr(expo));
					break;
				case 4:
					eos(cell,eidx.P) += (a + b/denom)*udens + A*mew;
					break;
				case 5: // this is an interpolation of region 2 and 3
					P2 = (a + b/denom)*udens + A*mew + B*sqr(mew);
					expo = (rho0_cgs/rho_cgs) - 1.0;
					P3 = a*udens + ((b*udens/denom) + A*mew*exp(-beta*expo))*exp(-alpha*sqr(expo));
					eos(cell,eidx.P) += ((u_cgs - rho_cgs*EIV_cgs)*P3 + (rho_cgs*ECV_cgs - u_cgs)*P2)/(rho_cgs*(ECV_cgs-EIV_cgs));
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
			eos(cell,eidx.T) = eosRef(cell,eidx.T) + (IE1*(1+nm0/nm1) - IE0*(1+nm1/nm0))/nmcv_sum;
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
			eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
			nm(cell) = nm1;
			IE(cell) = IE1;
		}
	}
}
