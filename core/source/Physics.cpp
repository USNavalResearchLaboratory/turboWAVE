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

tw::Float sparc::CoulombCrossSectionCGS(tw::Float q1,tw::Float q2,tw::Float m12,tw::Float v12,tw::Float N1,tw::Float N2,tw::Float T1,tw::Float T2)
{
	// temperatures in ergs
	tw::Float rmin,rmax,rmin_alt,coulombLog;
	rmin = fabs(q1*q2)/(tw::small_pos + m12*v12*v12);
	rmin_alt = cgs::hbar/(tw::small_pos + 2*m12*v12);
	rmin = (rmin*rmin_alt)/(rmin+rmin_alt); // efficiently estimate the smaller of the two
	rmax = 1/sqrt(tw::small_pos + 4*pi*N1*q1*q1/T1 + 4*pi*N2*q2*q2/T2);
	coulombLog = log((3*rmin+rmax)/rmin); // well behaved approximation of log(rmax/rmin)
	return (32/pi)*pow(v12,-4)*sqr(q1*q2/m12)*coulombLog;
}

tw::Float sparc::ElectronPhononFrequencyCGS(tw::Float Ti,tw::Float EFermi,tw::Float ks)
{
	// temperatures in ergs
	const tw::Float vF = sqrt(2*EFermi/cgs::me);
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
	nstar = Z / sqrt(2*Uion);
}

void Ionizer::DeduceStaticCoeff()
{
	I1 = A1 * pow(2*Uion,0.75) / sqrt(3/pi);
	I2 = A2 - 0.5;
	I3 = A3;
}

void Ionizer::DeduceAveragedCoeff()
{
	A1 = I1 * sqrt(3/pi) / pow(2*Uion,0.75);
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
	return A1*w0*pow(fabs(E)/E_MPI,two*photons) / Factorial(photons-1);
}

tw::Float Tunneling::InstantRate(tw::Float w0,tw::Float E)
{
	const tw::Float Ea = (fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
	return I1*pow(Ea,I2)*exp(I3/Ea);
}

tw::Float Tunneling::AverageRate(tw::Float w0,tw::Float E)
{
	const tw::Float Ea = (fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
	return A1*pow(Ea,A2)*exp(A3/Ea);
}

void KYH::Initialize()
{
	Ionizer::Initialize();
	const tw::Float alpha = 0.0072973525693;
	const tw::Float Ua2 = Uion*alpha*alpha;
	const tw::Float Ea = pow(2*Uion,1.5);
	A1 = multiplier * pow(2.0,3.0-4*Ua2) * sqrt(3/pi) * (1.0-7*Ua2/72) * exp(4*Ua2) * pow(2*Uion,1.75-3.0*Ua2);
	A1 /= tgamma(3.0-2*Ua2);
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
	m = fabs(m);
	const tw::Float e = exp(1.0);
	const tw::Float sn2l2 = sqrt(nstar*nstar-lstar*lstar);
	A1 = sqrt(3/cub(pi)) * (2*l+1) * tgamma(l+m+1) / tgamma(m+1) / tgamma(l-m+1);
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
	m = fabs(m);
	const tw::Float F0 = cub(Z/nstar);
	// First without the Coulomb factor
	A1 = Uion;
	A1 *= pow(2,2*nstar) / (nstar*tgamma(nstar+lstar+1)*tgamma(nstar-lstar)); // |C|^2
	A1 *= sqrt(6/pi);
	A1 *= (2*l+1) * tgamma(l+m+1) / pow(2,m) / tgamma(m+1) / tgamma(l-m+1);
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
	A1 /= nstar*tgamma(nstar+lstar+1)*tgamma(nstar-lstar);
	// use A3 to hold F0
	A3 = pow(2*Uion,tw::Float(1.5));
}

tw::Float PPT::FourierSum(tw::Float gam,tw::Float nu)
{
	tw::Float A0 = 0.0;
	const tw::Float alpha = 2*(asinh(gam) - gam/sqrt(1 + gam*gam));
	const tw::Float beta = 2*gam/sqrt(1 + gam*gam);
	const tw::Float dnu = MyCeil(nu)-nu;
	for (tw::Int n=0;n<terms;n++)
		A0 += exp(-alpha*(n+dnu))*tw::dawsoni(sqrt(beta*(n+dnu)));
	return A0;
}

tw::Float PPT::AverageRate(tw::Float w0,tw::Float E)
{
	tw::Float ans;
	const tw::Float wa = (w0*tw::dims::angular_frequency >> native >> atomic);
	const tw::Float Ea = (fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
	const tw::Float F = Ea/A3;
	const tw::Float gam = sqrt(2*Uion)*wa/Ea;
	const tw::Float g = (3/(2*gam))*((1 + 1/(2*gam*gam))*asinh(gam) - sqrt(1 + gam*gam)/(2*gam));
	const tw::Float nu = (Uion/wa) * (1 + 1/(2*gam*gam));
	ans = Uion*A1*sqrt(6/pi)*(2*l+1);
	ans *= pow(F*sqrt(1 + gam*gam)/2,tw::Float(1.5));
	ans *= (4/sqrt(3*pi)) * (gam*gam/(1 + gam*gam)) * FourierSum(gam,nu);
	ans *= exp(-2*g/(3*F));
	ans *= pow(2/F,2*nstar); // coulomb correction
	return ans*tw::dims::angular_frequency >> atomic >> native;
}

void PMPB::Initialize()
{
	Ionizer::Initialize();
	// use A1 to hold C_nl^2, n.b. PMPB has a different convention from PPT (factor of 4)
	A1 = pow(2,2*nstar-2);
	A1 /= nstar*tgamma(nstar+lstar+1)*tgamma(nstar-lstar);
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
	const tw::Float Ea = (fabs(E)*tw::dims::electric_field >> native >> atomic) + cutoff_field;
	const tw::Float F = Ea/A3;
	const tw::Float gam = sqrt(2*Uion)*wa/Ea;
	const tw::Float g = (3/(2*gam))*((1 + 1/(2*gam*gam))*asinh(gam) - sqrt(1 + gam*gam)/(2*gam));
	const tw::Float nu = (Uion/wa) * (1 + 1/(2*gam*gam));
	ans = (2/pi)*Uion*A1*(2*l+1);
	ans *= pow(Uion/wa,-1.5);
	ans *= sqrt(2*gam/sqrt(1+gam*gam)) * FourierSum(gam,nu);
	ans *= exp(-2*g/(3*F));
	// Following is the improved Coulomb correction factor
	ans *= pow(2/F,2*nstar) * pow(1+2*gam/exp(1),-2*nstar);
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
					eos(cell,eidx.P) += a*u + ((b*u/denom) + A*mew*exp(-beta*expo))*exp(-alpha*sqr(expo));
					break;
				case 4:
					eos(cell,eidx.P) += (a + b/denom)*u + A*mew;
					break;
				case 5: // this is an interpolation of region 2 and 3
					P2 = (a + b/denom)*u + A*mew + B*sqr(mew);
					expo = (rho0/rho) - 1.0;
					P3 = a*u + ((b*u/denom) + A*mew*exp(-beta*expo))*exp(-alpha*sqr(expo));
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
