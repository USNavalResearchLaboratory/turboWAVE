
module;

#include "tw_includes.h"
#include "tw_test.h"

/// Module to handle equation of state calculations.
///
/// The standard update pattern is
/// 1. Components load the nmcv array
/// 2. Mixture computes T, and returns IE and nm
/// 3. Components load P, K, and visc arrays
export module eos;

import input;
import driver;
import fields;
import functions;
import hydro_primitives;

/// EOS Component.
/// Maintains indexing information for accessing hydro and eos fields, and material parameters.
/// Defaults to an ideal gas.
export struct EOSComponent:ComputeTool
{
	sparc::hydro_set hidx;
	sparc::eos_set eidx;
	sparc::material mat;

	EOSComponent(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk) {}
	void SetupIndexing(tw::Int component_index,const sparc::hydro_set& h,const sparc::eos_set& e,const sparc::material& m)
	{
		hidx = h;
		hidx.ni = component_index;
		eidx = e;
		mat = m;
	}
    virtual void SetHeatCapacity(ScalarField& nm,Field& eos)
    {
        #pragma omp parallel
        {
            for (auto cell : EntireCellRange(*space,0))
                eos(cell,eidx.nmcv) = nm(cell) * mat.cvm / mat.mass;
        }
    }
    virtual void AddHeatCapacity(Field& hydro,Field& eos)
    {
        #pragma omp parallel
        {
            for (auto cell : EntireCellRange(*space,0))
                eos(cell,eidx.nmcv) += hydro(cell,hidx.ni) * mat.cvm;
        }
    }
    virtual void AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
    {
        #pragma omp parallel
        {
            for (auto cell : EntireCellRange(*space,0))
            {
                const tw::Float ngas = hydro(cell,hidx.ni);
                eos(cell,eidx.P) += ngas*eos(cell,eidx.T);
                eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * ngas;
                eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * ngas;
            }
        }
    }
};

/// Ideal gas component, just an alias for the base class
export struct EOSIdealGas:EOSComponent
{
	EOSIdealGas(const std::string& name,MetricSpace *m,Task *tsk) : EOSComponent(name,m,tsk) {
		// nothing to do, ideal gas is the default
	}
};

/// Braginskii model for plasma electrons, hot enough to ignore quantum effects
export struct EOSHotElectrons:EOSComponent
{
	EOSHotElectrons(const std::string& name,MetricSpace *m,Task *tsk) : EOSComponent(name,m,tsk) {}
	virtual void AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,0))
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
};

/// This is the most basic implementation of the MieGruneisen EOS
///
/// It assumes that GRUN = GRUN0 = const. at all times
/// This is typically not used in literature concerning MieGruneisen EOSs,
/// as the results are rarely physicaly accurate
/// If you produce sound waves with this model (for example with Cu),
/// you'll notice the sound speed is off. Qualitatively, it gives a broad picture.
export struct EOSSimpleMieGruneisen:EOSComponent
{
	tw::Float GRUN; // Gruneisen coefficient
	EOSSimpleMieGruneisen(const std::string& name,MetricSpace *m, Task *tsk) : EOSComponent(name,m,tsk)
	{
		GRUN = 2.0; // value for Cu on p. 257 of "Shock Wave Physics and Equation of State Modeling"
		// GRUN = 0.1; // value for water in the above book.
		directives.Add("gruneisen parameter",new tw::input::Float(&GRUN));
	}
	virtual void AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,0))
			{
				const tw::Float nion = hydro(cell,hidx.ni);
				const tw::Float partial_IE = IE(cell) * nion * mat.mass / nm(cell);
				eos(cell,eidx.P) += GRUN*partial_IE;
				eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * nion;
				eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * nion;
			}
		}
	}
};

/// This is a MieGruneisen EOS that assumes \rho * GRUN = const., and a linear Hugoniot fit
///
/// This is what is more typically what is found in literature regarding MieGruneisen models
/// You'll get a more accurate sound speed and shock speeds, assuming that the simulation is
/// within the range of relevant Hugoniot data and model assumptions are properly met
export struct EOSLinearMieGruneisen:EOSComponent
{
	tw::Float GRUN; // Gruneisen coefficient
	tw::Float n0;   // Reference density
	tw::Float c0;   // y - intercept of Hugoniot fit (usually appriximately speed of sound)
	tw::Float S1;   // coefficient of linear fit of Hugoniot data
	EOSLinearMieGruneisen(const std::string& name,MetricSpace *m, Task *tsk) : EOSComponent(name,m,tsk)
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
	virtual void AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,0))
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
};

/// Tillotson EOS for modeling vaporization, cavitation, and shocks.
/// Coefficients for water can be found at [A.L. Brundage, Procedia Engineering (2013)]
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

	EOSTillotson(const std::string& name,MetricSpace *m, Task *tsk) : EOSComponent(name,m,tsk)
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
	virtual void AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,0))
			{
				const tw::Float ndens = hydro(cell,hidx.ni);
				const tw::Float u = hydro(cell,hidx.u);
				const tw::Float rho = mat.mass * ndens;
				const tw::Float u0 = rho*E0 + tw::small_pos;

				const tw::Float eta = rho/rho0; // compression
				const tw::Float mew = eta - 1.0; // strain

				// Determine Region
				tw::Int region = 0; // 1, 2, 3, 4, or 5 : 0 is for error detection
				if ( rho >= rho0 && u > 0.0 ) region = 1; // eq. 1
				else if ( rho >= rhoIV && u <= rho*EIV ) region = 2; // eq . 2
				else if ( u >= rho*ECV ) region = 3; // eq. 3
				else if ( rho < rhoIV && u < rho*ECV ) region = 4; // eq. 6
				else if ( rho > rhoIV && u > rho*EIV && u < rho*ECV ) region = 5; // eq. 5

				if (region == 0) {
					std::stringstream err_mess;
					err_mess << "Unrecognized Region Detected in Tillotson EOS." << std::endl;
					err_mess << "rho = " << (rho*tw::dims::mass_density>>native>>cgs) << " [g/cm3]" << std::endl;
					err_mess << "E = " << ((u/rho)*tw::dims::specific_energy>>native>>cgs) << " [ergs/g]" << std::endl;
					throw tw::FatalError(err_mess.str());
				}

				// Pressure Calculation
				const tw::Float denom = (u/(u0*sqr(eta))) + 1.0; // this quantity is repeated in expressions
				const tw::Float expo = rho0/rho - 1;
				const tw::Float P4 = (a + b/denom)*u + A*mew;
				switch (region)
				{
					case 1:
					case 2:
						eos(cell,eidx.P) += P4 + B*sqr(mew);
						break;
					case 3:
						eos(cell,eidx.P) += a*u + ((b*u/denom) + A*mew*std::exp(-beta*expo))*std::exp(-alpha*sqr(expo));
						break;
					case 4:
						eos(cell,eidx.P) += P4;
						break;
					case 5: // this is an interpolation of region 2 and 3
						const tw::Float P2 = P4 + B*sqr(mew);
						const tw::Float P3 = a*u + ((b*u/denom) + A*mew*std::exp(-beta*expo))*std::exp(-alpha*sqr(expo));
						eos(cell,eidx.P) += ((u - rho*EIV)*P3 + (rho*ECV - u)*P2)/(rho*(ECV-EIV));
						break;
				}

				// keep rest of EOS the same
				eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * ndens;
				eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * ndens;
			}
		}
	}
};

/// Tool for handling non-additive quantities in a mixture.
/// An important task is calculating temperature given various inputs.
/// Additive quantities like pressure tend to be handled by components.
export struct EOSMixture:ComputeTool
{
	sparc::hydro_set hidx;
	sparc::eos_set eidx;
	sparc::material_set matset;
	EOSMixture(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk) {}
	void SetupIndexing(const sparc::hydro_set& h,const sparc::eos_set& e,const sparc::material_set& m)
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
	/// Initialize the temperature without using any reference data.
	/// Inputs - `hydro`, `eos` component nmcv
	/// Outputs - `IE`, `nm`, `eos` components T and Tv
	virtual void InitTemperature(ScalarField& IE, ScalarField& nm, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,0))
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
	/// Update the temperature using given reference data (e.g. previous time level).
	/// Inputs - `hydro`, `hydroRef`, `eos` component nmcv, `eosRef`
	/// Outputs - `IE`, `nm`, `eos` components T and Tv
	virtual void UpdateTemperature(ScalarField& IE, ScalarField& nm, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,0))
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
	/// Add energy corresponding to a change in temperature only.
	virtual void UpdateEnergy(ScalarField& nm,ScalarField& T0,Field& hydro,Field& eos)
	{
		// Not centered, because cv is not updated.
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,0))
				hydro(cell,hidx.u) += eos(cell,eidx.nmcv) * (eos(cell,eidx.T) - T0(cell));
		}
	}
};

/// Mixture that uses a polytropic ideal gas caloric EOS, reference states are ignored
export struct EOSIdealGasMix:EOSMixture
{
	EOSIdealGasMix(const std::string& name,MetricSpace *m,Task *tsk) : EOSMixture(name,m,tsk) {}
	virtual void UpdateTemperature(ScalarField& IE, ScalarField& nm, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,0))
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
};
