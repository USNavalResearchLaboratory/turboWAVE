module;

#include "tw_includes.h"
#include "tw_test.h"
#include "tw_logger.h"

export module eos:eos_component;
import input;
import driver;
import fields;
import functions;
import hydro_primitives;
import logger;

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
    // virtual void SetHeatCapacity(ScalarField& nm,Field& eos)
    // {
    //     #pragma omp parallel
    //     {
    //         for (auto cell : EntireCellRange(*space,1))
    //             eos(cell,eidx.nmcv) = nm(cell) * mat.cvm / mat.mass;
    //     }
    // }
    virtual void AddHeatCapacity(Field& hydro,Field& eos)
    {
        #pragma omp parallel
        {
            for (auto cell : EntireCellRange(*space,1))
                eos(cell,eidx.nmcv) += hydro(cell,hidx.ni) * mat.cvm;
        }
    }
	/**
	 * @brief Heat capacity (nmcv) for this component
	 * 
	 * @param n density to use in this calculation
	 * @param T temperature to use in this calculation
	 * @return heat capacity at constant volume (energy/volume/temperature)
	 */
	virtual tw::Float HeatCapacity(tw::Float n, tw::Float T) {
		return n * mat.cvm;
	}
	/**
	 * @brief Internal energy density at absolute zero
	 * 
	 * @param n density to use in this calculation
	 * @return internal energy density
	 */
	virtual tw::Float ColdCurve(tw::Float n) {
		return 0.0;
	}
	/**
	 * @brief Internal energy density at any temperature
	 * 
	 * @param n density to use in this calculation
	 * @param T temperature to use in this calculation
	 * @return internal energy density
	 */
	virtual tw::Float InternalEnergy(tw::Float n, tw::Float T) {
		return ColdCurve(n) + HeatCapacity(n,T) * T;
	}
	/**
	 * @brief Partial pressure for this component
	 * 
	 * @param IE internal energy density assigned to this component
	 * @param n number density of this component
	 * @return pressure
	 */
	virtual tw::Float Pressure(tw::Float IE, tw::Float n) {
		return IE / mat.cvm;
	}
	/**
	 * @brief Add pressure, heat conductivity, and viscosity to the EOS field
	 * 
	 * @param[in] IE aggregated internal energy density
	 * @param[in] nm aggregated mass density
	 * @param[in] nu_e collision frequency
	 * @param[in] hydro hydro data to use (n,np,u)
	 * @param[out] eos EOS field to update
	 */
    virtual void AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
    {
        #pragma omp parallel
        {
            for (auto cell : EntireCellRange(*space,1))
            {
				const tw::Float n = hydro(cell,hidx.ni);
				const tw::Float partial_IE = IE(cell) * n * mat.mass / nm(cell);
                eos(cell,eidx.P) += Pressure(partial_IE, n);
                eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * n;
                eos(cell,eidx.visc) += mat.kinematicViscosity * mat.mass * n;
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
	virtual void RegisterTests() {
		REGISTER(EOSIdealGas,PressureTest);
	}
	std::tuple<Field,Field,ScalarField,ScalarField,ScalarField> InitTest(tw::dnum n, tw::dnum T);
	void PressureTest();
};

/// Braginskii model for plasma electrons, hot enough to ignore quantum effects
export struct EOSHotElectrons:EOSComponent
{
	EOSHotElectrons(const std::string& name,MetricSpace *m,Task *tsk) : EOSComponent(name,m,tsk) {}
	virtual void AddPKV(ScalarField& IE, ScalarField& nm, ScalarField& nu_e, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				const tw::Float ne = hydro(cell,hidx.ni);
				const tw::Float partial_IE = IE(cell) * ne * mat.mass / nm(cell);
				eos(cell,eidx.P) += Pressure(partial_IE, ne);
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
	virtual tw::Float Pressure(tw::Float IE, tw::Float n) {
		return GRUN*IE;
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
	tw::Float c0;   // y - intercept of Hugoniot fit (usually approximately speed of sound)
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
	virtual tw::Float Pressure(tw::Float IE, tw::Float n) {
		// Pressure from <http://bluevistasw.com/2016/02/16/mie-gruneisen-eos-implementation/>
		const tw::Float mu = n/n0 - 1;
		const tw::Float sel = tw::Float(mu>=0.0);
		return (1-sel)*(mat.mass*n0*c0*c0*mu + GRUN*(mu+1)*IE) +
			sel*(mat.mass*n0*c0*c0*mu*(1 + (1 - GRUN*(mu+1)/2)*mu)/(1 - (S1-1)*mu) + GRUN*(mu+1)*IE);
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
	virtual tw::Float Pressure(tw::Float IE, tw::Float n) {
		const tw::Float rho = mat.mass * n;
		const tw::Float u0 = rho*E0 + tw::small_pos;

		const tw::Float eta = rho/rho0; // compression
		const tw::Float mew = eta - 1.0; // strain

		// Determine Region
		tw::Int region = 0; // 1, 2, 3, 4, or 5 : 0 is for error detection
		if ( rho >= rho0 && IE > 0.0 ) region = 1; // eq. 1
		else if ( rho >= rhoIV && IE <= rho*EIV ) region = 2; // eq . 2
		else if ( IE >= rho*ECV ) region = 3; // eq. 3
		else if ( rho < rhoIV && IE < rho*ECV ) region = 4; // eq. 6
		else if ( rho > rhoIV && IE > rho*EIV && IE < rho*ECV ) region = 5; // eq. 5
		if (region == 0) {
			std::stringstream err_mess;
			err_mess << "Unrecognized Region Detected in Tillotson EOS." << std::endl;
			err_mess << "rho = " << (rho*tw::dims::mass_density>>native>>cgs) << " [g/cm3]" << std::endl;
			err_mess << "E = " << ((IE/rho)*tw::dims::specific_energy>>native>>cgs) << " [ergs/g]" << std::endl;
			throw tw::FatalError(err_mess.str());
		}

		// Pressure Calculation
		const tw::Float denom = IE/(u0*sqr(eta)) + 1.0;
		const tw::Float expo = rho0/rho - 1;
		const tw::Float P4 = (a + b/denom)*IE + A*mew;
		switch (region)
		{
			case 1:
			case 2:
				return P4 + B*sqr(mew);
			case 3:
				return a*IE + ((b*IE/denom) + A*mew*std::exp(-beta*expo))*std::exp(-alpha*sqr(expo));
			case 4:
				return P4;
			case 5: // this is an interpolation of region 2 and 3
				const tw::Float P2 = P4 + B*sqr(mew);
				const tw::Float P3 = a*IE + ((b*IE/denom) + A*mew*std::exp(-beta*expo))*std::exp(-alpha*sqr(expo));
				return ((IE - rho*EIV)*P3 + (rho*ECV - IE)*P2)/(rho*(ECV-EIV));
		}
		return 0; // unreachable
	}

	virtual void RegisterTests() {
		REGISTER(EOSTillotson,PressureTest);
	}
	std::tuple<Field,Field,ScalarField,ScalarField,ScalarField> InitTest();
	void PressureTest();
};
