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
import numerics;
import logger;

/// EOS Component.
/// Maintains indexing information for accessing hydro and eos fields, and material parameters.
/// Defaults to an ideal gas.

/**
 * @brief EOS for one component of a mixture, its main task to to compute pressure
 * 
 */
export struct EOSComponent:ComputeTool
{
	tw::Float nm_ref; // reference mass density
	tw::Float E_ref; // reference specific energy (energy/mass)
	sparc::hydro_set hidx;
	sparc::eos_set eidx;
	sparc::material mat;

	EOSComponent(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk) {
		nm_ref = 0.0;
		E_ref = 0.0;
	}
	void SetupIndexing(tw::Int component_index,const sparc::hydro_set& h,const sparc::eos_set& e,const sparc::material& m)
	{
		hidx = h;
		hidx.ni = component_index;
		eidx = e;
		mat = m;
	}
	/**
	 * @brief Heat capacity (nmcv) for this component
	 * 
	 * @param n density to use in this calculation
	 * @param IE internal energy density to use in this calculation
	 * @return heat capacity (nmcv) at constant volume (energy/volume/temperature)
	 */
	virtual tw::Float HeatCapacity(tw::Float nm, tw::Float IE) {
		return nm * mat.cvm / mat.mass;
	}
	/**
	 * @brief Internal energy density at absolute zero
	 * 
	 * @param n density to use in this calculation
	 * @return internal energy density
	 */
	virtual tw::Float ColdCurve(tw::Float nm) {
		return 0.0;
	}
	/**
	 * @brief Internal energy density at any temperature
	 * 
	 * @param n density to use in this calculation
	 * @param T temperature to use in this calculation
	 * @return internal energy density
	 */
	virtual tw::Float InternalEnergy(tw::Float nm, tw::Float T) {
		return ColdCurve(nm) + HeatCapacity(nm,T) * T;
	}
	/**
	 * @brief Partial pressure for this component
	 * 
	 * @param IE internal energy density assigned to this component
	 * @param nm mass density of this component
	 * @return pressure
	 */
	virtual tw::Float Pressure(tw::Float nm, tw::Float IE) {
		return IE / mat.cvm;
	}
	/**
	 * @brief Add pressure, heat conductivity, and viscosity to the EOS field
	 * 
	 * @param[in] nm partial mass density
	 * @param[in] IE partial internal energy density
	 * @param[in] nu_e collision frequency
	 * @param[in] hydro hydro data to use (n,np,u)
	 * @param[out] eos EOS field to update
	 */
    virtual void AddPKV(ScalarField& nm, ScalarField& IE, ScalarField& nu_e, Field& hydro, Field& eos)
    {
        #pragma omp parallel
        {
            for (auto cell : EntireCellRange(*space,1))
            {
                eos(cell,eidx.P) += Pressure(nm(cell), IE(cell));
                eos(cell,eidx.K) += mat.thermometricConductivity * mat.cvm * nm(cell) / mat.mass;
                eos(cell,eidx.visc) += mat.kinematicViscosity * nm(cell);
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
	virtual void AddPKV(ScalarField& nm, ScalarField& IE, ScalarField& nu_e, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				const tw::Float ne = hydro(cell,hidx.ni);
				eos(cell,eidx.P) += Pressure(nm(cell), IE(cell));
				eos(cell,eidx.K) += 3.2*ne*eos(cell,eidx.T)/(mat.mass*nu_e(cell));
				//eos(cell,eidx.visc) += 0.0; // don't touch, may help caching.
				// Braginskii has for e-viscosity 0.73*ne*eos(cell,eidx.T)/nu_e(cell)
				// However, we are forcing electrons to move with ions and so should not diffuse velocity field
			}
		}
	}
};

/**
 * @brief EOS with built-in integrator to work out the cold curve
 * 
 */
export struct EOSCondensedMatter:EOSComponent
{
	EOSCondensedMatter(const std::string& name,MetricSpace *m,Task *tsk) : EOSComponent(name,m,tsk) {}
	virtual tw::Float ColdCurve(tw::Float nm) {
		auto dIEdnm = [this] (tw::Float nm,tw::Float IE) {
			// In the variables (P,V,E) we have P = dE/dV; in variables (P,n,IE) it becomes (P+IE)/n = du/dn
			// Here, IE = nE, V = 1/n
			return (Pressure(nm,IE) + IE) / nm;
		};
		return RK4Integrate<tw::Float>(E_ref*nm_ref, nm_ref, nm, (nm-nm_ref)/8, dIEdnm, 1e-7);
	}
};

/**
 * @brief Bare bones Mie-Gruneisen P = GRUN*IE
 * 
 */
export struct EOSSimpleMieGruneisen:EOSCondensedMatter
{
	tw::Float GRUN; // Gruneisen coefficient
	EOSSimpleMieGruneisen(const std::string& name,MetricSpace *m, Task *tsk) : EOSCondensedMatter(name,m,tsk)
	{
		GRUN = 2.0; // value for Cu on p. 257 of "Shock Wave Physics and Equation of State Modeling"
		// GRUN = 0.1; // value for water in the above book.
		directives.Add("gruneisen parameter",new tw::input::Float(&GRUN),true);
	}
	virtual tw::Float Pressure(tw::Float nm, tw::Float IE) {
		return GRUN*IE;
	}
};

/**
 * @brief Mie-Gruneisen EOS that uses a linear fit to the Hugoniot
 * 
 */
export struct EOSLinearMieGruneisen:EOSCondensedMatter
{
	tw::Float GRUN; // Gruneisen coefficient
	tw::Float c0;   // y - intercept of Hugoniot fit (usually approximately speed of sound)
	tw::Float S1;   // coefficient of linear fit of Hugoniot data
	EOSLinearMieGruneisen(const std::string& name,MetricSpace *m, Task *tsk) : EOSCondensedMatter(name,m,tsk)
	{
		GRUN = 2.0; // value for Cu on p. 257 of "Shock Wave Physics and Equation of State Modeling"
		// GRUN = 0.1; // value for water in the above book.

		directives.Add("gruneisen parameter",new tw::input::Float(&GRUN),true);
		directives.Add("reference mass density",new tw::input::Float(&nm_ref),true);
		directives.Add("hugoniot intercept",new tw::input::Float(&c0),true);
		directives.Add("hugoniot slope", new tw::input::Float(&S1),true);
	}
	virtual tw::Float Pressure(tw::Float nm, tw::Float IE) {
		// Pressure from <http://bluevistasw.com/2016/02/16/mie-gruneisen-eos-implementation/>
		const tw::Float mu = nm/(nm_ref) - 1;
		const tw::Float sel = tw::Float(mu>=0.0);
		return (1-sel)*(nm_ref*c0*c0*mu + GRUN*(mu+1)*IE) +
			sel*(nm_ref*c0*c0*mu*(1 + (1 - GRUN*(mu+1)/2)*mu)/(1 - (S1-1)*mu) + GRUN*(mu+1)*IE);
	}
};

/// Tillotson EOS for modeling vaporization, cavitation, and shocks.
/// Coefficients for water can be found at [A.L. Brundage, Procedia Engineering (2013)]
export struct EOSTillotson:EOSCondensedMatter
{
	// Tillotson rho0 and E0 are provided by inherited nm_ref and E_ref

	tw::Float a;   // Tillotson Coefficient
	tw::Float b;   // Tillotson Coefficient
	tw::Float A;   // Bulk Modulus [pressure]
	tw::Float B;   // Tillotson Parameter [pressure]
	tw::Float alpha;   // Tillotson Coefficient
	tw::Float beta;   // Tillotson Coefficient

	tw::Float rhoIV;   // Incipient vaporization density
	tw::Float EIV;   // Incipient vaporization specific energy
	tw::Float ECV;   // Complete vaporization specific energy

	EOSTillotson(const std::string& name,MetricSpace *m, Task *tsk) : EOSCondensedMatter(name,m,tsk)
	{
		// Tillotson parameters for H20
		nm_ref = tw::dnum("0.998 [g/cm3]") >> native;

		a = 0.7;   // Tillotson Coefficient
		b = 0.15;   // Tillotson Coefficient
		A = tw::dnum("21.8e3 [bar]") >> native;   // Tillotson Coefficient, pressure
		B = tw::dnum("132.5e3 [bar]") >> native;   // Tillotson Coefficient, pressure
		alpha = 10.0;   // Tillotson Coefficient
		beta = 5.0;   // Tillotson Coefficient

		rhoIV = tw::dnum("0.958 [g/cm3]") >> native;   // Incipient vaporization density
		E_ref = tw::dnum("0.07e12 [ergs/g]") >> native;   // Reference energy
		EIV = tw::dnum("0.00419e12 [ergs/g]") >> native;   // Incipient vaporization specific energy
		ECV = tw::dnum("0.025e12 [ergs/g]") >> native;   // Complete vaporization specific energy

		directives.Add("reference mass density",new tw::input::Float(&nm_ref));

		directives.Add("parameter a",new tw::input::Float(&a));
		directives.Add("parameter b",new tw::input::Float(&b));
		directives.Add("pressure A",new tw::input::Float(&A));
		directives.Add("pressure B",new tw::input::Float(&B));
		directives.Add("parameter alpha",new tw::input::Float(&alpha));
		directives.Add("parameter beta",new tw::input::Float(&beta));

		directives.Add("incipient vaporization mass density",new tw::input::Float(&rhoIV));
		directives.Add("reference specific energy",new tw::input::Float(&E_ref));
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
	virtual tw::Float Pressure(tw::Float nm, tw::Float IE) {
		const tw::Float rho = nm;
		const tw::Float u0 = rho*E_ref + tw::small_pos;

		const tw::Float eta = rho/nm_ref; // compression
		const tw::Float mew = eta - 1.0; // strain

		// Determine Region
		tw::Int region = 0; // 1, 2, 3, 4, or 5 : 0 is for error detection
		if ( rho >= nm_ref && IE > 0.0 ) region = 1; // eq. 1
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
		const tw::Float expo = nm_ref/rho - 1;
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
		REGISTER(EOSTillotson,ColdCurveTest);
	}
	std::tuple<Field,Field,ScalarField,ScalarField,ScalarField> InitTest();
	void PressureTest();
	void ColdCurveTest();
};
