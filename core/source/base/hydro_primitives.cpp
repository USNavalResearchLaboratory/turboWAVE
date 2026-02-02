module;

#include "tw_includes.h"

export module hydro_primitives;
import input;
import fields;

export namespace sparc
{
	enum laserModel { vacuum, isotropic };
	enum radiationModel { noRadiation, thin, thick };
	enum plasmaModel { neutral, quasineutral };

	struct hydro_set
	{
		// num = number of constituent chemicals
		// first = index of first hydro quantity in the set, must be the first chemical density
		// last = index of last hydro quantity in the set
		// ni = index to a particular chemical density (optionally used)
		// npx,npy,npz = index to momentum density
		// u = index to energy density
		// x = index to vibrational density
		tw::Int num,first,last,ni,npx,npy,npz,u,x;
		tw::Int Load(tw::Int i,tw::Int N)
		{
			num = N;
			ni = i;
			first = i;
			npx = i+N;
			npy = i+N+1;
			npz = i+N+2;
			u = i+N+3;
			x = i+N+4;
			last = x;
			return last+1;
		}
		tw::Float DensitySum(const Field& f,const tw::cell& cell) const
		{
			tw::Float ans = 0.0;
			for (tw::Int i=0;i<num;i++)
				ans += f(cell,first+i);
			return ans;
		}
		static const tw::Int count = 5; // gives count of state components not counting densities
	};

	struct eos_set
	{
		// T = index to translational and rotational temperature
		// Tv = vibrational temperature (if used)
		// P = pressure
		// K = heat conductivity
		// visc = dynamic viscosity
		// nmcv = mass density weighted heat capacity at constant volume
		tw::Int T,Tv,P,K,visc,nmcv;
		tw::Int Load(tw::Int i) { T=i;Tv=i+1;P=i+2;K=i+3;visc=i+4;nmcv=i+5;return i+6; }
		static const tw::Int count = 6;
	};

	/// Characteristics of a Chemical object
	struct material
	{
		tw::Float mass,charge,cvm,excitationEnergy,thermometricConductivity,kinematicViscosity,eps[2];
		/// The chemical has to call this to fully setup its input file parameters
		void AddDirectives(tw::input::DirectiveReader& directives) {
			directives.Add("mass",new tw::input::Float(&mass));
			directives.Add("charge",new tw::input::Float(&charge));
			directives.Add("cv",new tw::input::Float(&cvm));
			directives.Add("vibrational energy",new tw::input::Float(&excitationEnergy),false);
			directives.Add("thermometric conductivity",new tw::input::Float(&thermometricConductivity),false);
			directives.Add("kinematic viscosity",new tw::input::Float(&kinematicViscosity),false);
			directives.Add("permittivity",new tw::input::Numbers<tw::Float>(&eps[0],2),false);
		}
	};

	/// Allows packing of multiple materials in vector form
	struct material_set
	{
		std::valarray<tw::Float> mass,charge,cvm,excitationEnergy,thermo_cond_cvm,k_visc_m,eps_r,eps_i;
		void Allocate(tw::Int num)
		{
			mass.resize(num);
			charge.resize(num);
			cvm.resize(num);
			excitationEnergy.resize(num);
			thermo_cond_cvm.resize(num);
			k_visc_m.resize(num);
			eps_r.resize(num);
			eps_i.resize(num);
		}
		void AddMaterial(const material& mat,const tw::Int& i)
		{
			mass[i] = mat.mass;
			charge[i] = mat.charge;
			cvm[i] = mat.cvm;
			excitationEnergy[i] = mat.excitationEnergy;
			thermo_cond_cvm[i] = mat.thermometricConductivity * mat.cvm;
			k_visc_m[i] = mat.kinematicViscosity * mat.mass;
			eps_r[i] = mat.eps[0];
			eps_i[i] = mat.eps[1];
		}
	};

	tw::Float CoulombCrossSectionCGS(tw::Float q1,tw::Float q2,tw::Float m12,tw::Float v12,tw::Float N1,tw::Float N2,tw::Float T1,tw::Float T2)
	{
		// temperatures in ergs
		tw::Float rmin,rmax,rmin_alt,coulombLog;
		rmin = std::fabs(q1*q2)/(tw::small_pos + m12*v12*v12);
		rmin_alt = cgs::hbar/(tw::small_pos + 2*m12*v12);
		rmin = (rmin*rmin_alt)/(rmin+rmin_alt); // efficiently estimate the smaller of the two
		rmax = 1/std::sqrt(tw::small_pos + 4*pi*N1*q1*q1/T1 + 4*pi*N2*q2*q2/T2);
		coulombLog = std::log((3*rmin+rmax)/rmin); // well behaved approximation of std::log(rmax/rmin)
		return (32/pi)*std::pow(v12,-4)*sqr(q1*q2/m12)*coulombLog;
	}

	tw::Float ElectronPhononFrequencyCGS(tw::Float Ti,tw::Float EFermi,tw::Float ks)
	{
		// temperatures in ergs
		const tw::Float vF = std::sqrt(2*EFermi/cgs::me);
		// collision frequency in cgs units
		return 2*ks*sqr(cgs::qe)*Ti/(sqr(cgs::hbar)*vF);
	}
}
