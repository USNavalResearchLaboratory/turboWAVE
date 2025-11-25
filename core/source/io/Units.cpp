module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module units;
import base;
import navigate;

export namespace mks
{
	constexpr tw::Float c=2.99792458e8,qe=1.6021766208e-19,me=9.10938356e-31,eps0=8.854187818e-12,kB=1.38064852e-23,hbar=1.0545718001e-34;
	constexpr tw::Float alpha = qe*qe/(4*pi*eps0*hbar*c);

}
export namespace cgs
{
	constexpr tw::Float c=2.99792458e10,qe=4.8032047139e-10,me=9.10938356e-28,kB=1.38064852e-16,hbar=1.0545718001e-27;
	constexpr tw::Float alpha = qe*qe/(hbar*c);
}
export namespace tw
{
	class UnitConverter;

	enum class units
	{
		mks,cgs,plasma,atomic,natural
	};
	inline std::map<std::string,tw::units> get_unit_map()
	{
		return {{"mks",tw::units::mks},{"cgs",tw::units::cgs},{"plasma",tw::units::plasma},{"atomic",tw::units::atomic},{"natural",tw::units::natural}};
	}
	inline std::map<tw::units,std::string> get_unit_map_r()
	{
		return {{tw::units::mks,"mks"},{tw::units::cgs,"cgs"},{tw::units::plasma,"plasma"},{tw::units::atomic,"atomic"},{tw::units::natural,"natural"}};
	}
	enum class dims
	{
		none,
		angle,
		action,
		angular_frequency,
		frequency,
		time,
		length,
		velocity,
		number,
		mass,
		energy,
		momentum,
		angular_momentum,
		density,
		power,
		fluence,
		intensity,
		mass_density,
		energy_density,
		power_density,
		charge,
		current,
		current_density,
		charge_density,
		electric_field,
		magnetic_field,
		scalar_potential,
		vector_potential,
		diffusivity,
		thermal_conductivity,
		conductivity,
		rate_coefficient_2,
		rate_coefficient_3,
		mobility,
		temperature,
		pressure,
		specific_energy,
		cross_section
	};
	struct dnum_units
	{
		// Units of a dimensional number.
		// This is used to represent information in a dimensional post-fix such as "[W/cm2]"
		// This is used mainly as an argument to a dnum constructor.
		tw::dims unit_dimension;
		tw::units unit_system;
		tw::Float prefix;
		dnum_units() { unit_dimension=tw::dims::none; unit_system=tw::units::mks; prefix=1.0; }
		dnum_units(tw::dims dim,tw::units sys,tw::Float pre) { unit_dimension=dim; unit_system=sys; prefix=pre; }
	};
	std::map<std::string,tw::dnum_units> umap()
	{
		return
		{
			// time and space
			{"[rad]",tw::dnum_units(tw::dims::angle,tw::units::mks,1.0)},
			{"[mrad]",tw::dnum_units(tw::dims::angle,tw::units::mks,1e-3)},
			{"[urad]",tw::dnum_units(tw::dims::angle,tw::units::mks,1e-6)},
			{"[deg]",tw::dnum_units(tw::dims::angle,tw::units::mks,pi/180.0)},
			{"[um]",tw::dnum_units(tw::dims::length,tw::units::mks,1e-6)},
			{"[mm]",tw::dnum_units(tw::dims::length,tw::units::mks,1e-3)},
			{"[cm]",tw::dnum_units(tw::dims::length,tw::units::mks,1e-2)},
			{"[m]",tw::dnum_units(tw::dims::length,tw::units::mks,1.0)},
			{"[fs]",tw::dnum_units(tw::dims::time,tw::units::mks,1e-15)},
			{"[ps]",tw::dnum_units(tw::dims::time,tw::units::mks,1e-12)},
			{"[ns]",tw::dnum_units(tw::dims::time,tw::units::mks,1e-9)},
			{"[us]",tw::dnum_units(tw::dims::time,tw::units::mks,1e-6)},
			{"[s]",tw::dnum_units(tw::dims::time,tw::units::mks,1.0)},
			// thermodynamics
			{"[/m3]",tw::dnum_units(tw::dims::density,tw::units::mks,1.0)},
			{"[/cm3]",tw::dnum_units(tw::dims::density,tw::units::cgs,1.0)},
			{"[kg/m3]",tw::dnum_units(tw::dims::mass_density,tw::units::mks,1.0)},
			{"[g/cm3]",tw::dnum_units(tw::dims::mass_density,tw::units::cgs,1.0)},
			{"[J/m3]",tw::dnum_units(tw::dims::energy_density,tw::units::mks,1.0)},
			{"[J/cm3]",tw::dnum_units(tw::dims::energy_density,tw::units::mks,1e6)},
			{"[eV]",tw::dnum_units(tw::dims::temperature,tw::units::cgs,1.0)},
			{"[K]",tw::dnum_units(tw::dims::temperature,tw::units::mks,1.0)},
			{"[Pa]",tw::dnum_units(tw::dims::pressure,tw::units::mks,1.0)},
			{"[dynes/cm2]",tw::dnum_units(tw::dims::pressure,tw::units::cgs,1.0)},
			{"[bar]",tw::dnum_units(tw::dims::pressure,tw::units::mks,1e5)},
			{"[J/kg]",tw::dnum_units(tw::dims::specific_energy,tw::units::mks,1.0)},
			{"[ergs/g]",tw::dnum_units(tw::dims::specific_energy,tw::units::cgs,1.0)},
			// transport
			{"[cm2]",tw::dnum_units(tw::dims::cross_section,tw::units::cgs,1.0)},
			{"[m2]",tw::dnum_units(tw::dims::cross_section,tw::units::mks,1.0)},
			{"[cm2/s]",tw::dnum_units(tw::dims::diffusivity,tw::units::cgs,1.0)},
			{"[m2/s]",tw::dnum_units(tw::dims::diffusivity,tw::units::mks,1.0)},
			// elecrodynamics
			{"[V]",tw::dnum_units(tw::dims::scalar_potential,tw::units::mks,1.0)},
			{"[webers/m]",tw::dnum_units(tw::dims::vector_potential,tw::units::mks,1.0)},
			{"[G*cm]",tw::dnum_units(tw::dims::vector_potential,tw::units::cgs,1.0)},
			{"[V/m]",tw::dnum_units(tw::dims::electric_field,tw::units::mks,1.0)},
			{"[V/cm]",tw::dnum_units(tw::dims::electric_field,tw::units::mks,100.0)},
			{"[T]",tw::dnum_units(tw::dims::magnetic_field,tw::units::mks,1.0)},
			{"[G]",tw::dnum_units(tw::dims::magnetic_field,tw::units::cgs,1.0)}
		};
	}

	std::string plasma_label(tw::dims dim)
	{
		// These labels are inverse dimensions, meant to be concatenated with a dimensional variable.
		// E.g., the time label upon concatenation gives the label $\\omega t$.
		std::map<tw::dims,std::string> ans =
		{
			{dims::none,"None"},
			{dims::angle,"None"},
			{dims::action,"$\\omega/m_ec^2$"},
			{dims::angular_frequency,"$/\\omega$"},
			{dims::frequency,"$/\\omega$"},
			{dims::time,"$\\omega$"},
			{dims::length,"$\\omega/c$"},
			{dims::velocity,"$/c$"},
			{dims::number,"$\\omega^3/n_cc^3$"},
			{dims::mass,"$/m_e$"},
			{dims::energy,"$/m_ec^2$"},
			{dims::momentum,"$/m_ec$"},
			{dims::angular_momentum,"$\\omega/m_e$"},
			{dims::density,"$/n_c$"},
			{dims::power,"$/m_ec^2\\omega$"},
			{dims::fluence,"$4\\pi e/m_e^2c^3\\omega$"},
			{dims::intensity,"$4\\pi e/m_e^2c^3\\omega^2$"},
			{dims::mass_density,"$/m_en_c$"},
			{dims::energy_density,"$/m_ec^2n_c$"},
			{dims::power_density,"$/m_ec^2n_c\\omega$"},
			{dims::charge,"$/e$"},
			{dims::current,"$\\omega^2/n_cec^3"},
			{dims::current_density,"$/n_cec$"},
			{dims::charge_density,"$/n_ce$"},
			{dims::electric_field,"$e/m_ec\\omega$"},
			{dims::magnetic_field,"$e/m_ec\\omega$"},
			{dims::scalar_potential,"$e/m_ec^2$"},
			{dims::vector_potential,"$e/m_ec^2$"},
			{dims::diffusivity,"$\\omega/c^2$"}, // appropriate if temperature is regarded as energy
			//{dims::thermal_conductivity,"$\\omega/n_c c^2$"}, // appropriate if temperature is regarded as energy
			{dims::thermal_conductivity,"$\\omega/k_B n_c c^2$"}, // appropriate if temperature is in degrees
			{dims::conductivity,"$m_e\\omega/n_ce^2$"},
			{dims::rate_coefficient_2,"$n_c/\\omega$"},
			{dims::rate_coefficient_3,"$n_c^2/\\omega$"},
			{dims::mobility,"$m_e\\omega/e$"},
			{dims::temperature,"$k_B/m_ec^2$"},
			{dims::pressure,"$/m_ec^2n_c$"},
			{dims::specific_energy,"$/c^2$"},
			{dims::cross_section,"$n_cc/\\omega$"}
		};
		return ans[dim];
	}

	std::string mks_label(tw::dims dim)
	{
		std::map<tw::dims,std::string> ans =
		{
			{dims::none,"None"},
			{dims::angle,"rad"},
			{dims::action,"J$\\cdot$s"},
			{dims::angular_frequency,"rad$/$s"},
			{dims::frequency,"s$^{-1}$"},
			{dims::time,"s"},
			{dims::length,"m"},
			{dims::velocity,"m$/$s"},
			{dims::number,"particles"},
			{dims::mass,"kg"},
			{dims::energy,"J"},
			{dims::momentum,"kg$\\cdot$m$/$s"},
			{dims::angular_momentum,"kg$\\cdot$m$^2/$s"},
			{dims::density,"m$^{-3}$"},
			{dims::power,"W"},
			{dims::fluence,"J$/$m$^2$"},
			{dims::intensity,"W$/$m$^2$"},
			{dims::mass_density,"kg$/$m$^3$"},
			{dims::energy_density,"J$/$m$^3$"},
			{dims::power_density,"W$/$m$^3$"},
			{dims::charge,"C"},
			{dims::current,"A"},
			{dims::current_density,"A$/$m$^2$"},
			{dims::charge_density,"C$/$m$^3$"},
			{dims::electric_field,"V$/$m"},
			{dims::magnetic_field,"T"},
			{dims::scalar_potential,"V"},
			{dims::vector_potential,"weber$/$m"},
			{dims::diffusivity,"m$^2/$s"}, // if diffusing temperature in K we would multiply by J/K
			{dims::thermal_conductivity,"W$/$m$\\cdot$K"},
			{dims::conductivity,"S$/$m"},
			{dims::rate_coefficient_2,"m$^3/$s"},
			{dims::rate_coefficient_3,"m$^6/$s"},
			{dims::mobility,"m$^2/$V$\\cdot$s"},
			{dims::temperature,"K"},
			{dims::pressure,"Pa"},
			{dims::specific_energy,"J$/$kg"},
			{dims::cross_section,"m$^2$"},
		};
		return ans[dim];
	}

	std::string cgs_label(tw::dims dim)
	{
		std::map<tw::dims,std::string> ans =
		{
			{dims::none,"None"},
			{dims::angle,"rad"},
			{dims::action,"erg$\\cdot$s"},
			{dims::angular_frequency,"rad$/$s"},
			{dims::frequency,"s$^{-1}$"},
			{dims::time,"s"},
			{dims::length,"cm"},
			{dims::velocity,"cm$/$s"},
			{dims::number,"particles"},
			{dims::mass,"g"},
			{dims::energy,"erg"},
			{dims::momentum,"g$\\cdot$cm$/$s"},
			{dims::angular_momentum,"g$\\cdot$cm$^2/$s"},
			{dims::density,"cm$^{-3}$"},
			{dims::power,"erg$/$s"},
			{dims::fluence,"erg$/$cm$^2$"},
			{dims::intensity,"erg$/$s$\\cdot$cm$^2$"},
			{dims::mass_density,"g$/$cm$^3$"},
			{dims::energy_density,"erg$/$cm$^3$"},
			{dims::power_density,"erg$/$s$\\cdot$cm$^3$"},
			{dims::charge,"SC"},
			{dims::current,"SA"},
			{dims::current_density,"SA$/$cm$^2$"},
			{dims::charge_density,"SC$/$cm$^3$"},
			{dims::electric_field,"SV$/$cm"},
			{dims::magnetic_field,"G"},
			{dims::scalar_potential,"SV"},
			{dims::vector_potential,"G$\\cdot$cm"},
			{dims::diffusivity,"cm$^2/$s"}, // if diffusing temperature in K we would multiply by erg/K
			{dims::thermal_conductivity,"erg$/$cm$\\cdot$s$\\cdot$eV"},
			{dims::conductivity,"s$^{-1}$"},
			{dims::rate_coefficient_2,"cm$^3/$s"},
			{dims::rate_coefficient_3,"cm$^6/$s"},
			{dims::mobility,"cm$^2/$SV$\\cdot$s"},
			{dims::temperature,"eV"},
			{dims::pressure,"dyne$/$cm$^2$"},
			{dims::specific_energy,"ergs$/$g"},
			{dims::cross_section,"cm$^2$"},
		};
		return ans[dim];
	}

	std::string atomic_label(tw::dims dim)
	{
		if (dim==dims::none)
			return "None";
		else
			return "a.u.";
	}

	std::string natural_label(tw::dims dim)
	{
		if (dim==dims::none)
			return "None";
		else
			return "n.u.";
	}
}

tw::Float FactorizedMKSValue(tw::dims dim,tw::Float m1,tw::Float w1,tw::Float l1,tw::Float u1,tw::Float q1,tw::Float a1,tw::Float T1)
{
	const tw::Float v1 = l1*w1;
	switch (dim)
	{
		case tw::dims::none:
			return 1.0;
		case tw::dims::angle:
			return 1.0;
		case tw::dims::angular_frequency:
			return w1;
		case tw::dims::frequency:
			return w1;
		case tw::dims::time:
			return 1.0/w1;
		case tw::dims::length:
			return l1;
		case tw::dims::velocity:
			return v1;
		case tw::dims::number:
			return 1.0;
		case tw::dims::mass:
			return m1;
		case tw::dims::energy:
			return u1;
		case tw::dims::momentum:
			return m1*v1;
		case tw::dims::angular_momentum:
			return m1*l1*v1;
		case tw::dims::density:
			return std::pow(l1,-3.0);
		case tw::dims::power:
			return u1*w1;
		case tw::dims::fluence:
			return u1/(l1*l1);
		case tw::dims::intensity:
			return u1*w1/(l1*l1);
		case tw::dims::mass_density:
			return m1*std::pow(l1,-3.0);
		case tw::dims::energy_density:
			return u1*std::pow(l1,-3.0);
		case tw::dims::power_density:
			return u1*w1*std::pow(l1,-3.0);
		case tw::dims::charge:
			return q1;
		case tw::dims::current:
			return q1*w1;
		case tw::dims::charge_density:
			return q1*std::pow(l1,-3.0);
		case tw::dims::current_density:
			return q1*w1/(l1*l1);
		case tw::dims::scalar_potential:
			return u1/q1;
		case tw::dims::vector_potential:
			return a1;
		case tw::dims::electric_field:
			return u1/(q1*l1);
		case tw::dims::magnetic_field:
			return a1/l1;
		case tw::dims::diffusivity: // If diffusing heat, this assumes temperature is in joules
			return v1*v1/w1;
		case tw::dims::thermal_conductivity:
			return u1*w1/(l1*T1);
		case tw::dims::conductivity:
			return (q1*w1/(l1*l1))/(u1/(q1*l1)); // j/E
		case tw::dims::rate_coefficient_2:
			return w1*std::pow(l1,3.0);
		case tw::dims::rate_coefficient_3:
			return w1*std::pow(l1,6.0);
		case tw::dims::mobility:
			return v1/(u1/(q1*l1)); // v/E
		case tw::dims::temperature:
			return T1;
		case tw::dims::pressure:
			return u1*std::pow(l1,-3.0);
		case tw::dims::specific_energy:
			return u1/m1;
		case tw::dims::cross_section:
			return l1*l1;
		default:
			throw tw::FatalError("Unit could not be converted ("+tw::mks_label(dim)+")");
			return 1.0;
	}
	return 1.0;
}

/// This routine answers the question:
/// What is the mks value of the unit of <dim> in the <sys> system of units?
/// For example, if dim=length and sys=cgs, we get .01 meters.
tw::Float MKSValue(tw::dims dim,tw::units sys,tw::Float wp,tw::Float ne)
{
	tw::Float m1,u1,l1,q1,w1,a1,T1;
	tw::Float c = mks::c;
	tw::Float qe = mks::qe;
	tw::Float me = mks::me;
	tw::Float alpha = mks::alpha;
	tw::Float kB = mks::kB;
	tw::Float hbar = mks::hbar;
	tw::Float eps0 = mks::eps0;
	switch (sys)
	{
		case tw::units::mks:
			return 1.0;
		case tw::units::cgs:
			m1 = .001;
			u1 = 1e-7;
			l1 = 0.01;
			w1 = 1.0;
			q1 = 0.1/c;
			a1 = 1e-6;
			T1 = qe/kB;
			return FactorizedMKSValue(dim,m1,w1,l1,u1,q1,a1,T1);
		case tw::units::atomic:
			m1 = me;
			u1 = me*sqr(alpha*c);
			l1 = hbar/(me*alpha*c);
			w1 = u1/hbar;
			q1 = qe;
			a1 = u1/(q1*l1*w1);
			T1 = u1/kB;
			return FactorizedMKSValue(dim,m1,w1,l1,u1,q1,a1,T1);
		case tw::units::natural:
			m1 = me;
			u1 = me*c*c;
			l1 = hbar/(me*c);
			w1 = u1/hbar;
			q1 = qe/std::sqrt(alpha);
			a1 = u1/(q1*l1*w1);
			T1 = u1/kB;
			return FactorizedMKSValue(dim,m1,w1,l1,u1,q1,a1,T1);
		case tw::units::plasma:
			switch (dim)
			{
				case tw::dims::none:
					return 1.0;
				case tw::dims::angle:
					return 1.0;
				case tw::dims::angular_frequency:
					return wp;
				case tw::dims::frequency:
					return wp;
				case tw::dims::time:
					return 1.0/wp;
				case tw::dims::length:
					return c/wp;
				case tw::dims::velocity:
					return c;
				case tw::dims::number:
					return ne*cub(c/wp);
				case tw::dims::mass:
					return me;
				case tw::dims::energy:
					return me*c*c;
				case tw::dims::momentum:
					return me*c;
				case tw::dims::angular_momentum:
					return me*c*c/wp;
				case tw::dims::density:
					return ne;
				case tw::dims::power:
					return me*c*c*wp;
				case tw::dims::fluence:
					return sqr(me*c*wp/qe)*c*eps0/wp; // E^2/eta0/wp
				case tw::dims::intensity:
					return sqr(me*c*wp/qe)*c*eps0; // E^2/eta0
				case tw::dims::mass_density:
					return me*ne;
				case tw::dims::energy_density:
					return me*c*c*ne;
				case tw::dims::power_density:
					return me*c*c*ne*wp;
				case tw::dims::charge:
					return qe;
				case tw::dims::current:
					return ne*qe*c*sqr(c/wp);
				case tw::dims::current_density:
					return ne*qe*c;
				case tw::dims::charge_density:
					return ne*qe;
				case tw::dims::electric_field:
					return me*c*wp/qe;
				case tw::dims::magnetic_field:
					return me*wp/qe;
				case tw::dims::scalar_potential:
					return me*c*c/qe;
				case tw::dims::vector_potential:
					return me*c/qe;
				case tw::dims::diffusivity:
					return c*c/wp;
				case tw::dims::thermal_conductivity:
					return c*c*ne*kB/wp;
				case tw::dims::conductivity:
					return (ne*qe*c)/(me*c*wp/qe); // j/E
				case tw::dims::rate_coefficient_2:
					return wp/ne;
				case tw::dims::rate_coefficient_3:
					return wp/(ne*ne);
				case tw::dims::mobility:
					return c/(me*c*wp/qe); // c/E
				case tw::dims::temperature:
					return me*c*c/kB;
				case tw::dims::pressure:
					return me*c*c*ne;
				case tw::dims::specific_energy:
					return c*c;
				case tw::dims::cross_section:
					return wp/(ne*c);
				default:
					throw tw::FatalError("Plasma unit could not be converted ("+tw::plasma_label(dim)+")");
					return 1.0;
			}
			return 1.0;
	}
	return 1.0;
}

export namespace tw {
	/// Dimensional number class, often used in a `>>` operator pipeline.  This allows:
	/// 1. Reading of dimensional numbers from the input file
	/// 2. Creating dimensional numbers from a string
	/// 3. Conversion to ordinary floats in any supported system of units
	class dnum
	{
		tw::Float prefix,value;
		tw::dims unit_dimension;
		tw::units unit_system;
	public:
		dnum();
		dnum(const std::string& src);
		dnum(tw::Float v,const tw::dnum_units& d);
		dnum(TSTreeCursor *curs,const std::string &src);
		friend tw::Float operator >> (const tw::dnum& d,const tw::UnitConverter& uc);
		friend class tw::UnitConverter;
	};

	/// This is an intermediate product needed during `>>`-operator conversions.
	/// It is a number with a physical interpretation, but prior to attaching units.
	struct dnum_abstract
	{
		tw::Float value;
		tw::dims unit_dimension;
		dnum_abstract() { value=0.0; unit_dimension=tw::dims::none; }
		dnum_abstract(tw::Float v,tw::dims dims) { value=v; unit_dimension=dims; }
		friend dnum operator >> (const tw::dnum_abstract& a,const tw::UnitConverter& uc);
		friend class tw::UnitConverter;
	};

	inline tw::dnum_abstract operator * (const tw::Float& v,const tw::dims& dims)
	{
		return dnum_abstract(v,dims);
	}

	/// This class is mainly used as a node in a `>>`-operator pipeline.
	/// Each UnitConverter object is associated with a specific system of units.
	class UnitConverter
	{
	public:
		tw::units unit_system;
		tw::Float ne,wp;

	public:
		UnitConverter() {
			unit_system = tw::units::plasma;
			// Setup plasma normalization to fail, i.e.,
			// we require an assignment from another constructor.
			ne = 0.0;
			wp = std::sqrt(ne*sqr(mks::qe)/(mks::eps0*mks::me));
		}
		UnitConverter(tw::units sys,tw::Float unitDensityCGS) {
			unit_system = sys;
			// setup plasma normalization
			ne = unitDensityCGS*1e6;
			wp = std::sqrt(ne*sqr(mks::qe)/(mks::eps0*mks::me));
		}
		tw::Float UnitDensityCGS() const {
			return ne*1e-6;
		}
		/// This simply attaches units to a raw number and an abstract dimension, e.g.
		/// `1.0*tw::dims::length >> uc` would give a `dnum` representing 1 cm, if uc were a cgs converter.
		/// Here, the first product produces a `dnum_abstract` object.
		friend tw::dnum operator >> (const tw::dnum_abstract& a,const tw::UnitConverter& uc)
		{
			return tw::dnum(a.value,tw::dnum_units(a.unit_dimension,uc.unit_system,1.0));
		}
		/// Convert dimensional number to ordinary float using the units associated with `uc`.
		/// This is often cascaded with the other `>>` operator, i.e.,
		/// `1.0*tw::dims::length >> cgs_converter >> native_converter` would give 1 cm in
		/// whatever units are associated with `native_converter`.
		friend tw::Float operator >> (const tw::dnum& d,const tw::UnitConverter& uc)
		{
			return d.prefix*d.value*MKSValue(d.unit_dimension,d.unit_system,uc.wp,uc.ne)/
				MKSValue(d.unit_dimension,uc.unit_system,uc.wp,uc.ne);
		}
	};
}

std::istream& operator >> (std::istream& is,tw::dnum& d);

//////////////////////////
//                      //
//  Dimensional Number  //
//                      //
//////////////////////////

tw::dnum::dnum()
{
	value = 0.0;
	prefix = 1.0;
	unit_dimension = tw::dims::none;
	unit_system = tw::units::mks;
}

tw::dnum::dnum(tw::Float v,const tw::dnum_units& t)
{
	value = v;
	prefix = t.prefix;
	unit_dimension = t.unit_dimension;
	unit_system = t.unit_system;
}

/// @brief parse dimensional number at cursor, also works on bare decimal node
/// @param curs cursor, should be on either dimension node or decimal node
/// @param src the text of the source document
tw::dnum::dnum(TSTreeCursor *curs,const std::string& src) {
	std::map<std::string,tw::dnum_units> umap = tw::umap();
	if (tw::input::node_kind(curs) == "decimal") {
		value = std::stod(input::node_text(curs,src));
		prefix = 1.0;
		unit_dimension = tw::dims::none;
		unit_system = tw::units::mks;
	} else if (tw::input::node_kind(curs) == "dimension") {
		if (ts_tree_cursor_goto_first_child(curs)) {
			value = std::stod(input::node_text(curs,src));
			if (ts_tree_cursor_goto_next_sibling(curs)) {
				tw::dnum_units units = umap[input::node_text(curs,src)];
				prefix = units.prefix;
				unit_dimension = units.unit_dimension;
				unit_system = units.unit_system;
			}
			ts_tree_cursor_goto_parent(curs);
		} else {
			input::ThrowParsingError(curs,src,"something missing");
		}
	} else {
		input::ThrowParsingError(curs,src,"expected number");
	}
}

tw::dnum::dnum(const std::string& src) {
	std::map<std::string,tw::dnum_units> umap = tw::umap();
	auto s = src.find('[');
	if (s==std::string::npos) {
		value = std::stod(src);
		prefix = 1.0;
		unit_dimension = tw::dims::none;
		unit_system = tw::units::mks;
		return;
	}
	value = std::stod(src);
	tw::dnum_units units = umap[src.substr(s)];
	prefix = units.prefix;
	unit_dimension = units.unit_dimension;
	unit_system = units.unit_system;
}
