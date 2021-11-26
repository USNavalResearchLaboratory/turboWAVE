namespace mks
{
	static const tw::Float c=2.99792458e8,qe=1.6021766208e-19,me=9.10938356e-31,eps0=8.854187818e-12,kB=1.38064852e-23,hbar=1.0545718001e-34;
}
namespace cgs
{
	static const tw::Float c=2.99792458e10,qe=4.8032047139e-10,me=9.10938356e-28,kB=1.38064852e-16,hbar=1.0545718001e-27;
}
namespace tw
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

	class dnum
	{
		// Dimensional number class.  This allows:
		// 1. Reading of dimensional numbers from the input file
		// 2. Creating dimensional numbers from a string
		// 3. Conversion to ordinary floats in any supported system of units
		tw::Float prefix,value;
		tw::dims unit_dimension;
		tw::units unit_system;
	public:
		dnum();
		dnum(tw::Float v,const tw::dnum_units& d);
		dnum(const std::string& s);
		friend tw::Float operator >> (const tw::dnum& d,const tw::UnitConverter& uc);
		friend class tw::UnitConverter;
	};

	struct dnum_abstract
	{
		// This is an intermediate product needed during operator-based conversions.
		// It is a dimensional number prior to specifiying a particular system of units.
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

	class UnitConverter
	{
	public:
		tw::units native;
		tw::Float ne,wp;
		tw::Float c,qe,me,eps0,kB,hbar,alpha;

	public:
		UnitConverter();
		UnitConverter(tw::units sys,tw::Float unitDensityCGS);
		UnitConverter(tw::units sys,const tw::UnitConverter& uc);
		tw::Float FactorizedMKSValue(tw::dims dim,tw::Float m1,tw::Float w1,tw::Float l1,tw::Float u1,tw::Float q1,tw::Float a1,tw::Float T1) const;
		tw::Float MKSValue(tw::dims dim,tw::units sys) const;

		// The following two operators are typically cascaded, e.g., 1.0*tw::dims::length >> uc1 >> uc2;
		friend tw::dnum operator >> (const tw::dnum_abstract& a,const tw::UnitConverter& uc)
		{
			// Assign the native system of uc to the abstract dimensional number to produce a concrete dimensional number.
			// We do not use the assignment operator because it would (i) look strange and (ii) create issues of precedence.
			// First argument is typically an rvalue, e.g. (1.0*tw::dims::length)
			return tw::dnum(a.value,tw::dnum_units(a.unit_dimension,uc.native,1.0));
		}
		friend tw::Float operator >> (const tw::dnum& d,const tw::UnitConverter& uc)
		{
			// Convert the dimensional number d to an ordinary float in the native system of uc.
			// First argument is typically an rvalue formed using the other >> operator.
			return d.prefix*d.value*uc.MKSValue(d.unit_dimension,d.unit_system)/uc.MKSValue(d.unit_dimension,uc.native);
		}
	};

	inline std::map<std::string,tw::dnum_units> umap()
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

	inline std::string plasma_label(tw::dims dim)
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

	inline std::string mks_label(tw::dims dim)
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

	inline std::string cgs_label(tw::dims dim)
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

	inline std::string atomic_label(tw::dims dim)
	{
		if (dim==dims::none)
			return "None";
		else
			return "a.u.";
	}

	inline std::string natural_label(tw::dims dim)
	{
		if (dim==dims::none)
			return "None";
		else
			return "n.u.";
	}
}

std::istream& operator >> (std::istream& is,tw::dnum& d);
