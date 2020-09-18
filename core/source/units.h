namespace tw
{
	enum class units
	{
		mks,cgs,plasma,atomic,natural
	};
	inline std::map<std::string,tw::units> get_unit_map()
	{
		return {{"mks",tw::units::mks},{"cgs",tw::units::cgs},{"plasma",tw::units::plasma},{"atomic",tw::units::atomic},{"natural",tw::units::natural}};
	}
	enum class dimensions
	{
		none,
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
		cross_section
	};

	inline std::string plasma_label(tw::dimensions dim)
	{
		// These labels are inverse dimensions, meant to be concatenated with a dimensional variable.
		// E.g., the time label upon concatenation gives the label $\\omega t$.
		std::map<tw::dimensions,std::string> ans =
		{
			{dimensions::none,"None"},
			{dimensions::angular_frequency,"$/\\omega$"},
			{dimensions::frequency,"$/\\omega$"},
			{dimensions::time,"$\\omega$"},
			{dimensions::length,"$\\omega/c$"},
			{dimensions::velocity,"$/c$"},
			{dimensions::number,"$\\omega^3/n_cc^3$"},
			{dimensions::mass,"$/m_e$"},
			{dimensions::energy,"$/m_ec^2$"},
			{dimensions::momentum,"$/m_ec$"},
			{dimensions::angular_momentum,"$\\omega/m_e$"},
			{dimensions::density,"$/n_c$"},
			{dimensions::power,"$/m_ec^2\\omega$"},
			{dimensions::fluence,"$4\\pi e/m_e^2c^3\\omega$"},
			{dimensions::intensity,"$4\\pi e/m_e^2c^3\\omega^2$"},
			{dimensions::mass_density,"$/m_en_c$"},
			{dimensions::energy_density,"$/m_ec^2n_c$"},
			{dimensions::power_density,"$/m_ec^2n_c\\omega$"},
			{dimensions::charge,"$/e$"},
			{dimensions::current,"$\\omega^2/n_cec^3"},
			{dimensions::current_density,"$/n_cec$"},
			{dimensions::charge_density,"$/n_ce$"},
			{dimensions::electric_field,"$e/m_ec\\omega$"},
			{dimensions::magnetic_field,"$e/m_ec\\omega$"},
			{dimensions::scalar_potential,"$e/m_ec^2$"},
			{dimensions::vector_potential,"$e/m_ec^2$"},
			{dimensions::diffusivity,"$\\omega/c^2$"}, // appropriate if temperature is regarded as energy
			//{dimensions::thermal_conductivity,"$\\omega/n_c c^2$"}, // appropriate if temperature is regarded as energy
			{dimensions::thermal_conductivity,"$\\omega/k_B n_c c^2$"}, // appropriate if temperature is in degrees
			{dimensions::conductivity,"$m_e\\omega/n_ce^2$"},
			{dimensions::rate_coefficient_2,"$n_c/\\omega$"},
			{dimensions::rate_coefficient_3,"$n_c^2/\\omega$"},
			{dimensions::mobility,"$m_e\\omega/e$"},
			{dimensions::temperature,"$k_B/m_ec^2$"},
			{dimensions::pressure,"$/m_ec^2n_c$"},
			{dimensions::cross_section,"$n_cc/\\omega$"}
		};
		return ans[dim];
	}

	inline std::string mks_label(tw::dimensions dim)
	{
		std::map<tw::dimensions,std::string> ans =
		{
			{dimensions::none,"None"},
			{dimensions::angular_frequency,"rad$/$s"},
			{dimensions::frequency,"s$^{-1}$"},
			{dimensions::time,"s"},
			{dimensions::length,"m"},
			{dimensions::velocity,"m$/$s"},
			{dimensions::number,"particles"},
			{dimensions::mass,"kg"},
			{dimensions::energy,"J"},
			{dimensions::momentum,"kg$\\cdot$m$/$s"},
			{dimensions::angular_momentum,"kg$\\cdot$m$^2/$s"},
			{dimensions::density,"m$^{-3}$"},
			{dimensions::power,"W"},
			{dimensions::fluence,"J$/$m$^2$"},
			{dimensions::intensity,"W$/$m$^2$"},
			{dimensions::mass_density,"kg$/$m$^3$"},
			{dimensions::energy_density,"J$/$m$^3$"},
			{dimensions::power_density,"W$/$m$^3$"},
			{dimensions::charge,"C"},
			{dimensions::current,"A"},
			{dimensions::current_density,"A$/$m$^2$"},
			{dimensions::charge_density,"C$/$m$^3$"},
			{dimensions::electric_field,"V$/$m"},
			{dimensions::magnetic_field,"T"},
			{dimensions::scalar_potential,"V"},
			{dimensions::vector_potential,"weber$/$m"},
			{dimensions::diffusivity,"m$^2/$s"}, // if diffusing temperature in K we would multiply by J/K
			{dimensions::thermal_conductivity,"W$/$m$\\cdot$K"},
			{dimensions::conductivity,"S$/$m"},
			{dimensions::rate_coefficient_2,"m$^3/$s"},
			{dimensions::rate_coefficient_3,"m$^6/$s"},
			{dimensions::mobility,"m$^2/$V$\\cdot$s"},
			{dimensions::temperature,"K"},
			{dimensions::pressure,"Pa"},
			{dimensions::cross_section,"m$^2$"},
		};
		return ans[dim];
	}

	inline std::string cgs_label(tw::dimensions dim)
	{
		std::map<tw::dimensions,std::string> ans =
		{
			{dimensions::none,"None"},
			{dimensions::angular_frequency,"rad$/$s"},
			{dimensions::frequency,"s$^{-1}$"},
			{dimensions::time,"s"},
			{dimensions::length,"cm"},
			{dimensions::velocity,"cm$/$s"},
			{dimensions::number,"particles"},
			{dimensions::mass,"g"},
			{dimensions::energy,"erg"},
			{dimensions::momentum,"g$\\cdot$cm$/$s"},
			{dimensions::angular_momentum,"g$\\cdot$cm$^2/$s"},
			{dimensions::density,"cm$^{-3}$"},
			{dimensions::power,"erg$/$s"},
			{dimensions::fluence,"erg$/$cm$^2$"},
			{dimensions::intensity,"erg$/$s$\\cdot$cm$^2$"},
			{dimensions::mass_density,"g$/$cm$^3$"},
			{dimensions::energy_density,"erg$/$cm$^3$"},
			{dimensions::power_density,"erg$/$s$\\cdot$cm$^3$"},
			{dimensions::charge,"SC"},
			{dimensions::current,"SA"},
			{dimensions::current_density,"SA$/$cm$^2$"},
			{dimensions::charge_density,"SC$/$cm$^3$"},
			{dimensions::electric_field,"SV$/$cm"},
			{dimensions::magnetic_field,"G"},
			{dimensions::scalar_potential,"SV"},
			{dimensions::vector_potential,"G$\\cdot$cm"},
			{dimensions::diffusivity,"cm$^2/$s"}, // if diffusing temperature in K we would multiply by erg/K
			{dimensions::thermal_conductivity,"erg$/$cm$\\cdot$s$\\cdot$eV"},
			{dimensions::conductivity,"s$^{-1}$"},
			{dimensions::rate_coefficient_2,"cm$^3/$s"},
			{dimensions::rate_coefficient_3,"cm$^6/$s"},
			{dimensions::mobility,"cm$^2/$SV$\\cdot$s"},
			{dimensions::temperature,"eV"},
			{dimensions::pressure,"dyne$/$cm$^2$"},
			{dimensions::cross_section,"cm$^2$"},
		};
		return ans[dim];
	}

	inline std::string atomic_label(tw::dimensions dim)
	{
		if (dim==dimensions::none)
			return "None";
		else
			return "a.u.";
	}

	inline std::string natural_label(tw::dimensions dim)
	{
		if (dim==dimensions::none)
			return "None";
		else
			return "n.u.";
	}
}

struct UnitConverter
{
	tw::units native;
	tw::Float ne,wp;
	tw::Float c,qe,me,eps0,kB,hbar,alpha;

	UnitConverter(tw::units sys,tw::Float unitDensityCGS);
	tw::Float FactorizedMKSValue(tw::dimensions dim,tw::Float m1,tw::Float w1,tw::Float l1,tw::Float u1,tw::Float q1,tw::Float a1,tw::Float T1) const;
	tw::Float MKSValue(tw::dimensions dim,tw::units sys) const;
	tw::Float Convert(tw::Float val,tw::dimensions dim,tw::units convert_from,tw::units convert_to) const
	{
		return val*MKSValue(dim,convert_from)/MKSValue(dim,convert_to);
	}
	tw::Float ConvertToNative(tw::Float val,tw::dimensions dim,tw::units convert_from) const
	{
		return val*MKSValue(dim,convert_from)/MKSValue(dim,native);
	}
	tw::Float ConvertFromNative(tw::Float val,tw::dimensions dim,tw::units convert_to) const
	{
		return val*MKSValue(dim,native)/MKSValue(dim,convert_to);
	}
	[[deprecated]]
	tw::Float eV_to_sim(tw::Float val) const
	{
		return ConvertToNative(val*qe,tw::dimensions::energy,tw::units::mks);
	}
	[[deprecated]]
	tw::Float sim_to_eV(tw::Float val) const
	{
		return ConvertFromNative(val,tw::dimensions::energy,tw::units::mks)/qe;
	}
	[[deprecated]]
	tw::Float MKSToSim(tw::dimensions dim,tw::Float val) const
	{
		return ConvertToNative(val,dim,tw::units::mks);
	}
	[[deprecated]]
	tw::Float CGSToSim(tw::dimensions dim,tw::Float val) const
	{
		return ConvertToNative(val,dim,tw::units::cgs);
	}
	[[deprecated]]
	tw::Float SimToMKS(tw::dimensions dim,tw::Float val) const
	{
		return ConvertFromNative(val,dim,tw::units::mks);
	}
	[[deprecated]]
	tw::Float SimToCGS(tw::dimensions dim,tw::Float val) const
	{
		return ConvertFromNative(val,dim,tw::units::cgs);
	}
};
