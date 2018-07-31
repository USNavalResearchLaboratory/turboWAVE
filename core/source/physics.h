enum tw_ionization_model {noIonization,ADKTunneling,pptIonization,mpiSimple};

enum tw_dimensions
{
	time_dim,
	length_dim,
	number_dim,
	mass_dim,
	energy_dim,
	momentum_dim,
	angular_momentum_dim,
	density_dim,
	power_dim,
	fluence_dim,
	intensity_dim,
	energy_density_dim,
	power_density_dim,
	charge_dim,
	current_dim,
	current_density_dim,
	charge_density_dim,
	electric_field_dim,
	magnetic_field_dim,
	scalar_potential_dim,
	conductivity_dim,
	rate_coefficient_2_dim,
	rate_coefficient_3_dim,
	mobility_dim,
	temperature_dim,
	cross_section_dim,
	susceptibility_dim,
	susceptibility_2_dim,
	susceptibility_3_dim
};

struct UnitConverter
{
	tw::Float n1,wp;
	tw::Float c,qe,me,eps0,kB,hbar;

	UnitConverter(tw::Float unitDensityCGS);
	tw::Float MKSValue(tw_dimensions dim) const;
	tw::Float CGSValue(tw_dimensions dim) const;

	tw::Float SimToMKS(tw_dimensions dim,tw::Float val) const
	{
		return val*MKSValue(dim);
	}
	tw::Float MKSToSim(tw_dimensions dim,tw::Float val) const
	{
		return val/MKSValue(dim);
	}
	tw::Float SimToCGS(tw_dimensions dim,tw::Float val) const
	{
		return val*CGSValue(dim);
	}
	tw::Float CGSToSim(tw_dimensions dim,tw::Float val) const
	{
		return val/CGSValue(dim);
	}
	tw::Float sim_to_eV(tw::Float val) const
	{
		return val*MKSValue(energy_dim)/qe;
	}
	tw::Float eV_to_sim(tw::Float val) const
	{
		return val*qe/MKSValue(energy_dim);
	}
};

struct IonizationData
{
	tw::Float ionizationPotential;
	tw::Float electrons,protons;

	tw_ionization_model ionizationModel;
	tw::Float adkMultiplier,pptMultiplier;
	tw::Float C_PPT,C_ADK,C_ADK_AVG,nstar,lstar;
	tw::Float t_atomic,E_atomic,f_atomic_to_sim,E_sim_to_atomic;
	tw::Float photons,w0,E_MPI,max_rate;
	tw::Int terms;

	bool ionizeFromGas;
	tw::Int ionSpecies,electronSpecies;

	IonizationData();
	void Initialize(tw::Float unitDensity,tw::Float* carrierFrequency);
	tw::Float wfunc(tw::Float x);
	tw::Float Rate(tw::Float instant,tw::Float peak);
	void ReadInputFileTerm(std::stringstream& inputString,std::string& command);
};
