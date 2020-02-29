#include "meta_base.h"

/////////////////////////////
//                         //
//  Unit Conversion Tools  //
//                         //
/////////////////////////////

UnitConverter::UnitConverter(tw::Float unitDensityCGS)
{
	n1 = unitDensityCGS*1e6;
	wp = sqrt(n1*sqr(mks::qe)/(mks::eps0*mks::me));
	c = mks::c;
	qe = mks::qe;
	me = mks::me;
	eps0 = mks::eps0;
	kB = mks::kB;
	hbar = mks::hbar;
	alpha = qe*qe/(4*pi*eps0*hbar*c);
}

tw::Float UnitConverter::AtomicValue(tw_dimensions dim) const
{
	switch (dim)
	{
		case angular_frequency_dim:
			return hbar*wp/(me*sqr(alpha*c));
		case time_dim:
			return me*sqr(alpha*c)/(hbar*wp);
		case length_dim:
			return me*alpha*c*c/(hbar*wp);
		case mass_dim:
			return 1.0;
		case energy_dim:
			return 1.0/alpha/alpha;
		case momentum_dim:
			return 1.0/alpha;
		case charge_dim:
			return 1.0;
		case electric_field_dim:
			return hbar*wp/(me*c*c*cub(alpha));
		case scalar_potential_dim:
			return 1.0/alpha/alpha;
		default:
			return 0.0;
	}
}

tw::Float UnitConverter::MKSValue(tw_dimensions dim) const
{
	switch (dim)
	{
		case angular_frequency_dim:
			return wp;
		case time_dim:
			return 1.0/wp;
		case length_dim:
			return c/wp;
		case number_dim:
			return n1*cub(c/wp);
		case mass_dim:
			return me;
		case energy_dim:
			return me*c*c;
		case momentum_dim:
			return me*c;
		case angular_momentum_dim:
			return me*c*c/wp;
		case density_dim:
			return n1;
		case power_dim:
			return me*c*c*wp;
		case fluence_dim:
			return sqr(me*c*wp/qe)*c*eps0/wp; // E^2/eta0/wp
		case intensity_dim:
			return sqr(me*c*wp/qe)*c*eps0; // E^2/eta0
		case energy_density_dim:
			return me*c*c*n1;
		case power_density_dim:
			return me*c*c*n1*wp;
		case charge_dim:
			return qe;
		case current_dim:
			return n1*qe*c*sqr(c/wp);
		case current_density_dim:
			return n1*qe*c;
		case charge_density_dim:
			return n1*qe;
		case electric_field_dim:
			return me*c*wp/qe;
		case magnetic_field_dim:
			return me*wp/qe;
		case scalar_potential_dim:
			return me*c*c/qe;
		case diffusivity_dim:
			return c*c/wp;
		case conductivity_dim:
			return (n1*qe*c)/(me*c*wp/qe); // j/E
		case rate_coefficient_2_dim:
			return wp/n1;
		case rate_coefficient_3_dim:
			return wp/(n1*n1);
		case mobility_dim:
			return c/(me*c*wp/qe); // c/E
		case temperature_dim:
			return me*c*c/kB;
		case cross_section_dim:
			return wp/(n1*c);
		case susceptibility_dim:
			return 1.0;
		case susceptibility_2_dim:
			return 1.0/(me*c*wp/qe); // 1/E
		case susceptibility_3_dim:
			return 1.0/sqr(me*c*wp/qe); // 1/E^2
		default:
			return 0.0;
	}
}

tw::Float UnitConverter::CGSValue(tw_dimensions dim) const
{
	switch (dim)
	{
		case angular_frequency_dim:
			return wp;
		case time_dim:
			return 1.0/wp;
		case length_dim:
			return 1e2*c/wp;
		case number_dim:
			return n1*cub(c/wp);
		case mass_dim:
			return 1e3*me;
		case energy_dim:
			return 1e7*me*c*c;
		case momentum_dim:
			return 1e5*me*c;
		case angular_momentum_dim:
			return 1e7*me*c*c/wp;
		case density_dim:
			return 1e-6*n1;
		case power_dim:
			return 1e7*me*c*c*wp;
		case fluence_dim:
			return 1e3*sqr(me*c*wp/qe)*c*eps0/wp; // E^2/eta0/wp
		case intensity_dim:
			return 1e3*sqr(me*c*wp/qe)*c*eps0; // E^2/eta0
		case energy_density_dim:
			return 1e1*me*c*c*n1;
		case power_density_dim:
			return 1e1*me*c*c*n1*wp;
		case charge_dim:
			return 3e9*qe;
		case current_dim:
			return 3e9*n1*qe*c*sqr(c/wp);
		case current_density_dim:
			return 3e5*n1*qe*c;
		case charge_density_dim:
			return 3e3*n1*qe;
		case electric_field_dim:
			return 0.333333e-4*me*c*wp/qe;
		case magnetic_field_dim:
			return 1e4*me*wp/qe;
		case scalar_potential_dim:
			return 0.333333e-2*me*c*c/qe;
		case diffusivity_dim:
			return 1e4*c*c/wp;
		case conductivity_dim:
			return 9e9*(n1*qe*c)/(me*c*wp/qe); // j/E
		case rate_coefficient_2_dim:
			return 1e6*wp/n1;
		case rate_coefficient_3_dim:
			return 1e12*wp/(n1*n1);
		case mobility_dim:
			return 1e2*c/(0.333333e-4*me*c*wp/qe); // c/E
		case temperature_dim:
			return me*c*c/kB;
		case cross_section_dim:
			return 1e4*wp/(n1*c);
		case susceptibility_dim:
			return 1.0/(4.0*pi);
		case susceptibility_2_dim:
			return 1.0/(4.0*pi*0.333333e-4*me*c*wp/qe); // 1/4*pi*E
		case susceptibility_3_dim:
			return 1.0/sqr(4.0*pi*0.333333e-4*me*c*wp/qe); // 1/4*pi*E^2
		default:
			return 0.0;
	}
}
