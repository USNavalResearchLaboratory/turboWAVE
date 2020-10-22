#include "meta_base.h"

//////////////////////////
//                      //
//  Dimensional Number  //
//                      //
//////////////////////////

tw::dnum::dnum()
{
	prefix = 1.0;
	value = 0.0;
	unit_dimension = tw::dimensions::none;
	unit_system = tw::units::mks;
}

tw::dnum::dnum(tw::dimensions dim,tw::units sys,tw::Float pre)
{
	prefix = pre;
	value = 0.0;
	unit_dimension = dim;
	unit_system = sys;
}

tw::dnum::dnum(tw::Float v,const tw::dnum& d)
{
	*this = d;
	value = v;
}

std::istream& tw::operator >> (std::istream& is,tw::dnum& d)
{
	// Function to read a dimensional number off the input stream.
	// Accept several forms, e.g., -%1.0V, %-1.0[V], -1.0[V], -1.0 [V].
	// The '%' prefix and unbracketed unit are deprecated.
	std::size_t endpos;
	std::string word,val,units;
	std::map<std::string,tw::dnum> umap = tw::umap();
	std::map<std::string,tw::dnum> umap_alt = tw::umap_alt();

	is >> word;
	val = word;
	units = "";
	// If there is a '%' prefix remove it
	if (word.size()>1)
	{
		if (word[0]=='-' && word[1]=='%')
		{
			word[1] = '-';
			val = word.substr(1);
		}
		if (word[0]=='%')
		{
			val = word.substr(1);
		}
	}
	try { d.value = std::stod(val,&endpos); }
	catch (std::invalid_argument) { throw tw::FatalError("Invalid number : " + word); }
	if (endpos<val.size())
		units = val.substr(endpos);
	else
	{
		is >> units;
		if (umap.find(units)==umap.end() && umap_alt.find(units)==umap_alt.end())
		{
			is.seekg(-units.size(),std::ios::cur);
			units = "";
		}
	}
	if (units=="")
	{
		d.unit_dimension = tw::dimensions::none;
		d.unit_system = tw::units::mks;
	}
	else
	{
		if (umap.find(units)==umap.end() && umap_alt.find(units)==umap_alt.end())
			throw tw::FatalError("Unrecognized Units " + units);
		if (umap.find(units)!=umap.end())
			d = dnum(d.value,umap[units]);
		if (umap_alt.find(units)!=umap_alt.end())
			d = dnum(d.value,umap_alt[units]);
	}
	return is;
}

/////////////////////////////
//                         //
//  Unit Conversion Tools  //
//                         //
/////////////////////////////

UnitConverter::UnitConverter(tw::units sys,tw::Float unitDensityCGS)
{
	native = sys;
	// set up mks constants
	c = mks::c;
	qe = mks::qe;
	me = mks::me;
	eps0 = mks::eps0;
	kB = mks::kB;
	hbar = mks::hbar;
	alpha = qe*qe/(4*pi*eps0*hbar*c);
	// setup plasma normalization
	ne = unitDensityCGS*1e6;
	wp = sqrt(ne*sqr(mks::qe)/(mks::eps0*mks::me));
}

tw::Float UnitConverter::FactorizedMKSValue(tw::dimensions dim,tw::Float m1,tw::Float w1,tw::Float l1,tw::Float u1,tw::Float q1,tw::Float a1,tw::Float T1) const
{
	const tw::Float v1 = l1*w1;
	switch (dim)
	{
		case tw::dimensions::none:
			return 1.0;
		case tw::dimensions::angle:
			return 1.0;
		case tw::dimensions::angular_frequency:
			return w1;
		case tw::dimensions::frequency:
			return w1;
		case tw::dimensions::time:
			return 1.0/w1;
		case tw::dimensions::length:
			return l1;
		case tw::dimensions::velocity:
			return v1;
		case tw::dimensions::number:
			return 1.0;
		case tw::dimensions::mass:
			return m1;
		case tw::dimensions::energy:
			return u1;
		case tw::dimensions::momentum:
			return m1*v1;
		case tw::dimensions::angular_momentum:
			return m1*l1*v1;
		case tw::dimensions::density:
			return pow(l1,-3.0);
		case tw::dimensions::power:
			return u1*w1;
		case tw::dimensions::fluence:
			return u1/(l1*l1);
		case tw::dimensions::intensity:
			return u1*w1/(l1*l1);
		case tw::dimensions::mass_density:
			return m1*pow(l1,-3.0);
		case tw::dimensions::energy_density:
			return u1*pow(l1,-3.0);
		case tw::dimensions::power_density:
			return u1*w1*pow(l1,-3.0);
		case tw::dimensions::charge:
			return q1;
		case tw::dimensions::current:
			return q1*w1;
		case tw::dimensions::charge_density:
			return q1*pow(l1,-3.0);
		case tw::dimensions::current_density:
			return q1*w1/(l1*l1);
		case tw::dimensions::scalar_potential:
			return u1/q1;
		case tw::dimensions::vector_potential:
			return a1;
		case tw::dimensions::electric_field:
			return u1/(q1*l1);
		case tw::dimensions::magnetic_field:
			return a1/l1;
		case tw::dimensions::diffusivity: // If diffusing heat, this assumes temperature is in joules
			return v1*v1/w1;
		case tw::dimensions::thermal_conductivity:
			return u1*w1/(l1*T1);
		case tw::dimensions::conductivity:
			return (q1*w1/(l1*l1))/(u1/(q1*l1)); // j/E
		case tw::dimensions::rate_coefficient_2:
			return w1*pow(l1,3.0);
		case tw::dimensions::rate_coefficient_3:
			return w1*pow(l1,6.0);
		case tw::dimensions::mobility:
			return v1/(u1/(q1*l1)); // v/E
		case tw::dimensions::temperature:
			return T1;
		case tw::dimensions::pressure:
			return u1*pow(l1,-3.0);
		case tw::dimensions::cross_section:
			return l1*l1;
		default:
			throw tw::FatalError("Unit could not be converted ("+tw::mks_label(dim)+")");
			return 1.0;
	}
	return 1.0;
}

tw::Float UnitConverter::MKSValue(tw::dimensions dim,tw::units sys) const
{
	// This routine answers the question:
	// What is the mks value of the unit of <dim> in the <sys> system of units?
	// For example, if dim=length and sys=cgs, we get .01 meters.
	tw::Float m1,u1,l1,q1,w1,a1,T1;
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
			q1 = qe/sqrt(alpha);
			a1 = u1/(q1*l1*w1);
			T1 = u1/kB;
			return FactorizedMKSValue(dim,m1,w1,l1,u1,q1,a1,T1);
		case tw::units::plasma:
			switch (dim)
			{
				case tw::dimensions::none:
					return 1.0;
				case tw::dimensions::angle:
					return 1.0;
				case tw::dimensions::angular_frequency:
					return wp;
				case tw::dimensions::frequency:
					return wp;
				case tw::dimensions::time:
					return 1.0/wp;
				case tw::dimensions::length:
					return c/wp;
				case tw::dimensions::velocity:
					return c;
				case tw::dimensions::number:
					return ne*cub(c/wp);
				case tw::dimensions::mass:
					return me;
				case tw::dimensions::energy:
					return me*c*c;
				case tw::dimensions::momentum:
					return me*c;
				case tw::dimensions::angular_momentum:
					return me*c*c/wp;
				case tw::dimensions::density:
					return ne;
				case tw::dimensions::power:
					return me*c*c*wp;
				case tw::dimensions::fluence:
					return sqr(me*c*wp/qe)*c*eps0/wp; // E^2/eta0/wp
				case tw::dimensions::intensity:
					return sqr(me*c*wp/qe)*c*eps0; // E^2/eta0
				case tw::dimensions::mass_density:
					return me*ne;
				case tw::dimensions::energy_density:
					return me*c*c*ne;
				case tw::dimensions::power_density:
					return me*c*c*ne*wp;
				case tw::dimensions::charge:
					return qe;
				case tw::dimensions::current:
					return ne*qe*c*sqr(c/wp);
				case tw::dimensions::current_density:
					return ne*qe*c;
				case tw::dimensions::charge_density:
					return ne*qe;
				case tw::dimensions::electric_field:
					return me*c*wp/qe;
				case tw::dimensions::magnetic_field:
					return me*wp/qe;
				case tw::dimensions::scalar_potential:
					return me*c*c/qe;
				case tw::dimensions::vector_potential:
					return me*c/qe;
				case tw::dimensions::diffusivity:
					return c*c/wp;
				case tw::dimensions::thermal_conductivity:
					return c*c*ne*kB/wp;
				case tw::dimensions::conductivity:
					return (ne*qe*c)/(me*c*wp/qe); // j/E
				case tw::dimensions::rate_coefficient_2:
					return wp/ne;
				case tw::dimensions::rate_coefficient_3:
					return wp/(ne*ne);
				case tw::dimensions::mobility:
					return c/(me*c*wp/qe); // c/E
				case tw::dimensions::temperature:
					return me*c*c/kB;
				case tw::dimensions::pressure:
					return me*c*c*ne;
				case tw::dimensions::cross_section:
					return wp/(ne*c);
				default:
					throw tw::FatalError("Plasma unit could not be converted ("+tw::plasma_label(dim)+")");
					return 1.0;
			}
			return 1.0;
	}
	return 1.0;
}
