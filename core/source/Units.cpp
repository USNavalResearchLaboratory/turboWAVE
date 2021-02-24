#include "meta_base.h"

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

tw::dnum::dnum(const std::string& s)
{
	std::istringstream temp(s);
	temp >> *this;
}

std::istream& tw::operator >> (std::istream& is,tw::dnum& d)
{
	// Function to read a dimensional number off the input stream.
	// Accept several forms, e.g., %-1.0[V], -1.0[V], -1.0 [V].
	// The '%' prefix is deprecated but still allowed, old specifiers no longer allowed.
	std::size_t endpos;
	std::string word,val,units;
	std::map<std::string,tw::dnum_units> umap = tw::umap();

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
		if (umap.find(units)==umap.end())
		{
			is.seekg(-units.size(),std::ios::cur);
			units = "";
		}
	}
	if (units=="")
	{
		d.unit_dimension = tw::dims::none;
		d.unit_system = tw::units::mks;
	}
	else
	{
		if (umap.find(units)==umap.end())
			throw tw::FatalError("Unrecognized Units " + units);
		if (umap.find(units)!=umap.end())
			d = dnum(d.value,umap[units]);
	}
	return is;
}

/////////////////////////////
//                         //
//  Unit Conversion Tools  //
//                         //
/////////////////////////////

tw::UnitConverter::UnitConverter()
{
	native = tw::units::plasma;
	// set up mks constants
	c = mks::c;
	qe = mks::qe;
	me = mks::me;
	eps0 = mks::eps0;
	kB = mks::kB;
	hbar = mks::hbar;
	alpha = qe*qe/(4*pi*eps0*hbar*c);
	// Setup plasma normalization to fail, i.e.,
	// we require an assignment from another constructor.
	ne = 0.0;
	wp = sqrt(ne*sqr(mks::qe)/(mks::eps0*mks::me));
}

tw::UnitConverter::UnitConverter(tw::units sys,tw::Float unitDensityCGS)
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

tw::UnitConverter::UnitConverter(tw::units sys,const tw::UnitConverter& uc)
{
	*this = uc;
	native = sys;
}

tw::Float tw::UnitConverter::FactorizedMKSValue(tw::dims dim,tw::Float m1,tw::Float w1,tw::Float l1,tw::Float u1,tw::Float q1,tw::Float a1,tw::Float T1) const
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
			return pow(l1,-3.0);
		case tw::dims::power:
			return u1*w1;
		case tw::dims::fluence:
			return u1/(l1*l1);
		case tw::dims::intensity:
			return u1*w1/(l1*l1);
		case tw::dims::mass_density:
			return m1*pow(l1,-3.0);
		case tw::dims::energy_density:
			return u1*pow(l1,-3.0);
		case tw::dims::power_density:
			return u1*w1*pow(l1,-3.0);
		case tw::dims::charge:
			return q1;
		case tw::dims::current:
			return q1*w1;
		case tw::dims::charge_density:
			return q1*pow(l1,-3.0);
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
			return w1*pow(l1,3.0);
		case tw::dims::rate_coefficient_3:
			return w1*pow(l1,6.0);
		case tw::dims::mobility:
			return v1/(u1/(q1*l1)); // v/E
		case tw::dims::temperature:
			return T1;
		case tw::dims::pressure:
			return u1*pow(l1,-3.0);
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

tw::Float tw::UnitConverter::MKSValue(tw::dims dim,tw::units sys) const
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
