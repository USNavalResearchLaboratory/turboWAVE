
module;

#include "tw_includes.h"

/**
 * @brief Module for handling equation of state (EOS)
 *
 * @details
 * Although this module's native units can be in any supported system, in practice
 * it is still necessary to use plasma units to interact correctly with other modules.
 * Literal values are often coded in cgs and then converted to native units.
 * The user can use any preferred units in the input file.
 * 
 * ### Variables and units specific to EOS
 * * E - specific internal energy in ergs/g
 * * nmE - internal energy density in ergs/cm^3
 * * nm - mass density in g/cm^3
 * * T - temperature in eV
 * * Tv - vibrational temperature in eV
 * * P - pressure in dynes/cm^2
 * * K - heat conductivity in ergs/eV/cm/s
 * * cvm - specific heat times characteristic mass in ergs/eV
 * * nmcv - specific heat times mass density in ergs/eV/cm^3
 *
 * ### Variables and units from Hydro module
 * * n - particle density in 1/cm^3
 * * u - total energy density in ergs/cm^3
 * * x - vibrational energy density in ergs/cm^3
 * * np - momentum density in g/cm^2/s
 */
export module eos;
export import :eos_component;
export import :eos_mix;
