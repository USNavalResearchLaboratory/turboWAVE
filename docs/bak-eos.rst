Equations of State (EOS)
=========================

The Equation of State module (EOS) suite is designed to calculate themodynamic quantities (such as pressure, :math:`P`, and Temperature, :math:`T` ) from hydrodynamically calculated, thermodynamically meaningful quantities (such as energy density, :math:`u`, and particle density, :math:`n`). In other words, each EOS module takes in an array of hydrodynamic values and returns a corresponding array of thermodynamic values, :math:`P(n,u)` and :math:`T(n,u)`. The functional relationship between all of these quantities is called an Equation of State. Each EOS class within the turbowave structure represents a different material model, or source data, from which these quantities are calculated. The model used may be specified by the user at input.

It is important to mention that there are, broadly speaking, two types of EOS relations that need to be specified to describe the material completely; The 'Thermal' EOS, which specifies the pressure :math:`P(n,T)`, and the 'Caloric' EOS, which specifies the energy :math:`U(n,T)` both together define a material's complete EOS. In many applications the caloric EOS is implicitly defined by a constant specific heat, :math:`C_v`, and the relation :math:`dU = C_v dT`, in which case only the thermal EOS is explicitly specified.

The thermodynamic quantities :math:`P` and :math:`T` are integral in closing the set of the fluid dynamic equations describing a material - particularly the momentum density equation and the energy density equation. It is important that in a hydro simulation these quantities are calculated from a self-consistent model or source data.

The Principle Hugoniot
----------------------

The Rankine-Hugoniot relation is discussed here because of its relevance to EOS models that use a point along the Hugoniot as a reference point for a broader calculation of EOS. They describe the relationship between the states on both sides of a abrupt density and momentum discontinuity (shock wave) traveling through a material at a shock speed :math:`D`. The collection of possible final states resulting from these relations, given zero initial particle speed and pressure of a material, is called the principle Hugoniot. There are an infinite number of possible Hugoniots for a material depending on its initial conditions, but there is only one principle Hugoniot for any given material, and the principle Hugoniot is, physically speaking, determined by the EOS of the material. For this reason experiments are conducted which measure the shock velocity vs. fluid velocity for various materials as a means to better understand their equations of state.

The Hugoniot relations are fundamentally an expression of the convervation laws for mass, momentum, and energy. Let's suppose that a shock front, or a planar discontinuity in momentum (also pressure), energy, and density, moves along a fluid at speed :math:`D` across a material. The pressure and density of the material prior to the passing of the shock front is :math:`P_0` and :math:`\rho_0`, respectively. Immediately after the shock passes the pressure and density is :math:`P` and :math:`\rho`. Similarly, :math:`V` and :math:`V_0` here are the final and initial specific volumes, respectively, and :math:`E` and :math:`E_0` here are the final and initial specific internal energies.

A detailed derivation will not be presented here (an example of one can be found in Ref. [1]), but the priciple behind it is simple. By setting the quantity of mass entering the shock as equal to that which is leaving it (conservation of mass) we can argue

:math:`\rho ( D - u ) = \rho_0 ( D - u_0 ).`

We can express the change in momentum density across the shock front as caused by impulse applied by the pressure difference. This expression turns out to be

:math:`P - P_0 = \rho_0 (D - u_0) (u - u_0).`

The pressure, together with the velocity terms form an expression for the total amount of work done by the medium over time, and an expression for the change in specific internal energy emerges such that

:math:`E - E_0 = \frac{1}{2} (P + P_0) (V - V_0).`

Experimentally speaking, data for :math:`D` and :math:`u` are collected by inducing a shock on a substance and measuring the outcome (and one may suppose :math:`P_0 = 0` and :math:`u_0 = 0` so that we are referring to the principle Hugoniot). The density, pressure, and energy are calculated with the above equations and the EOS values along the shock front is found. In practice, for most solids the relationship between the shock velocity and particle velocity is found to be linear, and approximated by

:math:`D = c_0 + S_1 u.`

:math:`c_0` here is the y-intercept of the Hugoniot data, and is usually close to but not always the same as the normal sound speed of the material. :math:`S_1` is the slope of the linear interpolation of the Hugoniot for that given set of data. 

Mie Gruneisen EOS Theory
-------------------------

The simplest EOS implemented in Turbowave, other than the default ideal gas model, is the Mie-Gruneisen EOS, which is often used to describe a shock-compressed solid. It is an example of an incomplete EOS, as it only provides the pressure, and not the temperature, as a function of the specific internal energy. A basic assumption of the Mie Gruneisen model is that the thermal energy of a material is adequately described as the sum of the energies of a collection of simple harmonic oscillators with frequencies :math:`\nu_i(V)` are functions of specific volume only. It also neglects the electronic contributions to the total internal energy :math:`E` of the material. 

A detailed derivation of this EOS is given in Ref. [1], and here we give a very basic outline. Given the physical description detailed above the Helmholtz free energy for the system is

:math:`F = \phi(V) + \sum_{i=1}^{3 N} \frac{h \nu_i}{2} + k T \sum_{i=1}^{3 N}\ln(1 - e^{-h \nu_i/kT}),`

where :math:`\phi(V)` is the potential energy of the material with :math:`N` atoms in a total volume, :math:`V`. The sum is over the :math:`3 N` normal modes of the atoms. We can determine an expression for pressure by :math:`P = - ( \partial F/\partial V)_T`, which results in various derivatives of :math:`\nu_i (V)`. We define

:math:`\gamma_i \equiv -\frac{d \ln v_i}{d \ln V} = -\frac{V}{\nu_i} \left( \frac{d v_i}{d V} \right),`

and substitute all the derivatives of :math:`\nu_i` with expressions of :math:`\gamma_i`. If we assume that all :math:`\gamma_i` is the same ( or :math:`\gamma_i = \gamma_G` ), :math:`\gamma_G` falls outside the sum over :math:`i`, leaving an expression that can be expressed as a function of internal vibrational energy rather than a sum of frequencies. If we understand the specific internal energy to be a sum of the vibrational energy and the potential energy (:math:`E = E_\text{vib} + \phi(V)`), we result in an expression that can relate a difference in pressure with a difference in the associated internal energies,

:math:`P - P_R = \frac{\gamma_G}{V} (E - E_R),`

where the reference internal energy, :math:`E_R`, can be on any reference curve. This expression is useful for calculating states off the Hugoniot relative to measured states on the Hugoniot. The quantity :math:`\gamma_G` is called the Gruneisen parameter, and expresses how the energy changes with respect to internal energy at constant volume, or

:math:`\gamma_G = V \left( \frac{\partial P}{\partial E} \right)_V.`

The numerical value of :math:`\gamma_G` for a given material at a given density can be experimentally measured and may be looked up from available tables and databases (See for example Refs. [2] and [3]). 

Broadly speaking, the Gruneisen parameter describes the effect that changing the volume (and spacing) of a crystal lattice has on its vibrational proterties, and subsequently, on its pressure/size dynamics.

Mie Gruneisen EOS Implementation
---------------------------------

Due to varying assumptions implicit in the functional dependence of :math:`\gamma` on volume, choice of reference pressure, :math:`P_R`, and energy, :math:`E_r`, and the fit of the Hugoniot that may be used for reference, several different implementations of the Mie-Gruneisen pressure calculation are possible. Here we present two different Mie Gruneisen pressure calculations currently available on Turbowave, alongside their associated input parameters. Both of these models are currently implemented as a purely thermal EOS. In other words, the assumption :math:`dU = C_v dT` is made in which :math:`C_v` is a constant scalar value specified at the input. The models differ in the pressure calculation only.

		.. py:function:: EOS = mie-grunseisen

			As a very rough implementation of the mie-gruneisen pressure law, we might treat :math:`\gamma_G` as a constant value specified at the input. Values for :math:`\gamma_G` can be looked up on tables for specific reference densities, and is approximately correct if the density does not significantly deviate from the reference density. In addition, this minimal implementation simply lets :math:`E_R = 0` and :math:`P_R = 0` rather than referring to a experimentally known reference point or Hugoniot curve. It qualitatively describes Mie Gruneisen - like pressure behavior, but will not in most cases quantitatively recreate known physical results, such as the sound speed for the material. 

		.. py:function:: EOS = mie-grunseisen2

			A better approximation that is usually good is to take the Gruneisen coefficient as proportional to the specific volume ( :math:`\rho \gamma_G = \text{Const.}` ). This implementation uses this assumption, and in addition, implicitly derives reference energies and pressure from a linear Hugoniot fit. As a result, two additional parameters are needed to be specified; The y-intercept of the Hugoniot, :math:`c_0`, and the slope of the linear fit, :math:`S_1` must be given in the input deck in addition to the reference Gruneisen parameter. These may be acquired from published Hugoniot measurements (See Refs. [2] and [3]).

In principle other, more complex implementations of the Mie Gruneisen EOS may be added in the future, such as potentiall the one in Ref. [1], which uses a cubic interpolation of the Hugoniot and allows for a small linear deviation from the assumption that ( :math:`\rho \gamma_G = \text{Const.}` ). These seem to not be necessary for typical solids, for which the linear interpolation is sufficient. In addition, a different caloric EOS may be implimented in combination with these models. This background will be extended as needed to illustrate the underlying differences in such models. 


References
-----------
[1] Gathers, R. G., "Selected Topics in Shock Wave Physics and Equation of State Modeling", World Scientific (1994).

[2] Marsh, S. P., ed., "LASL Shock Hugoniot Data", University of California Press (1980)

[3] McQueen, R. G., Marsh, S. P, *Equation of State for Nineteen Elements from Shock-Wave Measurements to Two Megabars*, J. Appl. Phys **31**, 1253-1269 (1960)


