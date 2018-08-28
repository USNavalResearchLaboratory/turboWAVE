Input File: Fluids
==================

Relativistic Cold Fluid
-----------------------

The cold fluid module has priority 1 in the update sequence.

The following object generates a module that computes the motion of a cold, relativistic, electron fluid.  It is assumed there is an immobile background of ions.  Temperature only comes into the computation of the Coulomb collision frequency.  The electron temperature is controlled by the last ``generate`` block (see :ref:`matter-loading`).  There can only be one cold fluid module in a simulation.  Electron, ion, and neutral density are tracked in this one module.  Any profiles installed for the fluid refer to a gas with constant fractional ionization.

.. py:function:: new fluid name { directives }

	:param str name: assigns a name to the fluid
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`ionization`

		.. py:function:: charge = q

			:param float q: charge of the fluid constituent, usually electrons

		.. py:function:: mass = m

			:param float m: mass of the fluid constituent, usually electrons

		.. py:function initial ionization fraction = f

			:param float f: ratio of initial electron density to neutral + ion density

		.. py:function:: neutral cross section = sigma

		 	:param float sigma: electron-neutral collision cross section in :math:`cm^2`

		.. py:function:: coulomb collisions = cc

		 	:param bool cc: if true, add collision frequency based on Coulomb cross section to neutral collision frequency

		.. note::
			the ion species and electron species variables (see :ref:`ionization`) have no meaning here


SPARC Hydro Modules
-------------------

.. py:function:: new chemistry { directives }

	This is the top level SPARC module.  Internally it contains and manages all other SPARC modules.

	:param block directives: The following directives are supported:

		.. py:function:: epsilon factor = eps

			:param float eps: error tolerance for adaptive time step

		.. py:function:: radiation model = rad

		 	:param enum rad: takes values ``thin``, ``none``

		.. py:function:: laser model = las

		 	:param enum las: takes values ``vacuum``, ``isotropic``

		.. py:function:: plasma model = plas

			:param enum plas: takes values ``neutral``, ``quasineutral``

		.. py:function:: dipole center = Dx Dy Dz

		 	Reference point for dipole moment diagnostic

		.. py:function:: external potential = ( V1 , V2 )

			:param float V1: fixed potential at lower z boundary
			:param float V2: fixed potential at upper z boundary

		.. py:function:: tolerance = tol

		 	:param float tol: Elliptical solver error tolerance

		.. py:function:: overrelaxation = w

		 	:param float w: overrelaxation parameter

		.. py:function:: iterations = N

		 	:param int N: elliptical solver maximum iterations


.. _chemical:

.. py:function:: new chemical name { directives }

	:param str name: name given to the chemical species
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`ionization`, :ref:`eos`

		.. py:function:: mass = m0

			:param float m0: mass of the constituent particles, default = 1.0

		.. py:function:: charge = q0

			:param float q0: charge of the constituent particles, default = -1.0

		.. py:function:: cv = cv0

		 	:param float cv0: normalized specific heat at constant volume, :math:`mc_v/k_B`, a typical value is 1.5 for species with no internal degrees of freedom.

		.. py:function:: vibrational energy = epsv

		 	:param float epsv: energy (eV) between vibrational levels, default = 0 = no vibrations

		.. py:function:: implicit = tst

		 	:param bool tst: whether to use implicit electron advance for this species

		.. py:function:: thermometric conductivity = k

		 	:param float k: Thermometric conductivity in units of cm^2/s. Thermometric conductivity is :math:`K/\rho c_p`, where K = heat conductivity.  For air, k = 2e-5 m^2/s = 0.2 cm^2/s, and K = 2.5e-4 W/(cm*K). SPARC solves the heat equation :math:`\rho c_v dT/dt - \nabla\cdot (K \nabla T) = 0`.  For electrons the Braginskii conductivity is used.

		.. py:function:: kinematic viscosity = x

		 	:param float x: Kinematic viscosity in units of cm^2/s. Kinematic viscosity is :math:`X/\rho`, where X = dynamic viscosity. For air, kinematic viscosity is about 0.15 cm^2/s. SPARC solves the momentum diffusion equation :math:`\rho dv/dt - \nabla\cdot (X \nabla v) = 0`.

		.. py:function:: effective mass = meff

		 	:param float meff: effective mass for density >> 1.0 for electrons moving through this chemical

		.. py:function:: permittivity = (epsr,epsi)

		 	:param float epsr: real part of permittivity relative to free space permittivity
		 	:param float epsi: imaginary part of permittivity relative to free space permittivity


.. py:function:: new group { directives }

	Create an equilibrium group module.  This is a container for chemical species that are assumed to be in equilibrium with one another, and therefore have a common temperature and velocity.  All chemicals are part of a group.  If a chemical is declared outside any group, one is automatically created.

	:param block directives: The following directives are supported:

		.. py:function:: new chemical name { directives }

			see :ref:`chemical <chemical>` for description of parameters.  Can be repeated to associate multiple chemicals with the group.

		.. py:function:: mobile = tst

			:param bool tst: whether chemicals in this group are mobile or immobile


SPARC Collision Directives
--------------------------

SPARC collisions broadly include elastic and inelastic collisions, as well as chemical reactions.  All such processes have to explicitly resolved.

.. py:function:: new collision = sp1 <-> sp2 , cross section = sigma

	:param str sp1: name of chemical species 1 in two-body collision
	:param str sp2: name of chemical species 2 in two-body collision
	:param float sigma: cross section normalized to :math:`\omega_p/n_1c`

.. py:function:: new collision = sp1 <-> sp1 , coulomb

	Uses the Coulomb collision cross section, derived from local conditions.

	:param str sp1: name of chemical species 1 in two-body collision
	:param str sp2: name of chemical species 2 in two-body collision

.. py:function:: new collision = sp1 <-> sp1 , metallic , ks = ks0 , fermi_energy_ev = ef , ref_density = nref

	Uses the harmonic mean of electron-phonon and coulomb collision rates

	:param float ks0: some lattice constant, see K. Eidmann et al., Phys. Rev. E 62, 1202 (2000)
	:param float ef: the Fermi energy in eV
	:param float nref: the density at which the formula directly applies

.. py:function:: new reaction = { eq1 : eq2 : eq3 : ... } rate = c0 c1 c2 cat(range)

	Sets up a chemical reaction between arbitrary species using a modified Arrhenius rate

	:math:`{\cal R} = c_0 T^{c_1} \exp(-c_2/T)`

	over a range of temperatures.  Piecewise rate constructions can be created by using multiple reactions which have the same equation but different rates and different temperature ranges.

	:param str eq1: chemical equation, or subreaction, in the form ``r1 + r2 + ... -> p1 + p2 + ... + eps``, where ``r1`` etc. are replaced by names of reactants, ``p1`` etc. are replaced by names of products, and ``eps`` is a heat of reaction in eV.  Breaking the reaction into subreactions can be used to control the flow of energy from reactants to products.

	:param float c0: rate coefficient in cm^(3(N-1))/s, where N is the number of reactants
	:param float c1: dimensionless exponent in rate law
	:param float c2: temperature reference appearing in rate law in eV
	:param str cat: name of the chemical to be considered the catalyst, i.e., the one whose temperature affects the rate
	:param numpy_range range: range of temperatures specified as in numpy, i.e., T1:T2, where T1 and T2 are floating point literals, given in eV.  Also as in numpy, :T2 means 0-T2, while T1: means T1-infinity.

.. py:function:: new reaction = { eq1 : eq2 : eq3 : ... } janev_rate = c0 c1 c2 c3 c4 c5 c6 c7 c8 cat(range)

	Alternative way of specifying the reaction rate:

		:math:`\ln {\cal R} = \sum_{n=0}^{8} c_n (\ln T)^n`

.. py:function:: new excitation = sp1 -> sp2 level = n rate = c0 c1 c2

	Vibrational excitation of one species by another.  If level = n the transition is between ground and level n.  If level = 0 the transition is between adjacent levels, where it is assumed the rate for transitions from n to n+1 is the same for all n.

.. py:function:: new excitation = sp1 -> sp2 level = n janev_rate = c0 c1 c2 c3 c4 c5 c6 c7 c8

	Alternative way of specifying the excitation rate.
