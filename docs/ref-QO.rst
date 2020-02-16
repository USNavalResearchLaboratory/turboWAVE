Input File: Quantum
===================

.. _quantum-shared:

Quantum Shared Directives
-------------------------

The following may be used in any quantum optics module.

.. py:function:: orbiting charge = q0

	:param float q0: Charge of the particle in the external potential, in units appropriate for the given module.  See :doc:`bak-quantum` regarding units.

.. py:function:: orbiting mass = m0

	:param float m0: mass of the particle in the external potential, in units of electronic mass

.. py:function:: soft core potential , charge = Q , radius = dr

	The soft core potential is given by

	:math:`\Phi = \frac{Q}{\sqrt{r^2 + \delta r^2}}`

	:param float Q: charge associated with the soft core potential
	:param float dr: radius associated with the soft core potential

.. py:function:: bachelet potential = c1 , c2 , a1 , a2 , 0 0 0 0 0 0 0 0 0

	This is a more elaborate model for the effective potential in the single electron approximation.


.. _wavefunction:

.. py:function:: new wavefunction { directives }

	Create an initial state.  Can be repeated to form a superposition state.

	:param block directives: The following directives are supported:

		.. py:function:: type = wv_type

			:param enum wv_type: selected from the following:

				:samp:`lookup` --- read the wavefunction from a file.  See :ref:`state-file`.

				:samp:`free` --- use a wave packet that is nearly a definite momentum state, with definite spin in the rest frame.

				:samp:`helicity` --- use a wave packet that is nearly a definite momentum state, with definite helicity.

				:samp:`bound` --- use a stationary state solution assuming a Coulomb potential.

				:samp:`random` --- use a random wavefunction with a given maximum amplitude and spatial envelope.

		.. py:function:: file = fname

			:param str fname: name of the file to use if :samp:`type = lookup`

		.. py:function:: amplitude = ( re , im )

			:param float re: real part of amplitude to use for relative scaling and phasing
			:param float im: imaginary part of amplitude to use for relative scaling and phasing

		.. py:function:: k_e_s = ( kx , ky , kz , E , sz )

			Set quantum numbers defining a free state.

			:param float kx: x-component of momentum
			:param float ky: y-component of momentum
			:param float kz: z-component of momentum
			:param float E: Energy of the free state, only the sign is used.  The magnitude of the energy is always computed from the momentum.  Ignored for non-relativistic equations.
			:param float sz: If :samp:`type=helicity`, this sets the spin projected onto the momentum axis, which must be either -0.5 or 0.5.  If :samp:`type=free`, this sets the z-component of the spin in the rest frame, which must also be either -0.5 or 0.5.  Ignored for scalar equations.

		.. py:function:: size = ( sx , sy , sz )

			Determines the size of the wave packet envelope in the cases :samp:`type=free`, :samp:`type=helicity`, and :samp:`type=random`.

		.. py:function:: cylindrical = tst

		 	:param bool tst: if true use bound states appropriate for cylindrical atoms.

		.. py:function:: nr_j_l_m = nr j l m

			Set quantum numbers defining a bound state, see :doc:`bak-quantum` for full discussion.

		 	:param int nr: radial quantum number
			:param float j: total angular momentum quantum number
			:param int l: parity quantum number
			:param float m: angular momentum projection

			.. tip::
				The principle quantum number from Schroedinger theory is

				:math:`n = n_r + l + 1`


.. py:function:: new reference { directives }

	Create a wavefunction to use as a reference state.  Directives are exactly as in :ref:`wavefunction <wavefunction>`.



Quantum Propagation Modules
---------------------------

.. py:function:: new schroedinger equation module { directives }

	Creates a module for solving the time dependent Schroedinger equation for a particle in an arbitrary external field.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`quantum-shared`

		.. py:function:: keep a2 term = tst

			:param bool tst: whether to keep the second order term from the Hamiltonian

		.. py:function:: dipole approximation = tst

		 	:param bool tst: if true, vector potential is uniform (always evaluated at origin)

		.. py:function:: relaxation time = tr

		 	:param float tr: Causes module to spend this much time relaxing to ground. This may help refine the initial condition, but can be omitted. If used, one may start with a random wavefunction in order to not prejudice the results.


.. py:function:: new klein gordon equation module { directives }

	Creates a module for solving the time dependent Klein-Gordon equation for a particle in an arbitrary external field.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`quantum-shared`

		.. py:function:: dipole approximation = tst

		 	:param bool tst: if true, vector potential is uniform (always evaluated at origin)


.. py:function:: new dirac equation module { directives }

	Creates a module for solving the time dependent Dirac equation for a particle in an arbitrary external field.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`quantum-shared`



Bohmian Trajectories
--------------------

Particle species defined as in :doc:`ref-PIC` can be used to model Bohmian trajectories that are guided by the quantum propagation modules.  All matter loading directives are available, see :ref:`matter-loading`.  In Bohmian mechanics, particles move based on the guidance condition:

	:math:`{\bf v} = {\bf j}/\rho`

where :math:`{\bf j}` is the probability current appropriate for the wave equation in question, and :math:`\rho` is the probability density.  These must satisfy the conservation law

	:math:`\partial \rho / \partial t + \nabla \cdot {\bf j} = 0`

In order to recover the statistical predictions of conventional quantum mechanics, the Bohmian particles should be loaded into a particle density that is commensurate with the probability density of the wavefunction.  TurboWAVE will do this if the matter loading parameters satisfy :samp:`loading = statistical` and :samp:`particle weight = fixed`.  Otherwise, the initial Bohmian density and wavefunction probability density will be treated as independent.


.. _state-file:

State File Format
------------------

The state file used to initialize bound states that are computed externally is an ASCII text file.  All white space is treated as equivalent.  White space before and after "=" is required.

.. highlight:: none

The file format is::

	energy = [#]
	pts = [#]
	components = [#]
	cell_width = [#]
	soft_core_radius = [#]
	nuclear_charge = [#]
	nr_j_l_m = [#] [#] [#] [#]
	Bz = [#]
	cylindrical = [tw_bool]

	start_data
	[#]
	[#]

The ``tw_bool`` type resolves to true if the text is ``true``, ``yes``, or ``on``, false otherwise.
Any information can be added before ``start_data`` provided there are no label collisions.
The data that follows start data is an alternating list of real and imaginary parts
of the radial function at each radial grid position.  Angular functions are implicit
in the quantum numbers.
