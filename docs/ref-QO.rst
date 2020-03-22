Input File: Quantum
===================


.. _qstate:

Quantum State Tools
-------------------

.. _qstate-shared:

Quantum State Shared Directives
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

The following may be used in any quantum state tool.

		.. py:function:: amplitude = ( re , im )

			:param float re: real part of amplitude to use for relative scaling and phasing
			:param float im: imaginary part of amplitude to use for relative scaling and phasing

		.. py:function:: cylindrical = tst

		 	:param bool tst: if true use bound states appropriate for cylindrical atoms.

Specific Quantum State Tools
,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. py:function:: new qstate free [<name>] [for <module>] { <directives> }

	Create a free state.

	:param string module: Name of the module that will use the quantum state.
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`qstate-shared`

		.. py:function:: k4 = ( E , kx , ky , kz )

			The four-momentum of the free state.

			:param float kx: x-component of momentum
			:param float ky: y-component of momentum
			:param float kz: z-component of momentum
			:param float E: Energy of the free state, only the sign is used.  The magnitude of the energy is always computed from the momentum.  Ignored for non-relativistic equations.

		.. py:function:: spin = ( sx , sy , sz )

			Define the orientation of the spin.  The magnitude of the vector is ignored.  Spin 1/2 is always assumed.

			:param float sx: x-component of the spin direction
			:param float sy: y-component of the spin direction
			:param float sz: z-component of the spin direction

		.. py:function:: size = ( Lx , Ly , Lz )

			Size of the wave packet envelope.

.. py:function:: new qstate random [<name>] [for <module>] { <directives> }

	Create a random state.

	:param string module: Name of the module that will use the quantum state.
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`qstate-shared`

		.. py:function:: size = ( Lx , Ly , Lz )

			Size of the wave packet envelope.

.. py:function:: new qstate bound [<name>] [for <module>] { <directives> }

	Create a bound state.

	:param string module: Name of the module that will use the quantum state.
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`qstate-shared`

		.. py:function:: nr_j_l_m = ( nr, j, l, m )

			Set quantum numbers defining a bound state, see :doc:`bak-quantum` for full discussion.

		 	:param int nr: radial quantum number
			:param float j: total angular momentum quantum number
			:param int l: parity quantum number
			:param float m: angular momentum projection

			.. tip::
				The principle quantum number from Schroedinger theory is

				:math:`n = n_r + l + 1`

.. py:function:: new qstate tabulated [<name>] [for <module>] { <directives> }

	Create a state using data from a file.

	:param string module: Name of the module that will use the quantum state.
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`qstate-shared`

		.. py:function:: filename = fname`

			Name of the file containing the data describing the state.  The format is given :ref:`here <state-file>`.

Quantum Modules
---------------

.. _quantum-shared:

Quantum Module Shared Directives
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

The following may be used in any quantum propagation module.

.. py:function:: orbiting charge = q0

	:param float q0: Charge of the particle in the external potential, in units appropriate for the given module.  See :doc:`bak-quantum` regarding units.

.. py:function:: orbiting mass = m0

	:param float m0: mass of the particle in the external potential, in units of electronic mass

.. py:function:: soft core potential , charge = Q , radius = dr

	The soft core potential is given by

	:math:`\Phi = \frac{Q}{\sqrt{r^2 + \delta r^2}}`

	:param float Q: charge associated with the soft core potential
	:param float dr: radius associated with the soft core potential

Specific Quantum Propagation Modules
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. py:function:: new schroedinger equation module { directives }

	Creates a module for solving the time dependent Schroedinger equation for a particle in an arbitrary external field.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`quantum-shared`

		Installable tools: :ref:`qstate`, :ref:`radiation`

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

		Installable tools: :ref:`qstate`, :ref:`radiation`


.. py:function:: new dirac equation module { directives }

	Creates a module for solving the time dependent Dirac equation for a particle in an arbitrary external field.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`quantum-shared`
		
		Installable tools: :ref:`qstate`, :ref:`radiation`

Quantum Diagnostics
-------------------

There is a diagnostic module for performing overlap integrals against reference states.  This is useful for tracking occupation probabilities.

.. py:function:: new population diagnostic <name> { <directives> }

	Creates the quantum population diagnostic.  This diagnostic module uses an ``energy diagnostic`` tool to report the real and imaginary part of an overlap integral as a function of time.

	:param string name: The name of the diagnostic
	:param block directives: The directives block may contain declarations of quantum state tools or energy diagnostic tools.  These tools can also be attached using any little language syntax.

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
