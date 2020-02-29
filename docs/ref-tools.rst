Input File: Tools
=================

General Information
-------------------

TurboWAVE objects come in two flavors, Modules and Tools.  Modules are the high level objects that orchestrate a simulation.  Tools are lower level objects dedicated to specific computations. Tools are intended to be attached to one or more modules (modules can share the same tool).  Modules can be placed in a hierarchy.  Tools can only be attached to a module.

Tools are created like any object, see :doc:`ref-input`.

Many modules are able to create their tools automatically.  This page emphasizes those that tend to benefit from explicit management by the user.

Elliptic Solvers
----------------

All elliptic solvers share the following directives:

.. py:function:: poisson boundary condition <coord> = ( <bc1> , <bc2> )

	:param enum coord: can be ``x``, ``y``, ``z``
	:param enum bc1: boundary condition on lower side, can be ``open``, ``dirichlet``, ``neumann``.
	:param enum bc2: boundary condition on lower side, can be ``open``, ``dirichlet``, ``neumann``.

Iterative Solver
,,,,,,,,,,,,,,,,

.. py:function:: new iterative elliptic [<optional keys>] [for] <name> { directives }

	Uses successive over-relaxation to iteratively solve the elliptic equation.  This solver is slow, but flexible.  There is no limit on the topology of the boundary conditions, and arbitrary coordinates are supported.  The following directives are supported:

		Shared directives: see base elliptic solver

		.. py:function:: tolerance = <tol>

			:param float tol: iterate until the residual is reduced to this level

		.. py:function:: overrelaxation = <ov>

			:param float ov: overrides the default overrelaxation parameter (not generally recommended)

FACR Solver
,,,,,,,,,,,,,,,,

.. py:function:: new facr elliptic [<optional keys>] [for] <name> { directives }

	Uses Fourier analysis is in the transverse directions.  This solver is fast, but boundary conditions can only be imposed on constant z-surfaces, and Cartesian coordinates are required.  The following directives are supported:

		Shared directives: see base elliptic solver

Eigenmode Solver
,,,,,,,,,,,,,,,,

.. py:function:: new eigenmode elliptic [<optional keys>] [for] <name> { directives }

	Uses generalized spectral resolution of the transverse coordinates.  This solver works in arbitrary coordinates, and is fast as long as the transverse modes are truncated.  Boundary conditions can only be imposed on constant z-surfaces.  The following directives are supported:

		Shared directives: see base elliptic solver

		.. py:function:: modes = <N>

			:param int N: The number of transverse modes to keep.  The modes are taken from an ordered list, sorted by magnitude of the eigenvalue.

Laser Propagator
----------------

Eigenmode Propagator
,,,,,,,,,,,,,,,,,,,,

.. py:function:: new eigenmode propagator [<optional keys>] [for] <name> { directives }

	Uses generalized spectral resolution of the transverse coordinates.  This propagator works in arbitrary coordinates, and is fast as long as the transverse modes are truncated.  It has superior fidelity for highly dispersive systems.  The following directives are supported:

	.. py:function:: modes = <n>

		:param int n: maximum number of radial modes to keep (eigenmode propagator only)

	.. py:function:: damping time = <t>

		:param float t: e-folding time in the absorbing layers

	.. py:function:: absorbing layers = <l>

		:param int l: number of absorbing layers

ADI Propagator
,,,,,,,,,,,,,,,,,,,,

.. py:function:: new adi propagator [<optional keys>] [for] <name> { directives }

	Uses alternating direction implicit method.  This is a fast propagator that works in arbitrary coordinates.  It has poor fidelity for highly dispersive systems.  There are no directives.

PhotoIonization
---------------

Ionization Shared Directives
,,,,,,,,,,,,,,,,,,,,,,,,,,,,

All the photoionization tools support the following directives:

.. py:function:: ionization potential = ip

	:param float ip: Ionization potential, units are specified as usual, e.g., ``ionization potential = %13.6eV``

.. py:function:: saturated rate = sr

 	:param float sr: saturate the ionization rate at this value

.. py:function:: protons = np

 	:param int np: number of protons in nucleus (not needed for mpi model ; currently used to form residual charge only)

.. py:function:: electrons = ne

 	:param int ne: number of bound electrons (not needed for mpi model ; currently used to form residual charge only)

.. py:function:: ion species = is_name

	:param str is_name: name of a species to add a particle to upon ionization (usually positive charge)

.. py:function:: electron species = es_name

	:param str es_name: name of a species to add a particle to upon ionization (usually negative charge)

Multi-photon Ionization
,,,,,,,,,,,,,,,,,,,,,,,,

Model appropriate for low fields or high frequencies.

.. py:function:: new mpi ionization [<optional keys>] [for] <name> { directives }

	The following directives are supported:

		Shared directives: see above

		.. py:function:: reference field = E0

		 	:param float E0: :math:`E_0`, where the MPI rate is proportional to :math:`(E/E_0)^{2l}`

ADK Tunneling Ionization
,,,,,,,,,,,,,,,,,,,,,,,,,

Model appropriate for high fields or low frequencies.

.. py:function:: new adk ionization [<optional keys>] [for] <name> { directives }

	The following directives are supported:

		Shared directives: see above

PPT Photoionization
,,,,,,,,,,,,,,,,,,,,

Cycle-averaged model that works across multi-photon and tunneling regimes.  Cannot be used for ionization due to carrier-resolved fields, i.e., must be used with an enveloped field solver.

.. py:function:: new ppt ionization [<optional keys>] [for] <name> { directives }

	The following directives are supported:

		Shared directives: see above

		.. py:function:: terms = n

		 	:param int n: number of terms to keep in the PPT expansion.
