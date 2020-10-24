Input File: Tools
=================

General Information
-------------------

TurboWAVE objects come in two flavors, Modules and Tools (not to be confused with post-processing tools).  Modules are the high level objects that orchestrate a simulation.  Tools (internally called ``ComputeTool`` objects) are lower level objects dedicated to specific computations. Tools are intended to be attached to one or more modules (modules can share the same tool).  Modules can be placed in a hierarchy.  Tools can only be attached to a module.

Tools are created like any object, see :doc:`ref-input`.

Many modules are able to create their tools automatically.  This page emphasizes those that tend to benefit from explicit management by the user.

.. _radiation:

Radiation Injection
-------------------

To inject radiation, you specify a type of electromagnetic mode, directives defining its particular parameters, and attach it to a field solver.  For an in depth description of the available radiation modes see :doc:`bak-em-modes`. If the wave starts inside the box, an elliptical solver may be used to refine the initial divergence. If the wave starts outside the box, it will be coupled in, provided the field solver supports this. Each wave object has its own basis vectors :math:`({\bf u},{\bf v},{\bf w})`, with :math:`{\bf u}` the electric field polarization direction and :math:`{\bf w}` the propagation direction. All the available modes respond to the same set of directives. These are as follows:

.. _wave-obj:
.. py:function:: new <key> [<name>] [for <module_name>] { <directives> }

	:param str key: The key determines the type of mode.  Valid keys are ``plane wave``, ``hermite gauss``, ``laguerre gauss``, ``bessel beam``, ``airy disc``, and ``multipole``.
	:param str module_name: Name given to a previously defined field solver module.
	:param block directives: The following directives are supported:

		.. py:function:: direction = ( nx , ny, nz )

			:param float nx: first component of :math:`{\bf w}` in standard basis.
			:param float ny: second component of :math:`{\bf w}` in standard basis.
			:param float nz: third component of :math:`{\bf w}` in standard basis.

		.. py:function:: a = ( ax , ay , az )

			If the peak vector potential is :math:`a_0`, then :math:`{\bf a} = a_0{\bf u}`.
			TurboWAVE will force transversality by making the replacement :math:`{\bf a} \rightarrow {\bf w}\times{\bf a}\times{\bf w}`

			:param float ax: first component of :math:`{\bf a}` in standard basis
			:param float ay: second component of :math:`{\bf a}` in standard basis
			:param float az: third component of :math:`{\bf a}` in standard basis

		.. py:function:: focus position = ( fx , fy , fz )

			:param float fx: first focal position coordinate in standard basis
			:param float fy: second focal position coordinate in standard basis
			:param float fz: third focal position coordinate in standard basis

		.. py:function:: w = w0

			:param float w0: central frequency of the wave

		.. py:function:: refractiveindex = n0

			:param float n0: refractive index in the starting medium

		.. py:function:: chirp = c0

			:param float c0: creates a chirp :math:`\exp (-ic_0 t^2)`, with time referenced so that the center frequency occurs at the end of the risetime.  Up-chirp results from :math:`c_0>0`.

		.. py:function:: phase = p0

			:param float p0: phase shift in radians

		.. py:function:: delay = t0

			:param float t0: Front of wave reaches focus position after this amount of time

		.. py:function:: risetime = t1

		.. py:function:: holdtime = t2

		.. py:function:: falltime = t3

		.. py:function:: r0 = ( u0 , v0 )

			:param float u0: spot size in the :math:`{\bf u}` direction.  Note this is **not necessarily** the spot size in the first coordinate of the standard basis. Spot size is measured at :math:`1/e` point of the field amplitude.
			:param float v0: spot size in the :math:`{\bf v}` direction.

		.. py:function:: mode = ( mu , mv )

			Transverse mode numbers, different meanings depending on the mode type.

			:param int mu: mode number in the :math:`{\bf u}` direction
			:param int mv: mode number in the :math:`{\bf v}` direction

		.. py:function:: exponent = ( m , n )

			This directive applies only to the paraxial beam modes, Hermite and Laguerre.

			:param int m: exponent to use in transverse profile, default is 2 (standard Gaussian). If even induces order *m* supergaussian, if odd induces order *m+1* cosine.
			:param int n: If the mode is Hermite then *n* applies to the v-direction.  If it is Laguerre then *n* is ignored.

		.. py:function:: shape = pulse_shape

			:param enum pulse_shape: determines the shape of the pulse envelope, can be ``quintic`` (default), ``sin2``, ``sech``

		.. py:function:: boosted frame gamma = g

			:param float g: relativistic Lorentz factor of the boosted frame (default=1).  If g>1, turboWAVE will transform the wave into the boosted frame.  The parameters describing the wave should all be given in lab frame coordinates.  The grid coordinates are taken as the boosted frame.  At present this feature should only be used for paraxial modes propagating along the z-axis.

		.. py:function:: zones = z

			:param int z: create a superposition of transversely periodic modes across ``z`` zones in each dimension.  The number of zones should be odd, and large enough to span several spot sizes.  This is useful for boosted frame simulations where the Rayleigh length is much shorter than the pulse duration.

.. note::

	In the past there was a distinction between carrier resolved and enveloped radiation injection objects.  This distinction has been retired.  Envelope treatment is triggered automatically by attaching any radiation injection object to a enveloped field solver.

.. _matter-loading:

Matter Loading
--------------

The loading of matter into the simulation box is done using ``generate`` blocks.  These take the same form whether we are loading particles or fluid elements.  In loading matter it is important to distinguish the clipping region from the profile:

.. glossary::

	clipping region
		A clipping region is a filter that multiplies a physical quantity by zero outside the region, and unity inside.

	profile
		A profile is a spatial distribution of some intrinsic parameter such as density.

.. note::
	Our definition of thermal velocity is :math:`f(v) = f_0\exp(-v^2/2v_{th}^2)`

.. note::
	For isotropic distributions we have :math:`kT = mv_{th}^2`, :math:`v_i^{rms} = v_{th}`, and :math:`v_{tot}^{rms} = \sqrt{3}v_{th}`.

.. _matter-loading-shared:

Matter Loading Shared Directives
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

The following directives may be used with any profile type

.. py:function:: clipping region = name

 	Load the matter only within the specified geometric region.  See :doc:`ref-geometry` for documentation on creating complex geometric regions.

	:param str name: the name of the geometric region to use

.. py:function:: position = ( x , y , z )

 	Specify where to put profile’s reference point, typically extremum of profile.  For piecewise profiles this is interpreted as a translation.

	.. tip::
		This does not affect the position of the clipping region, only the profile.

.. py:function:: euler angles = ( qx , qy , qz )

	Rotation of the profile about the profile position.

	.. tip::
		This does not affect the rotation of the clipping region, only the profile.

.. py:function:: temperature = T

 	:param float T: initial temperature of the matter

.. py:function:: thermal momentum = (pthx,pthy,pthz)

.. py:function:: drift momentum = (px,py,pz)

.. py:function:: loading = lmethod

 	:param enum lmethod: loading method.  takes values ``deterministic``, ``statistical``

.. py:function:: particle weight = wscheme

 	:param enum wscheme: takes values ``variable``, ``fixed``

.. py:function:: type = profile_type

	Matter loading encompasses mass, energy, and momentum.  The type of profile determines which quantity is loaded.

 	:param enum profile_type: takes values ``density``, ``energy``, ``px``, ``py``, ``pz``

.. py:function:: timing = timing_type

	:param enum timing_type: takes values ``triggered`` or ``maintained`` (default = triggered). Triggered profiles are additive.  Maintained profiles try to hold fixed conditions.

.. py:function:: t0 = start_time

	:param float start_time: time at which matter loading begins.

.. py:function:: t1 = stop_time

	:param float stop_time: time at which matter loading ends.  If timing is ``triggered`` this is ignored.

.. py:function:: boosted frame gamma = g

	:param float g: relativistic Lorentz factor of the boosted frame (default=1).  If g>1, turboWAVE will transform the profile and clipping region into the boosted frame.  The parameters describing the profile and region should all be given in lab frame coordinates.  The grid coordinates are taken as the boosted frame.  When using a boosted frame the ``neutralize`` top level directive must be ``false``.

Specific Matter Loading Profiles
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

.. py:function:: generate uniform <name> { <directives> }

	Generate uniform density within the clipping region.

	:param str name: name of module defining type of matter to load.
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`matter-loading-shared`

		.. py:function:: density = n0

			:param float n0: density to load


.. py:function:: generate piecewise <name> { <directives> }

	Generate piecewise varying density within the clipping region.  The total density is the product of 3 piecewise functions:

		:math:`n(x,y,z) = X(x)Y(y)Z(z)`

	:param str name: name of module defining type of matter to load.
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`matter-loading-shared`

		.. py:function:: xpoints = x_list

			:param list x_list: Variable length list of floating point numbers giving the points at which :math:`X(x)` is known, e.g., ``{ 0 , 1.5 , 3.4 , 5.1 }``.

		.. py:function:: ypoints = y_list

			:param list y_list: Variable length list of floating point numbers giving the points at which :math:`Y(y)` is known, e.g., ``{ 0 , 1.5 , 3.4 , 5.1 }``.

		.. py:function:: zpoints = z_list

			:param list z_list: Variable length list of floating point numbers giving the points at which :math:`X(x)` is known, e.g., ``{ 0 , 1.5 , 3.4 , 5.1 }``.

		.. py:function:: xdensity = xd_list

			:param list xd_list: Variable length list of floating point numbers giving the values of :math:`X(x)` at the points listed with ``xpoints``.

		.. py:function:: ydensity = yd_list

			:param list yd_list: Variable length list of floating point numbers giving the values of :math:`Y(y)` at the points listed with ``ypoints``.

		.. py:function:: zdensity = zd_list

			:param list zd_list: Variable length list of floating point numbers giving the values of :math:`Z(z)` at the points listed with ``zpoints``.

		.. py:function:: shape = my_shape

			:param enum my_shape: ``quintic``, ``quartic``, ``triangle``

		.. py:function:: symmetry = sym

		 	:param enum sym: ``none``, ``cylindrical``, ``spherical``.  If cylindrical, x-profile is interpreted as radial, z-profile is axial, y is only used to define origin. If spherical, x-profile is radial, y and z are used only to define the origin.

		.. py:function:: mode number = nx ny nz

		 	Multiply final profile by :math:`\left[\cos(n_x x/2)\cos(n_y y/2)\cos(n_z z/2)\right]^2`

.. py:function:: generate channel <name> { <directives> }

	Generate density channel within the clipping region.  The defining formula is

		:math:`n(x,y,z) = Z(z)\left(n_0 + n_2\rho^2 + n_4\rho^4 + n_6\rho^6\right)`

		:math:`\rho = \sqrt{x^2 + y^2}`

		The matched beam condition for spot size :math:`\rho_0` is

		:math:`n_2 = 1/\pi r_e \rho_0^4`

		where :math:`r_e` is the classical electron radius, :math:`n_0` is arbitrary, and higher terms vanish.  The normalization is

		:math:`n_i \rightarrow \frac{n_i}{n} \left(\frac{c}{\omega}\right)^i`

		where :math:`\omega` is the unit frequency and :math:`n` is the unit density.  This leads to the matched beam condition in normalized units as

		:math:`n_2 = 4/\rho_0^4`

	:param str name: name of module defining type of matter to load.
	:param block directives: The following directives are supported:

		Shared directives:
			see :ref:`matter-loading-shared`

			piecewise profile :math:`Z(z)` function

			piecewise profile ``shape`` directive.

		.. py:function:: coefficients = n0 n2 n4 n6

			:param float n0: see :math:`n_0` in defining formula
			:param float n2: see :math:`n_2` in defining formula
			:param float n4: see :math:`n_4` in defining formula
			:param float n6: see :math:`n_6` in defining formula


.. py:function:: generate column <name> { <directives> }

	Generate density column within the clipping region.

		:math:`n(x,y,z) = Z(z)\exp(-x^2/\sigma_x^2 - y^2/\sigma_y^2)`

	:param str name: name of module defining type of matter to load.
	:param block directives: The following directives are supported:

		Shared directives:
			see :ref:`matter-loading-shared`

			piecewise profile :math:`Z(z)` function

			piecewise profile ``shape`` directive.

		.. py:function:: size = ( sx , sy , sz )

			:param float sx: radius of column, per :math:`\sigma_x` in defining formula.
			:param float sy: radius of column, per :math:`\sigma_y` in defining formula.
			:param float sz: ignored.

.. py:function:: generate gaussian <name> { <directives> }

	Generate a Gaussian ellipsoid within the clipping region.

		:math:`n(x,y,z) = n_0 \exp(-x^2/\sigma_x^2 - y^2/\sigma_y^2 - z^2/\sigma_z^2)`

	:param str name: name of module defining type of matter to load.
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`matter-loading-shared`

		.. py:function:: density = n0

			:param float n0: peak density, per defining formula.

		.. py:function:: size = ( sx , sy , sz )

			:param float sx: :math:`\sigma_x` in defining formula.
			:param float sy: :math:`\sigma_y` in defining formula.
			:param float sz: :math:`\sigma_x` in defining formula.


.. _conductor:

Conducting Regions
------------------

Conducting regions serve the following purposes:

	1. Perfect conductors filling arbitrary cells in electromagnetic simulations
	2. Antenna objects in electromagnetic simulations
	3. Impermeable objects filling arbitrary cells in hydrodynamic simulations
	4. Fixed potential objects filling arbitrary cells in electrostatic simulations
	5. Fixed temperature objects filling arbitrary cells in hydrodynamic simulations

.. py:function:: new conductor [<name>] [for <module_name>] { <directives> }

	The electrostatic potential can be fixed within the conductor as

		:math:`\Phi(t) = \Phi_0 S(t) \cos(\omega t + \varphi)`

	The dipole radiator elements oscillate according to

		:math:`{\bf P}(t,x,y,z) = {\bf P}_0 S[T(t,x,y)] \sin[\omega T(t,x,y) + \varphi + {\bf k}_s \cdot {\bf r}]`

		:math:`T(t,x,y) = t + \frac{x^2+y^2}{2f}`

	:param str name: Name given to the conductor
	:param block directives: The following directives are supported:

		Shared directives:
			Temporal envelope :math:`S(t)` is derived from pulse shape parameters per :ref:`wave object <wave-obj>`

		.. py:function:: clipping region = name

			Rotation of clipping region also rotates current distribution

			:param str name: name of geometric region to use

		.. py:function:: temperature = T

			:param float T: The temperature of the conductor can be fixed at T, serving as a dirichlet boundary condition for heat equations. Note that parabolic solver tools default to a homogeneous neumann condition (no heat flow).

		.. py:function:: enable electrostatic = tst

			:param bool tst: this conductor will fix the potential

		.. py:function:: enable electromagnetic = tst

			:param bool tst: this conductor will reflect EM waves

		ANTENNA DIRECTIVES:
		Currents are driven with dipole oscillators.  This avoids problems with static field generation.  All the lists must be of equal length.  Each list element is an oscillator. The total current is the superposition of the current of each oscillator.

		.. py:function:: current type = curr_typ

		 	:param enum curr_typ: takes values ``electric``, ``magnetic``, or ``none``

		.. py:function:: potential = lst

			Determines :math:`\Phi_0` for each oscillator.

			:param list lst: variable length list of scalar potentials, e.g., ``{ 1.0 , 2.0 }``

		.. py:function:: px = lst1 , py = lst2 , pz = lst3

			Determines :math:`{\bf P}_0` for each oscillator.

		.. py:function:: w = w0

			Determines :math:`\omega` for each oscillator.

		.. py:function:: phase = p0

			Determines :math:`\varphi` for each oscillator.

		.. py:function:: f = f0

			:param float f0: Determines :math:`f` parameter that appears in :math:`T(t,x,y)`.  This is supposed to produce a focus at the corresponding distance from the antenna (default = infinity).

		.. py:function:: ks = ksx ksy ksz

		 	Apply linear phase variation to create tilted wave (default = 0).

		.. py:function:: gaussian size = ( sx , sy , sz )

			Apply a gaussian spatial weight to the oscillator amplitudes.

.. _elliptic:

Elliptic Solvers
----------------

All elliptic solvers share the following directives:

.. py:function:: <coord>boundary = ( <bc1> , <bc2> )

	:param enum coord: can be ``x``, ``y``, ``z``
	:param enum bc1: boundary condition on lower side, can be ``open``, ``dirichlet``, ``neumann``.
	:param enum bc2: boundary condition on lower side, can be ``open``, ``dirichlet``, ``neumann``.

Inhomogeneous boundary conditions are implemented by using conductor tools to fix the potential.  For external boundaries, the conductor must occupy the far ghost cells.

Iterative Solver
,,,,,,,,,,,,,,,,

.. py:function:: new iterative elliptic [<name>] [for <module_name>] { <directives> }

	Uses successive over-relaxation to iteratively solve the elliptic equation.  This solver is slow, but flexible.  There is no limit on the topology of the boundary conditions, and arbitrary coordinates are supported.  The following directives are supported:

		Shared directives: see base elliptic solver

		.. py:function:: tolerance = <tol>

			:param float tol: iterate until the residual is reduced to this level

		.. py:function:: overrelaxation = <ov>

			:param float ov: overrides the default overrelaxation parameter (not generally recommended)

FACR Solver
,,,,,,,,,,,,,,,,

.. py:function:: new facr elliptic [<name>] [for <module_name>] { <directives> }

	Uses Fourier analysis is in the transverse directions.  This solver is fast, but boundary conditions can only be imposed on constant z-surfaces, and Cartesian coordinates are required.  The following directives are supported:

		Shared directives: see base elliptic solver

Eigenmode Solver
,,,,,,,,,,,,,,,,

.. py:function:: new eigenmode elliptic [<name>] [for <module_name>] { <directives> }

	Uses generalized spectral resolution of the transverse coordinates.  This solver works in arbitrary coordinates, and is fast as long as the transverse modes are truncated.  Boundary conditions can only be imposed on constant z-surfaces.  The following directives are supported:

		Shared directives: see base elliptic solver

		.. py:function:: modes = <N>

			:param int N: The number of transverse modes to keep.  The modes are taken from an ordered list, sorted by magnitude of the eigenvalue.

.. _propagator:

Laser Propagator
----------------

Eigenmode Propagator
,,,,,,,,,,,,,,,,,,,,

.. py:function:: new eigenmode propagator [<name>] [for <module_name>] { <directives> }

	Uses generalized spectral resolution of the transverse coordinates.  This propagator works in arbitrary coordinates, and is fast as long as the transverse modes are truncated.  It has superior fidelity for highly dispersive systems.  The following directives are supported:

	.. py:function:: modes = <n>

		:param int n: maximum number of radial modes to keep (eigenmode propagator only)

	.. py:function:: damping time = <t>

		:param float t: e-folding time in the absorbing layers

	.. py:function:: absorbing layers = <l>

		:param int l: number of absorbing layers

ADI Propagator
,,,,,,,,,,,,,,,,,,,,

.. py:function:: new adi propagator [<name>] [for <module_name>] { <directives> }

	Uses alternating direction implicit method.  This is a fast propagator that works in arbitrary coordinates.  It has poor fidelity for highly dispersive systems.  There are no directives.

.. _ionization:

PhotoIonization
---------------

Ionization Shared Directives
,,,,,,,,,,,,,,,,,,,,,,,,,,,,

All the photoionization tools support the following directives:

.. py:function:: ionization potential = ip

	:param float ip: Ionization potential, units are specified as usual, e.g., ``ionization potential = 13.6 [eV]``

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

.. py:function:: new mpi ionization [<name>] [for <module_name>] { <directives> }

	The following directives are supported:

		Shared directives: see above

		.. py:function:: reference field = E0

		 	:param float E0: :math:`E_0`, where the MPI rate is proportional to :math:`(E/E_0)^{2l}`

ADK Tunneling Ionization
,,,,,,,,,,,,,,,,,,,,,,,,,

Model appropriate for high fields or low frequencies.

.. py:function:: new adk ionization [<name>] [for <module_name>] { <directives> }

	The following directives are supported:

		Shared directives: see above

PPT Photoionization
,,,,,,,,,,,,,,,,,,,,

Cycle-averaged model that works across multi-photon and tunneling regimes.  Cannot be used for ionization due to carrier-resolved fields, i.e., must be used with an enveloped field solver.

.. py:function:: new ppt ionization [<name>] [for <module_name>] { <directives> }

	The following directives are supported:

		Shared directives: see above

		.. py:function:: terms = n

			:param int n: number of terms to keep in the PPT expansion.

.. _eos:

Equation of State Tools
-----------------------

:doc:`Equation of State <bak-eos>` (EOS) models are needed for hydrodynamics simulation.  EOS models are encapsulated in tool objects that can be attached to appropriate modules in the usual way.

.. py:function:: new eos ideal gas tool [<name>] [for <module_name>] { <directives> }

	Directs a module to use the ideal gas equation of state.  No directives.

.. py:function:: new eos hot electrons [<name>] [for <module_name>] { <directives> }

	Directs a module to use the ideal gas equation of state along with Braginskii electron transport coefficients.  No directives.

.. py:function:: new eos simple mie gruneisen [<name>] [for <module_name>] { <directives> }

	Directs a module to use the simplified mie-gruneisen equation of state.  The following directives are supported:

		.. py:function:: gruneisen parameter = grun

			:param float grun: the gruneisen parameter relating density, temperature, and pressure

.. py:function:: new eos linear mie gruneisen [<name>] [for <module_name>] { <directives> }

	Directs a module to use the linear Hugoniot-based mie-gruneisen equation of state.

	The following directives are supported:

		.. py:function:: gruneisen parameter = grun

			:param float grun: the gruneisen parameter relating density, temperature, and pressure

		.. py:function:: reference density = nref

			:param float nref: the reference density for the Hugoniot data

		.. py:function:: hugoniot intercept = c0

			:param float c0: y-intercept of the Hugoniot curve, typically the speed of sound

		.. py:function:: hugoniot slope = s1

			:param float s1: slope of the Hugoniot curve at the reference density

Diagnostics
------------

.. Note::

	If a diagnostic tool is not explicitly attached to any module, it will be automatically attached to all modules.  This is an optimization based on the observation that one would often like a similar type of output from all modules.

Diagnostic Formats
,,,,,,,,,,,,,,,,,,

TurboWAVE binaries are in numerical Python (numpy) format (extension ``.npy``).  They can be easily read into a Python program using ``numpy.load``.

All metadata is in the file ``tw_metadata.json``, which can easily be read into a Python dictionary.  One then looks up the file of interest to expose further dictionaries pertaining to the file.

Text files are generally tab delimited tables of ASCII data, with a one-line header containing column labels.

.. highlight:: none

Box diagnostics produce a ``.npy`` file containing a four-dimensional array with axes (t,x,y,z).  Phase space diagnostics write a similar array, but the axes can have different meanings.  The grid data is found by looking up the ``'grid'`` key.  This returns a string with the name of the grid file containing the mesh point and time level information::

	t = <t0>
	axis1 = <x1> <x2> ... <xN>
	axis2 = <y1> <y2> ... <yN>
	axis3 = <z1> <z2> ... <zN>
	<...more frames>

For static grids, the additional frames only have the line giving the elapsed time.

The orbit diagnostic produces a ``.npy`` file containing the detailed particle information::

	Particle record 1 at time level 1
	Particle record 2 at time level 1
	...more particle records at time level 1
	Particle record N at time level 1
	Time level separator
	Particle record 1 at time level 2
	Particle record 2 at time level 2
	...more particle records at time level 2
	Particle record M at time level 2
	...more time level separators and particle records

Each particle record is an 8 element vector (x,px,y,py,z,pz,aux1,aux2).
The order of the particles within a time level is not significant.
Particles must be identified by unique values of aux1 and aux2.
The time level separator is a record with all zeros.
Valid particles can never have aux1 = aux2 = 0.

.. _diagnostics-shared:

Diagnostics Shared Directives
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

The following directives may be used with any diagnostic

.. py:function:: filename = f

	:param str f: name of the file to write. Actual file names may be prepended with the name of some subset of the overall data associated with the diagnostic (some diagnostics write multiple files).  This may be postpended with a filename extension such as ``.txt`` or ``.npy``.  The special name ``full`` causes the files to have only the prepended string and the extension in their names.  This is the default.

.. py:function:: clipping region = name

 	write data only within the specified geometric region.  See :doc:`ref-geometry` for documentation on creating complex geometric regions.  For some diagnostics there is a restriction on the complexity of the region.

	:param str name: the name of the geometric region to use

.. py:function:: t0 = start_time

	:param float start_time: time at which diagnostic write-out begins (default=0).

.. py:function:: t1 = stop_time

	:param float stop_time: time after which diagnostic write-out ends (default=infinity).

.. py:function:: period = steps

	:param int steps: number of simulation cycles between write-outs.

.. py:function:: time period = duration

	:param float duration: simulated time between write-outs, overrides ``period`` if specified.  If an adaptive time step is in use, this can approximate uniform spacing of write-outs.

.. py:function:: galilean velocity = (vx,vy,vz)

	Transform output to a Galilean frame, i.e., :math:`{\bf r}' = {\bf r} - {\bf v}t`.

	:param float vx: x-component of the galilean transformation velocity
	:param float vy: y-component of the galilean transformation velocity
	:param float vz: z-component of the galilean transformation velocity

.. _specific-diagnostics:

Specific Diagnostics
,,,,,,,,,,,,,,,,,,,,

.. py:function:: new box diagnostic [<name>] [for <module_name>] { <directives> }

	Write out grid data as sequence of frames.  Clipping region must be a simple box.
	This diagnostic produces several files per module, by default.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`diagnostics-shared`

		.. py:function:: average = tst

			:param bool tst: average over sub-grid, or not.  If not, diagnose lower corner cell only.

		.. py:function:: skip = ( sx , sy , sz )

			Defines a reduced grid produced by downsampling the full grid.  The reduction factor is the product of the three skipping parameters.  Note the centroid of the sampling points is shifted.

			:param int sx: advance this many cells in the x-direction between writes
			:param int sy: advance this many cells in the y-direction between writes
			:param int sz: advance this many cells in the z-direction between writes

		.. py:function:: reports = { <fields> }

			:param list fields: Put a list of fields to get restricted output.  If omitted then all available fields are written.


.. py:function:: new energy diagnostic [<name>] [for <module_name>] { <directives> }

	Diagnostic of volume integrated quantities.  Normalization includes the unit of particle number.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`diagnostics-shared`

		.. py:function:: precision = digits

		 	:param int digits: number of digits used to represent each result


.. py:function:: new point diagnostic [<name>] [for <module_name>] { <directives> }

	Diagnostic to write out grid data at a specific point.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`diagnostics-shared`

		.. py:function:: point = (Px,Py,Pz)

			Coordinates of the point to diagnose.  This is subject to the Galilean transformation (see :ref:`diagnostics-shared`).

.. py:function:: new phase space diagnostic [<name>] [for <module_name>] { <directives> }

	Diagnostic to write out up to 3D phase space projections.  Setting a dimension to 1 produces a lower dimensional projection.

	:param str species_name: the name of the species to diagnose
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`diagnostics-shared`

		.. py:function:: axes = (ax1,ax2,ax3)

			Determines the axes of the phase space.  Axes can be any of ``t``, ``x``, ``y``, ``z``, ``mass``, ``px``, ``py``, ``pz``, ``g``, ``gbx``, ``gby``, ``gbz``.  If a dimension collapses to 1, the axis is ignored, but still must be specified.

			:param enum ax1: the phase space variable to associate with axis 1
			:param enum ax2: the phase space variable to associate with axis 2
			:param enum ax3: the phase space variable to associate with axis 3

		.. py:function:: dimensions = (N1,N2,N3)

			:param int N1: cells along axis 1
			:param int N2: cells along axis 2
			:param int N3: cells along axis 3

		.. py:function:: bounds = (x0,x1,y0,y1,z0,z1)

			:param float x0: lower bound for axis 1
			:param float x1: upper bound for axis 1
			:param float y0: lower bound for axis 2
			:param float y1: upper bound for axis 2
			:param float z0: lower bound for axis 3
			:param float z1: upper bound for axis 3


.. py:function:: new orbit diagnostic [<name>] [for <module_name>]

	Diagnostic to write out full phase space data of the particles.

	.. caution::
		Orbit diagnostics can create excessively large files if not used carefully.  To avoid this, define a species with a small number of test particles and use this on them.

	:param str species_name: the name of the species to diagnose
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`diagnostics-shared`

		.. py:function:: minimum gamma = gmin

			:param float gmin: only save data for particles with gamma greater than this
