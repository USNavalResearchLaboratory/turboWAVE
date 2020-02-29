Input File: Base
================

General Rules
--------------

#.	Comments may be in C or C++ style
#.	Case matters
#.	Keywords are ``new``, ``generate``, ``get``, ``for``, and ``open``.
#.	Names assigned to objects must not contain white space, or any of the characters :samp:`/={}()`
#.  Any word starting with the :samp:`%` character triggers a macro substitution (see :ref:`preprocessor`)
#.	Any uninterrupted sequence of white space is equivalent to any other (same as in C)
#.	Parenthesis or commas may be freely exchanged with any white space character

	* :samp:`( 1 , 2 , 3 )` is the same as :samp:`1 2 3`
	* Do not confuse parenthesis with curly braces, cf. :samp:`()` and :samp:`{}`

#.	Spaces around an equals sign or curly braces may be omitted

 	* :samp:`x=1` is the same as :samp:`x = 1`

#.	All boolean variables can be set to an affirmative value using either :samp:`true`, :samp:`yes`, or :samp:`on`, or a negative value using :samp:`false`, :samp:`no`, or :samp:`off`.
#.	Variable length lists are enclosed by curly braces

	* ``{ 0 , 1 , 2 , 3 }`` is a list that can have any number of elements

.. _preprocessor:

Input Preprocessor
------------------

The turboWAVE input file is run through a preprocessor.  The features provided are as follows.

File Substitution
,,,,,,,,,,,,,,,,,

While preprocessing the input file, the contents of another file can be inserted.  The format is similar to that used by the C preprocessor.  For example,

.. code-block:: c

	#include myfile.tw

would substitute the contents of ``myfile.tw`` at the point in the file where the ``#include`` directive appears.  This can be done recursively. The ``#include`` directive may appear anywhere in the input file, except where it would interrupt another directive.

User Defined Macros
,,,,,,,,,,,,,,,,,,,

The effect of user variables can be achieved via macro substitution.  The format is the same as that used by the C preprocessor.  For example,

.. code-block:: c

	#define $r0 2.5

causes every subsequent occurrence of ``$r0`` to be replaced with ``2.5``.  The use of the ``$`` prefix is optional, but highly recommended, as it helps prevent unintended substitutions, and improves readability (including syntax highlights in supported editors).

The analogy with the C preprocessor is limited.  Function-like macros are not supported.  The substitution value cannot contain any white space characters.  The substitution is unconditional, e.g., if the key occurs as a word in a string it is replaced.

User macros can be defined at any point in an input file, except where they would interrupt another directive. Attempting to redefine a macro throws an error.

Unit Conversion
,,,,,,,,,,,,,,,

Almost all input parameters are in normalized units.  However, there are several pre-defined macros that make it simple to use physical units.  These are triggered by the ``%`` character. The format is :samp:`%{n}{u}`, where :samp:`{n}` is a number and :samp:`{u}` is a string identifying the units.  An example is :samp:`%10ps`, which means 10 picoseconds. No spaces may appear in the macro.  Supported units and identifier string are:

	* micrometers = um
	* millimeters = mm
	* centimeters = cm
	* meters = m
	* femtoseconds = fs
	* picoseconds = ps
	* nanoseconds = ns
	* microseconds = us
	* seconds = s
	* degrees = deg
	* radians = rad
	* milliradians = mrad
	* microradians = urad
	* Particles per cubic meter = m-3
	* Particles per cubic centimeter = cm-3
	* Joules per cubic centimeter = Jcm3
	* Joules per cubit meter = Jm3
	* electron volts = eV
	* Kelvin = K
	* CGS cross section = cm2
	* MKS cross section = m2
	* CGS diffusivity = cm2s
	* MKS diffusivity = m2s

Preprocessor Order
,,,,,,,,,,,,,,,,,,

The order of preprocessor operations is as follows:

	#. Strip comments
	#. Recursive file substitution

		* Comments are stripped at each level

	#. Clean white space
	#. Process user defined macros

		* At present keys must be unique across all included files.

	#. Process predefined macros

Top Level Directives
--------------------

Top level directives tend to come first in an input file.  They are not contained in any other block or structure.

.. py:function:: hardware acceleration device string = dev

	Use hardware accelerators having the given substring in their name

	:param str dev: the substring to search for in the device name, e.g., ``radeon``.  Case doesn't matter.

.. py:function:: hardware acceleration device numbers = dev_list

	Optional specification of preferred OpenCL device numbers.  If specified these take precedence over name search.

	:param list dev_list: variable length list of integers, e.g., ``{ 0 , 1 , 2 }``

.. py:function:: hardware acceleration platform string = platform

	Use only OpenCL platforms having the given substring in their name

	:param str platform: the substring to search for in the platform name, e.g., ``cuda``.  Case doesn't matter.

.. py:function:: open restart file name

	opens the named restart file.  If this command is used, all others are optional. Opening the restart file automatically creates the grid and any waves, pulses, species, or particles that were present when the restart file was saved. Subsequent directives override or add to the data in the restart file.

	:param str name: the name of the restart file to load

.. py:function:: unit density = dens

	Sets the unit density.  Fixes the normalization if needed.

	:param float dens: the density in particles per cubic centimeter

.. py:function:: steps = s

	:param int s: the number of simulation cycles to execute before terminating

.. py:function:: timestep = dt

	:param float dt: the timestep in units of :math:`\omega_p^{-1}`

.. py:function:: dtmin = dtm

	:param float dtm: if adaptive timestep in use, don't let it become less than this

.. py:function:: dtmax = dtx

	:param float dtx: if adaptive timestep in use, don't let it become greater than this

.. py:function:: dtcrit = dtc

	:param float dtc: if adaptive timestep falls below this value, switch to a fixed timestep.  The fixed timestep is taken from the ``timestep`` directive.

.. py:function:: maxtime = tm

	:param float tm: stop simulation after this much simulated time (useful with adaptive timestep)

.. py:function:: neutralize = n

	:param bool n: if yes, this causes an equal and opposite fixed charge to be added to the grid for every particle created.

.. py:function:: window speed = v

	:param float v: If moving window = yes, speed that lab frame quantities move back.  If moving window = no, speed that light frame quantities move forward.

.. py:function:: moving window = mv

	:param bool mv: Whether or not to move the lab frame quantities backward at the window speed. If no, light frame quantities are moved forward at the window speed.

.. py:function:: dump period = dp

	:param int dp: steps before dumping restart file

.. py:function:: append mode = am

	if simulation is restarted, append data to diagnostic files rather than overwrite them

	:param bool am: Use append mode if true.  Default is false.

.. py:function:: output level = lvl

	:param int lvl: If 0 then only MPI rank 0 writes an output file (to stdout).  If lvl > 0 than every MPI process produces an output file.

.. _boundaries:
.. py:function:: xboundary = ( b1 , b2 )

	Boundary conditions for whole simulation at the extremities in the x-coordinate. Can be overridden by individual modules. Parameters take values ``absorbing``, ``periodic``, ``emitting``, ``reflecting``, ``axisymmetric``, ``ejecting``.

	:param enum b1: Boundary condition of the low side.
	:param enum b2: Boundary condition on the high side.

.. py:function:: yboundary = ( b1 , b2 )

	Boundary conditions for whole simulation at the extremities in the y-coordinate, see xboundary.

.. py:function:: zboundary = ( b1 , b2 )

	Boundary conditions for whole simulation at the extremities in the z-coordinate, see xboundary.

Object Creation
---------------

There is a general form for creating objects:

.. _block-create:
.. py:function:: new <key1> [<key2> <key3> ...] [<name>] { <directives> }
.. py:function:: new <key1> [<key2> <key3> ...] for <name> { <directives> }
.. py:function:: generate <key1> [<key2> <key3> ...] [for] <name> { <directives> }

This entire construct is called a block.  The start of the block is signaled by a keyword, either ``new`` or ``generate``.  The next several words are ordered keys.  The keys are used to identify the type of object requested.  The user is free to add any number of trailing keys.  In the first form, the last word is the user-defined name of the new object.  The second form allows the new object to be associated with a parent object.  In this case the given name refers to the parent object, while the name of the new object is assigned automatically based on the keys. The third form is merely an alternative syntax for the second form, in which the ``for`` keyword is not required. Finally, there is a set of directives enclosed by curly braces.

.. Tip::

	Enclosing the name in quotes is not required, but can be helpful in distinguishing your own names from internal names.  If you use quotes remember they must be used everywhere the name is referenced.  Using quotes does not change the naming rules, e.g., you still may not use spaces.

Some low level objects use ordered directives.  Then the form is typically

.. py:function:: new <key1> [<key2> <key3> ...] = <directives>

In this case the directives are all required and must be in the right order.

Associating Objects
-------------------

Objects may be related by a containment hierarchy.  There are three ways to express this.

Nested Declarations
,,,,,,,,,,,,,,,,,,,

To use nested declarations, simply create the new object using the ``new`` command from within the directives block of the higher level object:

.. code-block:: none

	new direct field solver 'em'
	{
		new hermite gauss 'HG00'
		{
			// fill in directives defining the mode
		}
	}

Pre-declaration
,,,,,,,,,,,,,,,

To use a predeclaration, create a named low level object.  Then add it to a higher level object with a directive:

.. code-block:: none

	new hermite gauss 'HG00'
	{
		// fill in directives defining the mode
	}
	new direct field solver 'em'
	{
		get 'HG00'
	}

Post-declaration
,,,,,,,,,,,,,,,,

To use a post-declaration use one of the two associative forms of object creation:

.. code-block:: none

	new species 'ions'
	{
		// fill in directives defining the species
	}
	generate uniform 'ions'
	{
		// fill in directives defining the profile
	}

Numerical Grid
--------------

.. py:function:: new grid { directives }

	There must be exactly one grid block, which defines the numerical grid for all modules.

	:param block directives: The following directives are supported:

		.. py:function:: geometry = g

			:param enum g: can be ``cartesian``, ``cylindrical``, ``spherical``

		.. py:function:: corner[ijk] = ( x0 , y0 , z0 )

			Coordinates of the given vertex of the grid region.  If the optional ``ijk`` are omitted the vertex is the one where all coordinates are minimum.  Otherwise ``ijk`` is a binary code identifying one of eight vertices. Only one vertex may be given, otherwise the geometry is over-specified.  The coordinates are not necessarily Cartesian, but rather in the coordinate system of the grid.

			:param binary ijk: three binary digits, 0 indicates low side, 1 indicates high side.  For example, 011 means low x-side, high y-side, and high z-side.  Can be omitted, defaults to 000.
			:param float x0: The first coordinate of the corner
			:param float y0: the second coordinate of the corner
			:param float z0: the third coordinate of the corner

		.. py:function:: dimensions = (Nx,Ny,Nz)

			Dimensions of the grid region in numbers of cells along the three coordinate axes.

			:param int Nx: cells along the first coordinate
			:param int Ny: cells along the second coordinate
			:param int Nz: cells along the third coordinate

		.. py:function:: cell size = (dx,dy,dz)

			The cell size is given in parameter space, i.e., it could be an arc length or an angular sweep.

			:param float dx: length of cell edge along first coordinate
			:param float dy: length of cell edge along second coordinate
			:param float dz: length of cell edge along third coordinate


		.. py:function:: decomposition = ( Dx , Dy , Dz )

			Number of cuts of the domain along each coordinate.  This determines how the domain is split across parallel tasks.  The number of MPI tasks should be set to the product of all three parameters.

			:param int Dx: cuts along the first coordinate
			:param int Dy: cuts along the second coordinate
			:param int Dz: cuts along the third coordinate

		.. py:function:: radial progression factor = rpf

			:param float rpf: radial cells start to increase by this factor after the first 1/3 of radial cells

		.. py:function:: region : start = s , end = e , length = l

			Create a non-uniform grid region along the z-coordinate.
			The z-width of the cells grows according to a quintic polynomial.

			:param int s: the cell where the non-uniform grid region starts
			:param int e: the cell where the non-uniform grid region ends
			:param float l: the total length of the non-uniform grid region.  The cells sizes are adjusted to give the requested length.

		.. py:function:: adaptive timestep = at

			:param bool at: whether or not to use an adaptive time stepping scheme.


Radiation Injection
-------------------

To inject radiation, you specify a type of electromagnetic mode, directives defining its particular parameters, and attach it to a field solver.  For an in depth description of the available radiation modes see :doc:`bak-em-modes`. If the wave starts inside the box, an elliptical solver may be used to refine the initial divergence. If the wave starts outside the box, it will be coupled in, provided the field solver supports this. Each wave object has its own basis vectors :math:`({\bf u},{\bf v},{\bf w})`, with :math:`{\bf u}` the electric field polarization direction and :math:`{\bf w}` the propagation direction. All the available modes respond to the same set of directives. These are as follows:

.. _wave-obj:
.. py:function:: new <key> for <field_solver_name> { directives }

	:param str key: The key determines the type of mode.  Valid keys are ``plane wave``, ``hermite gauss``, ``laguerre gauss``, ``bessel beam``, ``airy disc``, and ``multipole``.
	:param str field_solver_name: Name given to a previously defined field solver module.
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

.. note::

	In the past there was a distinction between carrier resolved and enveloped radiation injection objects.  This distinction has been retired.  Envelope treatment is triggered automatically by attaching any radiation injection object to a enveloped field solver.

.. _eos:

Equation of State Shared Directives
-----------------------------------

:doc:`Equation of State <bak-eos>` (EOS) models are needed for hydrodynamics simulation.  EOS tools may be created inside a module block on the fly, using the following shared directives.  They may also be created at the root level as named tools (not covered here).

.. note::
	As of this writing EOS is a moving target.  The interface may change.

.. py:function:: eos = ideal-gas

	Directs a module to use the ideal gas equation of state

.. py:function:: eos = hot-electrons

	Directs a module to use the ideal gas equation of state along with Braginskii electron transport coefficients

.. py:function:: eos = simple-mie-gruneisen , gruneisen parameter = grun

	Directs a module to use the simplified mie-gruneisen equation of state

	:param float grun: the gruneisen parameter relating density, temperature, and pressure

.. py:function:: eos = linear-mie-gruneisen , subdirectives

	Directs a module to use the linear Hugoniot-based mie-gruneisen equation of state. The sub-directives are processed by the enclosing module and can be treated as any other module directive, so long as they come after ``linear-mie-gruneisen``.

	:param block subdirectives: the following subdirectives are supported:

		.. py:function:: gruneisen parameter = grun

			:param float grun: the gruneisen parameter relating density, temperature, and pressure

		.. py:function:: reference density = nref

			:param float nref: the reference density for the Hugoniot data

		.. py:function:: hugoniot intercept = c0

			:param float c0: y-intercept of the Hugoniot curve, typically the speed of sound

		.. py:function:: hugoniot slope = s1

			:param float s1: slope of the Hugoniot curve at the reference density

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

.. py:function:: generate uniform <name> { directives }

	Generate uniform density within the clipping region.

	:param str name: name of module defining type of matter to load.
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`matter-loading-shared`

		.. py:function:: density = n0

			:param float n0: density to load


.. py:function:: generate piecewise <name> { directives }

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

.. py:function:: generate channel <name> { directives }

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


.. py:function:: generate column <name> { directives }

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

.. py:function:: generate gaussian <name> { directives }

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

.. py:function:: new conductor <name> { directives }

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

Diagnostics
------------

Diagnostic Formats
,,,,,,,,,,,,,,,,,,

TurboWAVE uses simple text and binary formats.  Text files are generally tab delimited tables of ASCII data, with a one-line header containing column labels.  There are two binary formats.

.. note::

	A correctly compiled TurboWAVE executable always writes binary data in big-endian format.  All data readers should assume every turboWAVE binary is big-endian.

.. highlight:: none

DataViewer Box Diagnostic Format::

	The string "DataViewer 2.0.0"
	32 bit integer : x dimension
	32 bit integer : y dimension
	32 bit integer : z dimension
	32 bit float : coordinate of lower bound in x
	32 bit float : coordinate of upper bound in x
	32 bit float : coordinate of lower bound in y
	32 bit float : coordinate of upper bound in y
	32 bit float : coordinate of lower bound in z
	32 bit float : coordinate of upper bound in z
	3D array of 32 bit floats: frame 1, written in FORTRAN order
	3D array of 32 bit floats: frame 2, written in FORTRAN order
	...more frames (size of files allows reader to determine frames)

DataViewer Orbit Diagnostic Format::

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

	:param str f: name of the file to write. Actual file names may be prepended with the name of some subset of the overall data associated with the diagnostic (some diagnostics write multiple files).  This may be postpended with a filename extension such as ``.txt``, ``.dvdat`` or ``.dvpar``.  The special name ``full`` causes the files to have only the prepended string and the extension in their names.  This is the default.

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

Specific Diagnostics
,,,,,,,,,,,,,,,,,,,,

.. py:function:: new box diagnostic { directives }

	Write out grid data as sequence of frames.  Clipping region must be a simple box.
	This diagnostic produces several files per module.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`diagnostics-shared`

		.. py:function:: average = tst

			:param bool tst: average over sub-grid, or not.  If not, diagnose lower corner cell only.

		.. py:function:: skip = ( sx , sy , sz )

			Defines a reduced grid produced by downsampling the full grid.  The reduction factor is the product of the three skipping parameters.  Note the centroid of the sampling points is shifted.

			:param int sx: advance this many cells in the x-direction between writes
			:param int sy: advance this many cells in the y-direction between writes
			:param int sz: advance this many cells in the z-direction between writes


.. py:function:: new energy series { directives }

	Diagnostic of volume integrated quantities.  Normalization includes the unit of particle number.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`diagnostics-shared`

		.. py:function:: precision = digits

		 	:param int digits: number of digits used to represent each result


.. py:function:: new point series { directives }

	Diagnostic to write out grid data at a specific point.

	:param block directives: The following directives are supported:

		Shared directives: see :ref:`diagnostics-shared`

		.. py:function:: point = (Px,Py,Pz)

			Coordinates of the point to diagnose.

		.. py:function:: move with window = tst

			:param bool tst: if true the point moves with the window


.. py:function:: new phase space plot for <species_name> { directives }

	Diagnostic to write out 2D phase space projections.
	Phase space variables include ``x``, ``y``, ``z``, ``px``, ``py``, ``pz``, ``mass``, ``energy``
	Here, ``mass`` is the total relativistic energy, ``energy`` is the kinetic part.

	:param str species_name: the name of the species to diagnose
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`diagnostics-shared`

		.. py:function:: abcissa = var

			:param enum var: the phase space variable to associate with the x axis

		.. py:function:: ordinate = var

			:param enum var: the phase space variable to associate with the y axis

		.. py:function:: minimum = (xmin,ymin)

			:param float xmin: the lower bound of the x axis
			:param float ymin: the lower bound of the y axis

		.. py:function:: maximum = (xmax,ymax)

			:param float xmax: the upper bound of the x axis
			:param float ymax: the upper bound of the y axis

		.. py:function:: dimensions = (Nx,Ny)

			:param int Nx: the number of cells in the x direction for the phase space grid
			:param int Ny: the number of cells in the y direction for the phase space grid


.. py:function:: new orbit diagnostic for <species_name>

	Diagnostic to write out full phase space data of the particles.

	.. caution::
		Orbit diagnostics can create excessively large files if not used carefully.  To avoid this, define a species with a small number of test particles and use this on them.

	:param str species_name: the name of the species to diagnose
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`diagnostics-shared`

		.. py:function:: minimum gamma = gmin

			:param float gmin: only save data for particles with gamma greater than this
