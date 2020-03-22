Input File: Base
================

.. _little_lang:

Little Language Quick Start
---------------------------

TurboWAVE input files are written using an intuitive "little language" which provides three basic functions:

	#. Create objects
	#. Associate objects
	#. Set parameters

There are only four keywords, **new**, **get**, **generate**, and **for**.  You can often get by with only **new**.  The following short example illustrates almost the full range of little language syntax::

	// For illustration only, do not try to run
	#define $dens %1.0e18cm-3 // (1) variables and units
	timestep = .01 // (2) floating point assignment

	new plane wave 'pw1' // (3) creating a tool with name given in quotes
	{
		// Several parameter assignments for the plane wave go here.
		// Some parameters are required, some are optional.
	}

	new coulomb field solver 'coul' // (4) creating a module
	{
		dipole center = ( 0.0 , 0.0 , 0.0 ) // (5) tuple assignment
		get 'pw1' // (6) attach the previously defined plane wave tool to this field solver
		new hermite gauss // (7) adding a tool in nested fashion
		{
			// several parameter assignments go here
		}
	}

	new plane wave 'pw2' for 'coul' // (8) simulataneously create and attach a tool
	{
		// several parameter assignments go here
	}

	new species 'electrons'
	{
		// several parameter assignments go here
	}

	generate uniform 'electrons' // (9) create a tool and attach to named module
	{
		density = $dens // (10) using the C-style macro as a numerical constant
	}

If you are familiar with C++ syntax you will recognize comments as preceded by ``//``.  The C-style ``/*`` and ``*/`` comment delimiters are supported also.  Numbered highlights in the above example are as follows:

 	#. C-style preprocessor macros are used to gain the effect of user defined constants.  Numbers can be given physical units using the % prefix followed by a unit specifier postfix.
	#. Simple assignment to a floating point parameter. Since no units are given, it is assumed dimensionless.
	#. Creating a plane wave tool with the user assigned name 'pw1'.  This name can be used later in the input file.
	#. Creating a field solver module; modules do the high level management of data.
	#. Assignment to a tuple, in this case a spatial 3-vector.  The strict syntax requires parenthesis, although the current parser ignores them.
	#. Attaching the previously defined plane wave to the field solver by name.
	#. Attaching a Hermite-Gauss wave tool to the field solver.  Since the tool is nested inside the module, the association is established automatically.  The user-assigned name is optional in this case.
	#. Create another plane wave tool and use the **for** keyword to attach it to the previously defined module.
	#. Create a uniform profile tool and associate with the species 'electrons'.  This illustrates the **generate** keyword, which is merely syntactic sugar that can be used in place of the **new**-**for** combination.
	#. The macro key ``$dens`` is used to achieve the effect of a user-defined constant.

A Little More Little Language
-----------------------------

The input file is a sequence of **directives**.  The directives are one of the following: **preprocessor-directives**, **statements**, or **assignments**.

Preprocessor directives are a limited and slightly modified version of C language preprocessor directives.

Statements start with **new**, **get**, or **generate**.  The **new** and **generate** statements are usually terminated by a curly-brace delimited **block** that may contain a further sequence of directives.  Blocks of directives can be nested in this way with arbitrary depth.

Assignments copy a value in the input file to a simulation parameter.  These values can be decimal numbers, physical quantities with units, or **identifiers**.  Identifiers are used to select from choices, for example, ``periodic`` is an identifier corresponding to a boundary condition.

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

Top level directives may include statements to create modules or tools, as well as assignments to parameters that are associated with the root ``Simulation`` object.  The ``Simulation`` parameter assignments are as follows.

.. py:function:: hardware acceleration device string = dev

	Use hardware accelerators having the given substring in their name

	:param str dev: the substring to search for in the device name, e.g., ``radeon``.  Case doesn't matter.

.. py:function:: hardware acceleration device numbers = dev_list

	Optional specification of preferred OpenCL device numbers.  If specified these take precedence over name search.

	:param list dev_list: variable length list of integers, e.g., ``{ 0 , 1 , 2 }``

.. py:function:: hardware acceleration platform string = platform

	Use only OpenCL platforms having the given substring in their name

	:param str platform: the substring to search for in the platform name, e.g., ``cuda``.  Case doesn't matter.

.. py:function:: unit density = CGS_density

	Sets the unit density and fixes the system of normalized units.

	:param float CGS_density: the density in particles per cubic centimeter.  Unit conversion macros must **not** be used.

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

	:param int dp: Write out checkpoint data every ``dp`` steps.  If zero do not save any checkpoints.

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

Objects (modules and tools) can be created using the following syntax:

.. _block-create:
.. py:function:: new <key1> [<key2> <key3> ...] [<name>] [for <name>] { <directives> }
.. py:function:: generate <key1> [<key2> <key3> ...] [for] <name> { <directives> }

Each form has a preamble followed by a curly-brace delimited block.  The start of the preamble is signaled by a keyword, either ``new`` or ``generate``.  The next several words are ordered keys.  The keys are used to identify the type of object requested.  The user is free to add any number of trailing keys.  In the first form, the first optional name is the user-defined name of the new object, and the second is the name of a previously defined parent object.  Giving a parent object is optional.

The second form allows the new object to be associated with a parent object without using the **for** keyword.  This can be more suggestive in some cases, e.g., ``generate uniform 'electrons'`` is perhaps more suggestive than ``new uniform 'profile' for 'electrons'``.

When optional names are not given, the turboWAVE parser will automatically choose a unique name for the object.

The strict little language syntax requires that the user-assigned name should be in quotes, and that the quoted string satisfy the usual rules for naming an identifier (no white space or special characters except ``_``).  At present the turboWAVE parser does not strictly enforce this, but is likely to do so in the future.

Objects which may have a high multiplicity use a more compact form with ordered directives.  The form is typically

.. py:function:: new <key1> [<key2> <key3> ...] = <directives>

In this case the directives are all required and must be in the right order.

.. _associations:

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

TurboWAVE uses only structured grids, at present.  Some of the effect of unstructured grids can be obtained by using :ref:`grid warps <warps>`.

TurboWAVE axes are labeled as ``x``, ``y``, or ``z`` regardless of coordinate system.  Internally these are often mapped as ``x=1``, ``y=2``, and ``z=3``.  In cylindrical coordinates, ``x`` is radial, ``y`` is azimuthal, and ``z`` is axial.  In spherical coordinates, ``x`` is radial, ``y`` is azimuthal, and ``z`` is polar.

.. py:function:: new grid { directives }

	There must be exactly one grid block, which defines the numerical grid for all modules.

	:param block directives: The following directives are supported:

		.. py:function:: geometry = g

			:param enum g: can be ``cartesian``, ``cylindrical``, ``spherical``.

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

		.. py:function:: adaptive timestep = at

			:param bool at: whether or not to use an adaptive time stepping scheme.

.. _warps:

Grid Warps
----------

Grid warps allow the user to ramp the cell size up or down, along a given axis, and through a given range of cell indices.  Any number of grid warps can be declared as follows:

.. py:function:: new warp { <directives> }

	Ramp the cell sizes along a given axis through the given range of cell indices.  The form of the ramp is a quintic polynomial that has continuous first and second derivatives.

	:param block directives: The following directives are supported:

		.. py:function:: axis = ax

			:param enum ax: The axis along which to create the warp, one of ``x``, ``y``, or ``z``. As usual these are merely labels for whatever coordinate system is in use.

		.. py:function:: increasing = inc

			:param bool inc: If affirmative, the cell size increases with increasing coordinate, otherwise the cell size decreases.

		.. py:function:: index range = ( i0 , i1 )

			:param int i0: cell index where the ramp begins
			:param int i1: cell index where the ramp ends

		.. py:function:: length = L

			:param float L: the length of the ramp
