What's New
//////////

What's New in 4.5
=================

#. Addition of units helpful for hydro, namely specific energy, mass density, and pressure in bars.

#. Internally, dimensional number class with operator based conversions streamlines notation.

Gotchas
-------

Old dimensional number specifiers are no longer allowed. Attempting to use them will result in ``unrecognized units`` or ``unexpected directive``.

What's New in 4.4
=================

#. Input file support for a subset of C-style preprocessor conditionals (``#ifdef``, ``#ifndef``, ``#else``, and ``#endif``).

#. Script for carrying out parallel-serial or parallel-parallel parameter scans is included in ``twutils`` package.

#. Support for additional ionization models.

#. Documentation for ``twutils`` modules.

What's New in 4.3
=================

Core configuration and installation is done using an installer program that is packaged with tools.  The installer uses a ``tk`` graphical interface by default, but can also use a textual interface (based on ``curses``) if needed for remote work.  The old installation workflow can still be used as well.

Gotchas
-------

#. Install workflow is reversed relative to previous versions.  First install tools to gain access to the installer.  Then use the installer to configure and install core.

What's New in 4.2
=================

The turboWAVE tools are available on both PyPi and conda.  Installation is simplified and certain scripts are made available as part of the environment.

Gotchas
-------

#. Scripts moved from ``extras`` to the ``twutils`` package.  Invocation is simpler, e.g., run ``twplot ...`` rather than ``python path/to/twplot.py ...`` (from within your conda or venv environment).

#. The ``twtest`` script also moved into the ``twutils`` package.

What's New in 4.1
=================

Python Friendly Outputs
-----------------------

Binaries are now written as ``numpy`` arrays using the standard ``.npy`` format. These can be read into your Python scripts directly, using ``numpy.load``.  The old turboWAVE specific ``.dvdat`` format is retired.

Metadata (units, axis information, labels, etc.) associated with each binary file is kept in a single file, ``tw_metadata.json``.  This can be easily loaded into a Python dictionary.

Further Improvements in Units
-----------------------------

TurboWAVE's internal comprehension of several systems of units is more fully integrated.  Unit conversion is taken out of the preprocessor and moved into the parser for better semantic control.  Input files can now specify a native system of units, subject to possible veto by individual modules.  This improvement allows information on the units to be reliably incorporated into the metadata, and allows visualization tools to perform unit translations very easily.

Better Jupyter DataViewer
-------------------------

The Jupyter DataViewer is greatly improved.  Axes are labeled clearly and units can be selected on the fly.  Animation is smoother, and interactions through ``ipympl`` are supported.  The native MacOS and Windows DataViewers are retired.

Gotchas
-------

#. Since the ``.dvdat`` format is retired, the CLI plotter ``plot-dvdat.py`` is renamed ``twplot.py``.  Similarly ``maya-dvdat.py`` is renamed ``twmaya.py``.

#. The metadata file ``tw_metadata.json`` and the grid files must be kept with the simulation outputs in order for post-processing tools to perform optimally.

#. Your Python environment needs to add the ``ipympl`` package for the Jupyter DataViewer.  Also don't forget to install the new ``twutils`` package.

#. Using explicitly dimensional numbers (unit macros) in SPARC reactions/collisions is still forbidden, although the latest improvements open a path to removing this exception.

What's new in 4.0
=================

TurboWAVE 4.0 is a major upgrade with many internal improvements and modernizations.  Internally it is streamlined, using 10% fewer lines of code while providing enhanced functionality.

Grid Warps
----------

The system of modulating the cell size is generalized to support arbitrary numbers of :ref:`warps` along any axis.  Each warp defines an upramp or downramp in the cell size through a specified range of cell indices along the given axis.  The form of the ramp is the usual :math:`{\cal C}^2` quintic polynomial.  N.b. grid warps are useful for SPARC hydro modules, not PIC.

Physical Units
--------------

Physical units are treated in a consistent way throughout.  Units are specified by the user in a simple intuitive way.  For example, entering ``%5cm`` is read by the parser as five centimeters, while ``%10deg`` is understood as an angular dimension in degrees.  If a raw decimal number is given, the quantity is assumed normalized, e.g. lengths in units of :math:`c/\omega_p` or angles in units of radians.

.. Note::

  The exception is the SPARC chemistry database, which still assumes the CGS-eV system (and likely always will).  Unit conversion macros should **not** be used when creating chemical reactions.

Comprehensive Error Checking
----------------------------

The tuboWAVE parser is much more sophisticated, providing useful error messages for almost any input file error.  This functionality is also provided for free to developers of ``Module`` and ``ComputeTool`` objects.  Internally, the code to setup input file interactions is much more streamlined.

The groundwork is laid for a professional syntax checker to be incorporated which can pinpoint the line and column of a syntax error.

Input Files
-----------

the TurboWAVE input file now has a strict language definition, and all internal objects conform to consistent semantics.  This enhances the scope and predictability of relationships the user can create between modules and tools.

Better Diagnostics
------------------

Standard input file semantics allow any diagnostic to be associated (or not) with any number of modules.  Storage can be saved by directing the box diagnostic to save only fields of interest.  The phase space diagnostic is more versatile, supporting up to three dimensions, and twelve possible axes.  Internally the system encourages the development of new, sophisticated diagnostic modules.

C++11 and C++17
------------------

TurboWAVE started as a C++98 code. We have been gradually incorporating C++11 style coding.  With version 4.0 the code is solidly C++11.  Improved special function support comes from C++17, although internal special functions are still kept around until compiler support is more consistent.

.. Note::

  We have our eye on C++20

Gotchas for 3.x Users
---------------------

#. First the good news, improved error checking will help you correct most input file errors.

#. Various directives are changed or retired.  If you get ``Unexpected directive`` or ``keys were not understood`` you must consult this documentation, or the examples, and find the appropriate replacement.  The parser is, by design, not as forgiving as before.

#. Now that units are treated in a consistent way, version 3.x input files, which have inconsistent treatment of units, can silently break.

	* Angular dimensions are radians by default, if you want degrees use a dimensional number, e.g., ``%45deg``.
	* Diffusivity units in SPARC are normalized by default, if you want to use dimensional numbers, you must do so explicitly, e.g., put ``%1.0cm2s``.
	* Due to their high multiplicity, SPARC reactions and collisions are an exception.  The raw numbers are expected to be in CGS-eV and always will be.  **Your collisions need to be edited** because in version 3.x the cross section was taken as normalized.  Put it in CGS.  **Do not use unit conversion macros** in reactions or collisions.


#. The ``open`` keyword for reading checkpoint data is retired.  To restart a simulation use the command line argument ``--restart`` and leave the input file the same.  The ``dump period`` parameter works the same as before.

#. The initial condition gets written to diagnostic files as the first frame, so there is typically one extra diagnostic frame relative to before.  This awareness is all you really need.

	* In detail, diagnostics are written at the beginning of a step just as before.  the first step is now numbered as step 0, which causes the diagnostic write-out evaluation to always be initially true.  In order to get the last step written out in the expected way, turboWAVE will take one extra step at the end, i.e., if you ask for n steps, the actual number of steps is n+1 (you will see this on ``stdout``).  For checkpointing, the ``dump period`` should still be some integer factor of n, the restart mechanism is aware of the extra step and will take care of everything.

#. The filename ``full`` is no longer treated specially.  If you want to eliminate the prefix on box diagnostic files simply do not assign a filename.  Trying to do this for more than one box throws an error.

#. In order to have consistent semantics, the syntax for injecting radiation needed to be slightly changed.  In brief, the radiation tool has to be explicitly associated with a module.  Please see :ref:`associations` and :ref:`radiation`.

#. The syntax for phase space diagnostics is changed, see :ref:`specific-diagnostics`.

#. Ionization models are now attached to modules using tools.  See :ref:`ionization` and the examples (search for ``ionization`` in examples folder).

#. Equation of state tools use standard syntax and semantics, see :ref:`eos`.

#. OpenCL platforms and devices are specified on the command line rather than in the input file.
