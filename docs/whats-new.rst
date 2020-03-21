What's New
============

TurboWAVE 4.0.0 is a major upgrade with many internal improvements, as well as improvements visible to users.

Grid Warps
----------

The system of modulating the cell size is generalized to support arbitrary numbers of "warps" along any axis.  Each warp defines an upramp or downramp in the cell size through a specified range of cell indices along the given axis.  The form of the ramp is the usual :math:`{\cal C}^2` quintic polynomial.  N.b. grid warps are useful for SPARC hydro modules, not PIC.

Physical Units
--------------

Physical units are treated in a consistent way throughout.  Units are specified by the user in a simple intuitive way.  For example, entering ``%5cm`` is read by the parser as five centimeters, while ``%10deg`` is understood as an angular dimension in degrees.  If a raw decimal number is given, the quantity is assumed normalized, e.g. lengths in units of :math:`c/\omega_p` or angles in units of radians.

.. Note::

  The exception is the SPARC chemistry database, which still assumes the CGS-eV system.  Unit conversion macros should **not** be used when creating chemical reactions.

Comprehensive Error Checking
----------------------------

The tuboWAVE parser is much more sophisticated, providing useful error messages for almost any input file error.  This functionality is also provided for free to developers of ``Module`` and ``ComputeTool`` objects.  Internally, the code to setup input file interactions is much more streamlined.

Input Files
-----------

the TurboWAVE input file now has a strict language definition, and all internal objects conform to consistent semantics.  This enhances the scope and predictability of relationships the user can create between modules and tools.

Better Diagnostics
------------------

Diagnostics generally benefit from consistent input file semantics.  Any diagnostic can now be associated (or not) with any number of modules.  The phase space diagnostic is more versatile, supporting up to three dimensions, and twelve possible axes.  This system encourages the development of new, sophisticated diagnostic modules.

C++11 and C++17
------------------

TurboWAVE started as a C++98 code. We have been gradually incorporating C++11 style coding.  With release 4.0.0 the code is solidly C++11.  Improved special function support comes from C++17, although internal special functions are still kept around until compiler support is more consistent.

.. Note::

  We have our eye on C++20
