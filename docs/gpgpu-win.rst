GPGPU on Windows 10
===================

.. caution::

	We assume core turboWAVE **and** tools have already been installed according to the documentation, with no steps omitted.

Driver
------

Update to the latest graphics drivers, following the guidance on the vendor's website.

Compile with OpenCL
-------------------

#. Edit :samp:`{twroot}\\core\\source\\winmakefile`
#. Change to makefile to support OpenCL compilation (TBD)
#. Edit :samp:`{twroot}/core/source/definitions.h`
#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`, and :samp:`#` is **not** a comment.  For this installation, only :samp:`#define USE_DESKTOP` and :samp:`#define USE_OPENCL` should be uncommented.
#. From the developer prompt navigate to :samp:`{twroot}`:samp:`\\core\\source`
#. :samp:`nmake /F winmakefile`
#. Copy :samp:`tw3d.exe` and the OpenCL kernel files (``*.cl``) to :samp:`{Run}`.
