GPGPU on Windows 10
===================

.. caution::

	We assume core turboWAVE **and** tools have already been installed according to the documentation, with no steps omitted.

Driver
------

Update to the latest graphics drivers, following the guidance on the vendor's website.

Compile for NVIDIA using LLVM
-----------------------------

#. Install LLVM if not already done, see :doc:`core-windows`.
#. Install NVIDIA CUDA SDK
#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = WIN` and :samp:`HARDWARE_ACCEL = CUDA`, should be uncommented.
#. Find ``CL_INCLUDE`` and ``CL_LIB`` in the makefile.  Check the paths against your actual file system and edit them if necessary.
#. Open a new PowerShell and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make clean`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.
