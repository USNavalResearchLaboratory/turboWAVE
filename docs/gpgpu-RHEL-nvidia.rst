NVIDIA GPGPU on RHEL 7.5
========================

.. caution::

	We assume core turboWAVE **and** tools have already been installed according to the documentation, with no steps omitted.

Driver
------

#. :samp:`sudo yum update`
#. :samp:`sudo yum install ocl-icd clinfo`
#. Perform internet search to find instructions for installing ``ELRepo``, and carry out.
#. :samp:`sudo yum install kmod-nvidia`
#. Reboot the system
#. :samp:`clinfo` should give a listing of platforms and devices, if the installation succeeded.
#. If there is a problem compiling (see below), you may need to explicitly create a symbolic link as follows:

	* :samp:`sudo ln -s /usr/lib64/libOpenCL.so.1.0.0 /usr/lib64/libOpenCL.so`

#. Close all terminal windows

Compile with OpenCL
--------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = CUDA`, and :samp:`COMPILER_PREF = GNU` should be uncommented.
#. Edit :samp:`{twroot}/core/source/definitions.h`
#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`, and :samp:`#` is **not** a comment.  For this installation, only :samp:`#define USE_DESKTOP` and :samp:`#define USE_OPENCL` should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`scl enable devtoolset-7 'make'`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.
