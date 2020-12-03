GPGPU Acceleration
//////////////////

General Notes
=============

As of this writing, GPGPU support is only useful for Quantum Optics modules.
To offload computations to a GPGPU, turboWAVE uses OpenCL.  Typically OpenCL can be used for any modern GPU from either AMD or NVIDIA.  In the case of NVIDIA devices, the OpenCL support is packaged with CUDA.

.. warning::

	Enabling GPGPU acceleration involves manipulating display drivers, which on Linux can be treacherous.  If a Linux display manager is lost, it can be difficult to recover without a full OS reinstallation.  If you cannot afford for this to happen, you should take steps to backup your system configuration.

To run on GPGPU you must prepare a special executable.  The procedure for several operating systems is below.  Once the executable is prepared, running is almost the same as it is for a CPU.  The primary differences are as follows:

	#. Use a single MPI processes per GPGPU device, i.e., if your system has dual video cards you can use 2 MPI processes.
	#. OpenMP threads cannot be used.  The number of OpenMP threads must be one.
	#. If you want to control the particular OpenCL platform and device, use the command line arguments.  Otherwise turboWAVE will select the first available.

GPGPU on MacOS
==============

Driver
------

For MacOS there is nothing to do, OpenCL is supported out of the box.

Compile with OpenCL
-------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = OSX` and :samp:`HARDWARE_ACCEL = APPLE_CL`, should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make clean`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.

NVIDIA GPGPU on RHEL/CentOS 8
=============================

.. Warning::

	The following is tested only for CentOS.

Driver
------

#. Prepare for EPEL:

	* For CentOS, type ``sudo dnf config-manager --set-enabled PowerTools``
	* For RHEL, type ``ARCH=$( /bin/arch )`` followed by ``sudo subscription-manager repos --enable "codeready-builder-for-rhel-8-${ARCH}-rpms"``

#. Go to `EPEL <https://fedoraproject.org/wiki/EPEL>`_ and install.  As of this writing there is a link, ``epel-release-latest-8``, that runs a graphical installer.
#. Go to `RPM Fusion <https://rpmfusion.org/Configuration>`_ and install the ``nonfree`` repository for RHEL 8 or compatible (there is no charge, ``nonfree`` refers to license restrictions).  As of this writing there is a link to run a graphical installer.
#. ``sudo dnf install akmod-nvidia``

	* This automatic kernel module recompiles automatically when a new Linux kernel is installed (e.g. during a system update).  After restarting you must allow extra time for the kernel module to compile.  There could be a long delay before the login screen appears.

#. Restart the system, allow extra time for this restart.
#. ``sudo dnf install xorg-x11-drv-nvidia-cuda``
#. ``sudo dnf install ocl-icd-devel``

Compile with OpenCL
--------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = CUDA`, and :samp:`COMPILER_PREF = GNU` should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make clean`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.

AMD GPGPU on RHEL/CentOS 8
==========================

Driver
-------

#. Install AMD ROCm

	* Perform internet search to find the installation instructions and carry out.  As of this writing the simplest way appears to be to use the script ``rocminstall.py``, see `<https://github.com/srinivamd/rocminstaller>`_.
	* Be sure to test the installation per the installation instructions.
	* This may involve multiple restarts.

#. Create a symbolic link to the ROCm installation

	* :samp:`cd /opt && ls`
	* The output should include a directory in the form :samp:`rocm-{x.y.z}`.
	* If there is no symbolic link :samp:`rocm` create it using :samp:`sudo ln -s rocm-{x.y.z} rocm`

Compile with OpenCL
-------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = RADEON_PRO`, and the compiler preference, should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. :samp:`make clean`
#. :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.

NVIDIA GPGPU on Ubuntu 20.04
============================

Driver
------

	#. :samp:`sudo apt update`
	#. :samp:`sudo apt install nvidia-driver-{XXX} libnvidia-compute-{XXX} nvidia-opencl-dev`

		* Replace :samp:`{XXX}` with version number, e.g., 450

	#. :samp:`sudo apt update`

Compile with OpenCL
-------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = CUDA`, and the compiler preference should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make clean`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.

AMD GPGPU on Ubuntu 20.04
=========================

Driver
-------

#. Install AMD ROCm

	* Perform internet search to find the installation instructions and carry out.  As of this writing the simplest way appears to be to use the script ``rocminstall.py``, see `<https://github.com/srinivamd/rocminstaller>`_.
	* Be sure to test the installation per the installation instructions.
	* This may involve multiple restarts.

#. Create a symbolic link to the ROCm installation

	* :samp:`cd /opt && ls`
	* The output should include a directory in the form :samp:`rocm-{x.y.z}`.
	* If there is no symbolic link :samp:`rocm` create it using :samp:`sudo ln -s rocm-{x.y.z} rocm`

Compile with OpenCL
-------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = RADEON_PRO`, and the compiler preference, should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make clean`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.

GPGPU on Windows 10
===================

Driver
------

Update to the latest graphics drivers, following the guidance on the vendor's website.

Compile for NVIDIA using LLVM
-----------------------------

#. Install LLVM if not already done.
#. Install NVIDIA CUDA SDK
#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = WIN` and :samp:`HARDWARE_ACCEL = CUDA`, should be uncommented.
#. Find ``CL_INCLUDE`` and ``CL_LIB`` in the makefile.  Check the paths against your actual file system and edit them if necessary.
#. Open a new PowerShell and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make clean`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.
