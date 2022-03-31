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

Apple deprecated OpenCL some years ago, situation is unclear at present.

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

#. Edit :samp:`{turboWAVE}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = CUDA`, and :samp:`COMPILER_PREF = GNU` should be uncommented.
#. Open a new terminal window and navigate to :samp:`{turboWAVE}/core/source`
#. Type :samp:`make clean`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.

AMD GPGPU on RHEL/CentOS 8
==========================

Driver
-------

#. Download the driver for your GPU from the AMD website.
#. Follow the instructions from AMD to install the driver **with OpenCL**.
	* As of this writing, the link is `here <https://amdgpu-install.readthedocs.io>`_.  You must specify a command line option ``--opencl=pal`` or ``--opencl=legacy`` depending on your GPU.  You must also use the ``amdgpu-pro-install`` variant of the install script.

Compile with OpenCL
-------------------

#. Edit :samp:`{turboWAVE}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = RADEON_PRO`, and the compiler preference, should be uncommented.
#. Open a new terminal window and navigate to :samp:`{turboWAVE}/core/source`
#. :samp:`make clean`
#. :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.

NVIDIA GPGPU on Ubuntu 20.04
============================

Driver
------

	#. :samp:`sudo apt update`
	#. :samp:`sudo apt install nvidia-opencl-icd nvidia-opencl-dev`
	#. :samp:`sudo apt update`

Compile with OpenCL
-------------------

#. Edit :samp:`{turboWAVE}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = CUDA`, and the compiler preference should be uncommented.
#. Open a new terminal window and navigate to :samp:`{turboWAVE}/core/source`
#. Type :samp:`make clean`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.

AMD GPGPU on Ubuntu 20.04
=========================

Driver
-------

#. Download the driver for your GPU from the AMD website.
#. Follow the instructions from AMD to install the driver **with OpenCL**.
	* As of this writing, the link is `here <https://amdgpu-install.readthedocs.io>`_.  You must specify a command line option ``--opencl=pal`` or ``--opencl=legacy`` depending on your GPU.  You must also use the ``amdgpu-pro-install`` variant of the install script.

Compile with OpenCL
-------------------

#. Edit :samp:`{turboWAVE}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = RADEON_PRO`, and the compiler preference, should be uncommented.
#. Open a new terminal window and navigate to :samp:`{turboWAVE}/core/source`
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
#. Edit :samp:`{turboWAVE}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = WIN` and :samp:`HARDWARE_ACCEL = CUDA`, should be uncommented.
#. Find ``CL_INCLUDE`` and ``CL_LIB`` in the makefile.  Check the paths against your actual file system and edit them if necessary.
#. Open a new PowerShell and navigate to :samp:`{turboWAVE}/core/source`
#. Type :samp:`make clean`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.

Troubleshooting OpenCL
======================

When OpenCL is not working there is a good chance it has to do with the ability of the software to find the chain of files used to interface with the device driver. The chain looks something like this::

	ICD Loader-->ICD Registry-->ICD File-->Device

For native OS applications, these files have homes in the root directory tree.  For software that runs in an Anaconda environment, they also have homes in the Anaconda directory tree.  A typical troubleshooting procedure is the see if these files are where they should be, and if not, manually copy them to the right location.

Installable Client Driver (ICD) Loader
--------------------------------------

**LINUX** - The ICD loader is a library file that interfaces with device driver libraries in an indirect way.  It is designed specifically to interface with widely varying device types from multiple vendors.  This is the library that is directly called by OpenCL enabled software.  The full path is typically :samp:`/usr/lib{64}/libOpenCL.so`.  In an Anaconda environment, this becomes :samp:`{anaconda}/envs/{my_env}/lib{64}/libOpenCL.so`.  The ICD loader will probably be duplicated somewhere in the vendor's directory tree, e.g. :samp:`/usr/local/cuda` or :samp:`/opt/amdgpu-pro`.  It may help to set the environment variable :samp:`LD_LIBRARY_PATH` to the path of any directory where the ICD loader is found.

**WINDOWS** - In the case of windows, ``opencl.dll`` should be in ``C:\Windows\System32`` or ``C:\Windows\SysWOW64``.

ICD Registry Files
--------------------------------------

**LINUX** - The ICD loader looks in a specific place for "registry" files containing the locations of the ICD libraries themselves (we adopted the terminology used `here <https://wiki.tiker.net/OpenCLHowTo#Installation>`_). The full path is typically :samp:`/etc/OpenCL/vendors/{specific_name}.icd`.  In an Anaconda environment, this becomes :samp:`{anaconda}/envs/{my_env}/etc/OpenCL/vendors/{specific_name}.icd`.  The contents of the registry files are readable ASCII strings with the path of the ICD file from the vendor.  Sometimes only the name of the file (without the path) is given.  If you are having problems you will want to navigate to :samp:`/etc/OpenCL/vendors/`, verify that the registry files are present, and type :samp:`cat *` to print the names of the ICD files.

**WINDOWS** - In the case of windows, the ICD registry is in the system registry.  You can use ``Registry Editor`` to check to see if appropriate entries are present.  There should be a key for each vendor in ``Computer\HKEY_LOCAL_MACHINE\SOFTWARE\Khronos\OpenCL\Vendors``. The name of the key should be the path to a particular vendor's ICD File (see below).

ICD Files
--------------------------------------

**LINUX** - The ICD files are the specific OpenCL implementation from a given vendor.  These are libraries that actually know how to interact with a specific device set. If the full path is given in the ICD registry file, then that is the location of the ICD file (if registry files in different locations point to different places, making these consistent may be the solution).  If only the name is given, then :samp:`/usr/lib{64}/` is a likely place the ICD loader will try.  For Anaconda, the loader might try :samp:`{anaconda}/envs/{my_env}/lib{64}/`.

**WINDOWS** - In the case of windows, the full path should be in the system registry (see above). The path is typically ugly, but the filename should be something like ``IntelOpenCL64.dll`` or ``amdocl64.dll``.

Device Drivers
--------------------------------------

Device drivers are not to be manipulated manually, but if you are curious about locations and names, you can try :samp:`lsmod` to get a list of kernel modules.  Look for something relating to your graphics card (e.g., :samp:`radeon`) and use :samp:`modinfo {radeon}` to see its properties.
