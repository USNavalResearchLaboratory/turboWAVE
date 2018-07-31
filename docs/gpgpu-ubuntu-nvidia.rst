NVIDIA GPGPU on Ubuntu 16.04
============================

.. caution::

	We assume core turboWAVE **and** tools have already been installed according to the documentation, with no steps omitted.

Driver
------

Graphics drivers can change rapidly, so internet searches may figure prominently into your installation effort.
NVIDIA GPGPU devices use CUDA software, which includes OpenCL.  This typically supports NVIDIA devices only.

#. Find the Debian local install file for NVIDIA CUDA (use internet search)
#. OS=linux, Arch=x86_64, Dist=Ubuntu, Vers=16.04, Installer=deb(local)
#. Do NOT use the runfile
#. Navigate to downloaded file in the terminal

	- :samp:`sudo dpkg -i {downloaded_file}`
	- :samp:`sudo apt update`
	- :samp:`sudo apt install cuda`
	- :samp:`sudo apt update`

#. The CUDA files should be in :samp:`/usr/local/cuda/`.  At the time of writing the paths below are correct, but NVIDIA can change them.  For the headers, you need the directory containing the :samp:`CL` directory.  For the libraries, you need the one containing :samp:`libOpenCL.so`.
#. Edit :samp:`~/.bashrc`

	- Add line :samp:`export PATH=/usr/local/cuda/bin:$PATH`
	- Add line :samp:`export CPATH=/usr/local/cuda/include:$CPATH`
	- Add line :samp:`export LIBRARY_PATH=/usr/local/cuda/lib64:$LIBRARY_PATH`
	- Add line :samp:`export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LIBRARY_PATH`

#. Remember, no spaces around equals sign.
#. Close all terminal windows.


Compile with OpenCL
-------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = CUDA`, and the compiler preference should be uncommented.
#. Edit :samp:`{twroot}/core/source/definitions.h`
#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`, and :samp:`#` is **not** a comment.  For this installation, only :samp:`#define USE_DESKTOP` and :samp:`#define USE_OPENCL` should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.
