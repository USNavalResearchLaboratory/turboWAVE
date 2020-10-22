Core Install for Ubuntu 20.04
=============================

Install GCC
-----------

#. Open a terminal window
#. :samp:`sudo apt update`
#. :samp:`sudo apt install g++ libomp-dev`

Install LLVM
------------

#. Open a terminal window
#. :samp:`sudo apt update`
#. :samp:`sudo apt install llvm clang`

Configure File System
---------------------

#. Get the turboWAVE components, see :doc:`getting-components`. You should end up with a new directory containing at least ``core`` and ``tools``.  This directory can be renamed if desired.  We refer to it generically as :samp:`{twroot}` throughout this documentation.
#. Open a terminal window
#. :samp:`cd ~`
#. :samp:`mkdir bin`
#. :samp:`mkdir Run`
#. Verify bash as the login shell

	* :samp:`echo $0` prints the currently running shell
	* :samp:`chsh -s /bin/bash` changes the default to bash

#. Edit :samp:`~/.bashrc`, adding the line :samp:`export PATH=~/bin/:$PATH`
#. When exporting variables in bash, do not put spaces around the equals sign.
#. Close all terminal windows.

Compile turboWAVE
-----------------

	#. Edit :samp:`{twroot}/core/source/makefile`
	#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = OMP`, and :samp:`COMPILER_PREF = LLVM_CLANG`, should be uncommented.  You may substitute :samp:`GNU` for :samp:`LLVM_CLANG`, per your preference.
	#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
	#. Type :samp:`make`
	#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.
