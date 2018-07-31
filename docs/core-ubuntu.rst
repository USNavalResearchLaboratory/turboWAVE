Core Install for Ubuntu 16.04
=============================

As of this writing, the default LLVM-clang compiler on Ubuntu fails to compile due to a bug in OpenMP support.  One must explicitly install LLVM 5.0 to resolve the issue.  One can choose to install either GCC, LLVM, or both.

Install GCC
-----------

#. Open a terminal window
#. :samp:`sudo apt update`
#. :samp:`sudo apt install g++`
#. :samp:`sudo apt install libomp-dev`

Install LLVM 5.0
----------------

#. Open a terminal window
#. :samp:`sudo apt update`
#. :samp:`sudo apt install llvm-5.0`
#. :samp:`sudo apt install clang-5.0`
#. :samp:`sudo apt install libc++-dev`

Configure File System
---------------------

#. Put the turboWAVE components (:samp:`core` and :samp:`tools`) into some directory, denoted :samp:`{twroot}`.
#. Open a terminal window
#. :samp:`cd ~`
#. :samp:`mkdir bin`
#. :samp:`mkdir Run`
#. Verify bash as the login shell

	* :samp:`echo $0` prints the currently running shell
	* :samp:`chsh -s /bin/bash` changes the default to bash

#. Edit :samp:`~/.bashrc`, adding the line :samp:`export PATH=/usr/lib/llvm-5.0/bin/:$PATH`
#. Edit :samp:`~/.bashrc`, adding the line :samp:`export PATH=~/bin/:$PATH`
#. When exporting variables in bash, do not put spaces around the equals sign.
#. Close all terminal windows.

Compile turboWAVE
-----------------

	#. Edit :samp:`{twroot}/core/source/makefile`
	#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = OMP`, and :samp:`COMPILER_PREF = LLVM_CLANG`, should be uncommented.  You may substitute :samp:`GNU` for :samp:`LLVM_CLANG`, per your preference.
	#. Edit :samp:`{twroot}/core/source/definitions.h`
	#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`.  For this installation, only :samp:`#define USE_DESKTOP` and :samp:`#define USE_OPENMP` should be uncommented.
	#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
	#. Type :samp:`make`
	#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.
	#. If you use :samp:`vim`, you may want to copy :samp:`{twroot}/tools/config-files/filetype.vim` to :samp:`~/.vim`.  This will enable syntax highlighting while editing turboWAVE input files.
