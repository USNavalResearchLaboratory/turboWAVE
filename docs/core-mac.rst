Core Install for MacOS 10.13
----------------------------

As of this writing, the default compiler for MacOS does not support OpenMP sufficiently for turboWAVE.  If you do not plan to use OpenMP, you can use the simple install, which avoids having to deal with Homebrew or MacPorts.

In general, beware of the effects of major OS upgrades.  After each upgrade you may need to re-install XCode.  Remember that the XCode command line tools must be installed in a separate step.  If you are using MacPorts or Homebrew for OpenMP, you should check their respective web sites for information on updating for the new OS.

Simple Install
,,,,,,,,,,,,,,

.. sidebar:: bash profiles

	MacOS treats :file:`~/.bash_profile` differently from linux. Linux runs it only upon login at the OS level.  MacOS runs it for each new interactive terminal session.  That is why for linux we edit :file:`~/.bashrc` rather than :file:`~/.bash_profile`.

#. Put the turboWAVE components (:samp:`core` and :samp:`tools`) into some directory, denoted :samp:`{twroot}`.
#. Install XCode from the App Store
#. Install the XCode command line tools.  Go to the terminal and type :samp:`xcode-select --install`.  Respond affirmatively to the prompts.
#. :samp:`cd ~`
#. :samp:`mkdir bin`
#. :samp:`mkdir Run`
#. Verify bash as the login shell

	* :samp:`echo $0` prints the currently running shell
	* :samp:`chsh -s /bin/bash` changes the default to bash

#. Edit :samp:`~/.bash_profile`, adding the line :samp:`export PATH=~/bin/:$PATH` (no spaces around equals sign)
#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = OSX` should be uncommented.
#. Edit :samp:`{twroot}/core/source/definitions.h`
#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`, and :samp:`#` is **not** a comment.  For this installation, only :samp:`#define USE_DESKTOP` should be uncommented.  Comment out :samp:`#define USE_OPENMP` and :samp:`#define USE_OPENCL` if they are not already.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.

Advanced Install with Homebrew
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

#. Put the turboWAVE components (:samp:`core` and :samp:`tools`) into some directory, denoted :samp:`{twroot}`.
#. Install XCode from the App Store
#. Install the XCode command line tools.  Go to the terminal and type :samp:`xcode-select --install`.  Respond affirmatively to the prompts.
#. Perform internet search to find Homebrew installation instructions and carry out
#. In the terminal type :samp:`brew install llvm`
#. :samp:`cd ~`
#. :samp:`mkdir bin`
#. :samp:`mkdir Run`
#. Verify bash as the login shell

	* :samp:`echo $0` prints the currently running shell
	* :samp:`chsh -s /bin/bash` changes the default to bash

#. Edit :samp:`~/.bash_profile`, adding the line :samp:`export PATH=~/bin/:$PATH`
#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = OSX`, :samp:`HARDWARE_ACCEL = OMP`, :samp:`COMPILER_PREF = LLVM_CLANG`, and :samp:`PACKAGE_PREF = HOMEBREW` should be uncommented.
#. Edit :samp:`{twroot}/core/source/definitions.h`
#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`.  For this installation, only :samp:`#define USE_DESKTOP` and :samp:`#define USE_OPENMP` should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.
#. If you use :samp:`vim`, you may want to copy :samp:`{twroot}/tools/config-files/filetype.vim` to :samp:`~/.vim`.  This will enable syntax highlighting while editing turboWAVE input files.

Advanced Install with MacPorts
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

#. Put the turboWAVE components (:samp:`core` and :samp:`tools`) into some directory, denoted :samp:`{twroot}`.
#. Install XCode from the App Store
#. Install the XCode command line tools.  Go to the terminal and type :samp:`xcode-select --install`.  Respond affirmatively to the prompts.
#. Perform internet search to find MacPorts installation instructions and carry out
#. In the terminal type :samp:`port selfupdate`
#. :samp:`sudo port install llvm-5.0`
#. :samp:`sudo port install clang-5.0`
#. :samp:`sudo port install clang_select`
#. :samp:`sudo port install libcxx`
#. :samp:`sudo port select clang mp-clang-5.0`
#. :samp:`cd ~`
#. :samp:`mkdir bin`
#. :samp:`mkdir Run`
#. Verify bash as the login shell

	* :samp:`echo $0` prints the currently running shell
	* :samp:`chsh -s /bin/bash` changes the default to bash

#. Edit :samp:`~/.bash_profile`, adding the line :samp:`export PATH=~/bin/:$PATH`
#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = OSX`, :samp:`HARDWARE_ACCEL = OMP`, :samp:`COMPILER_PREF = LLVM_CLANG`, and :samp:`PACKAGE_PREF = MACPORTS` should be uncommented.
#. Edit :samp:`{twroot}/core/source/definitions.h`
#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`.  For this installation, only :samp:`#define USE_DESKTOP` and :samp:`#define USE_OPENMP` should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.
#. If you use :samp:`vim`, you may want to copy :samp:`{twroot}/tools/config-files/filetype.vim` to :samp:`~/.vim`.  This will enable syntax highlighting while editing turboWAVE input files.
