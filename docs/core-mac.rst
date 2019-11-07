Core Install for MacOS 10.15
----------------------------

As of this writing, the default compiler for MacOS does not support OpenMP sufficiently for turboWAVE.  If you do not plan to use OpenMP, you can use the simple install, which avoids having to deal with Homebrew or MacPorts.

In general, beware of the effects of major OS upgrades.  After each upgrade you may need to re-install XCode.  Remember that the XCode command line tools must be installed in a separate step.  If you are using MacPorts or Homebrew for OpenMP, you should check their respective web sites for information on updating for the new OS.

.. _simple-install:

Simple Install
,,,,,,,,,,,,,,

#. If not already done, get the turboWAVE components, see :doc:`getting-components`. You should end up with a new directory containing at least ``core`` and ``tools``.  This directory can be renamed if desired.  We refer to it generically as :samp:`{twroot}` throughout this documentation.
#. Install XCode from the App Store
#. Install the XCode command line tools.  Go to the terminal and type :samp:`xcode-select --install`.  Respond affirmatively to the prompts.
#. :samp:`cd ~`
#. :samp:`mkdir bin`
#. :samp:`mkdir Run`
#. Verify ``zsh`` as the login shell

	* :samp:`echo $0` prints the currently running shell
	* :samp:`chsh -s /bin/zsh` changes the default to ``zsh``
	* Terminal preferences can also affect the default shell

#. Edit :samp:`~/.zshrc`, adding the line :samp:`export PATH=~/bin/:$PATH` (no spaces around equals sign)
#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = OSX` should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.

Advanced Install with Homebrew
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

#. Perform the :ref:`simple-install` steps.
#. Perform internet search to find Homebrew installation instructions and carry out.
#. In the terminal type :samp:`brew update`
#. In the terminal type :samp:`brew install llvm`
#. In the makefile, uncomment :samp:`HARDWARE_ACCEL = OMP` and :samp:`PACKAGE_PREF = HOMEBREW`.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.

Advanced Install with MacPorts
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

#. Perform the :ref:`simple-install` steps.
#. Perform internet search to find MacPorts installation instructions and carry out
#. In the terminal type :samp:`sudo port selfupdate`
#. Use :samp:`port search llvm` to find the latest version of LLVM. In the following this is denoted :samp:`{X}.0`.
#. :samp:`sudo port install llvm-{X}.0 clang-{X}.0`
#. :samp:`sudo port select clang mp-clang-{X}.0`
#. In the makefile, uncomment :samp:`HARDWARE_ACCEL = OMP` and :samp:`PACKAGE_PREF = MACPORTS`.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.
