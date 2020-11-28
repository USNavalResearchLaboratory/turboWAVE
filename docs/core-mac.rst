Core Install for MacOS 10.15
----------------------------

As of this writing, the default compiler for MacOS does not support OpenMP sufficiently for turboWAVE. Furthermore advanced special functions are not packaged with LLVM.  For this reason we recommend using GCC on Mac OS at present.  This is installed using either the Homebrew or MacPorts package manager.

In general, beware of the effects of software updates.  When XCode is updated it may be necessary to rebuild the command line tools, and update the compiler.  Remember that the XCode command line tools must be installed in a separate step.  After a major OS upgrade, check the package manager's web site for information on updating for the new OS.

.. _mac-prelim:

Preliminary Steps
,,,,,,,,,,,,,,,,,

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

Install with Homebrew
,,,,,,,,,,,,,,,,,,,,,

#. Perform the :ref:`mac-prelim`.
#. Perform internet search to find Homebrew installation instructions and carry out.
#. In the terminal type :samp:`brew update`
#. In the terminal type :samp:`brew install gcc`
#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = OSX`, :samp:`HARDWARE_ACCEL = OMP`, :samp:`COMPILER_PREF = GNU`, and :samp:`PACKAGE_PREF = HOMEBREW` should be uncommented.
#. In the makefile, set the constant ``VBITS`` to match the width of the available vector extensions, using the guidance in the nearby comments.  You can search the output of :samp:`sysctl machdep.cpu` for the available vector extensions.  If still in doubt set it to 256.
#. The ``g++`` version number is hard-coded in the makefile, e.g., as ``g++-10``.  If you have a different version you will have to edit this.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.

Install with MacPorts
,,,,,,,,,,,,,,,,,,,,,

#. Perform the :ref:`mac-prelim`.
#. Perform internet search to find MacPorts installation instructions and carry out
#. In the terminal type :samp:`sudo port selfupdate`
#. :samp:`sudo port install gcc10`
#. :samp:`sudo port select --set gcc mp-gcc10`
#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = OSX`, :samp:`HARDWARE_ACCEL = OMP`, :samp:`COMPILER_PREF = GNU`, and :samp:`PACKAGE_PREF = MACPORTS` should be uncommented.
#. In the makefile, set the constant ``VBITS`` to match the width of the available vector extensions, using the guidance in the nearby comments.  You can search the output of :samp:`sysctl machdep.cpu` for the available vector extensions.  If still in doubt set it to 256.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.
