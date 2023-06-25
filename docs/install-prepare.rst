Prepare to Install
//////////////////

The turboWAVE installer assumes Anaconda and a modern C++ compiler are available.  If you are sure you have this you can skip to :doc:`install-easy`.  Otherwise find your operating system below and follow the steps.

High Performance Computing (HPC)
================================

On HPC systems, the necessary components should already be prepared.  However, one sometimes has to load appropriate modules using the ``module load`` command.  Look for modules supporting Anaconda and a C++ compiler.

.. tip::

	Non-GCC compilers sometimes rely on GCC for the standard library.  If you are trying to load an updated non-GCC compiler, you may need to load GCC as well.

Linux - Ubuntu
==============

#. Open a terminal window
#. :samp:`sudo apt update`
#. :samp:`sudo apt install g++ libomp-dev`
#. :samp:`sudo apt install llvm clang`
#. Download Anaconda3 or Miniconda3 from the internet and install
#. :samp:`conda update conda`

For other Debian based distributions the procedure should be the same.


Mac OS
======

As of this writing, the default compiler for MacOS does not support OpenMP sufficiently for turboWAVE. Furthermore advanced special functions are not packaged with LLVM.  For this reason we recommend using GCC on Mac OS at present.  This is installed using either the Homebrew or MacPorts package manager.

In general, beware of the effects of software updates.  When XCode is updated it may be necessary to rebuild the command line tools, and update the compiler.  Remember that the XCode command line tools must be installed in a separate step.  After a major OS upgrade, check the package manager's web site for information on updating for the new OS.

Install XCode
-------------

#. Install XCode from the App Store
#. Install the XCode command line tools.  Go to the terminal and type :samp:`xcode-select --install`.  Respond affirmatively to the prompts.

Install Homebrew GCC
--------------------

#. If you want to use MacPorts skip this section.
#. Perform internet search to find Homebrew installation instructions and carry out.
#. In the terminal type :samp:`brew update`
#. In the terminal type :samp:`brew install gcc`

Install MacPorts GCC
--------------------

#. If you want to use Homebrew skip this section.
#. Perform internet search to find MacPorts installation instructions and carry out
#. In the terminal type :samp:`sudo port selfupdate`
#. :samp:`sudo port install gcc12`
#. :samp:`sudo port select --set gcc mp-gcc12`

In the above you can substitute a later version for ``12``, if available.

Install Anaconda
----------------

#. Download Anaconda3 or Miniconda3 from the internet and install
#. In the terminal type ``conda update conda``

Windows 11
==========

Visual Studio
-------------

Install the latest Visual Studio.  The Community Edition is free and should suffice.  Select at least the option ``Desktop Development with C++``.

PowerShell Setup
----------------

#. Enter ``powershell`` into the Cortana search field.  You should see the PowerShell as an option.  Right click this and select ``Run as Administrator``.
#. In your internet browser search for Chocolatey and follow the instructions to install it using a PowerShell command.

	* You should be guided through setting up the ExecutionPolicy prior to running the installation command

.. tip::

	The PowerShell supports the use of many UNIX style conventions, such as forward slashes as directory separators, the twiddle as a short-cut for the home directory, and short form commands like ``ls`` and ``cp``.

Install Anaconda
----------------

#. Run Anaconda3 or Miniconda3 installer from internet, accept defaults.
#. Open the special Anaconda PowerShell terminal
#. :samp:`conda update conda`

Install LLVM
------------

#. Open an administrator PowerShell window.
#. :samp:`choco install llvm`
