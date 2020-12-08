Prepare to Install
//////////////////

The turboWAVE installer assumes Anaconda and a modern C++ compiler are available.  If you are sure you have this you can skip to :doc:`install-easy`.  Otherwise find your operating system below and follow the steps.

High Performance Computing (HPC)
================================

On HPC systems, the necessary components should already be prepared.  However, one sometimes has to load appropriate modules using the ``module load`` command.  Look for modules supporting Anaconda and a C++ compiler.  The installer is designed to work with Cray systems.  For others it may be necessary to perform a manual install.

Linux - CentOS 8
================

#. Open a terminal window
#. :samp:`sudo dnf install gcc make libomp-devel`
#. :samp:`sudo dnf install llvm clang`
#. Download Anaconda3 or Miniconda3 from the internet and install
#. :samp:`conda update conda`

For other RPM based distributions the procedure should be the same.  Note Enterprise Linux 7 variants will not work without the developer toolset.

Linux - Ubuntu 20.04
====================

#. Open a terminal window
#. :samp:`sudo apt update`
#. :samp:`sudo apt install g++ make libomp-dev`
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
#. :samp:`sudo port install gcc10`
#. :samp:`sudo port select --set gcc mp-gcc10`

Install Anaconda
----------------

#. Download Anaconda3 or Miniconda3 from the internet and install
#. In the terminal type ``conda update conda``

Windows 10
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
#. Open a new PowerShell window
#. If you get an error you likely need to update the Execution Policy

	* :samp:`Set-ExecutionPolicy Bypass -Scope CurrentUser`
	* Respond affirmatively to the prompt, close and reopen the PowerShell

#. :samp:`conda update conda`
#. :samp:`conda init powershell`

Notes on Text Editors
----------------------

Most turboWAVE text files, such as input file examples, have UNIX line feeds.  This is no problem for WordPad (set word wrap to no wrap), but Notepad may not display them properly.  Installing a developer-oriented text editor (e.g. Atom, Sublime) might be useful.  You can also install terminal-style editors such as ``vim`` for use in the PowerShell::

	choco install vim

Install LLVM (free)
-------------------

#. Open an administrator PowerShell window.
#. :samp:`choco install llvm`
#. :samp:`choco install make`

Install Intel (may require purchase)
------------------------------------

#. Download and install Intel Parallel Studio.

	* The Intel compiler is a commercial product, but you may be able to use it freely on a trial basis.
