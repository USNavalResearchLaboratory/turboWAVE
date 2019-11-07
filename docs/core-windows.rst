Core Install for Windows 10
===========================

Visual Studio
-------------

Install the latest Visual Studio.  The Community Edition is free and should suffice.  Select at least the option ``Desktop Development with C++``.

PowerShell Setup
----------------

#. Enter ``powershell`` into the Cortana search field.  You should see the PowerShell as an option.  Right click this and select ``Run as Administrator``.
#. In your internet browser search for Chocolatey and follow the instructions to install it using a PowerShell command.

	* You should be guided through setting up the ExecutionPolicy prior to running the installation command

#. Close the PowerShell window and open a new one, again as administrator.
#. :samp:`choco install git`
#. Close the administrator PowerShell.
#. Open a new PowerShell window as a user (left click).
#. :samp:`mkdir ~/bin`
#. :samp:`mkdir ~/Run`

.. tip::

	The PowerShell supports the use of many UNIX style conventions, such as forward slashes as directory separators, the twiddle as a short-cut for the home directory, and short form commands like ``ls`` and ``cp``.

Notes on Text Editors
----------------------

Most turboWAVE text files, such as input file examples, have UNIX line feeds.  This is no problem for WordPad (set word wrap to no wrap), but Notepad may not display them properly.  Installing a developer-oriented text editor (e.g. Atom, Sublime) might be useful.  You can also install terminal-style editors such as ``vim`` for use in the PowerShell::

	choco install vim

Path Variable
-------------

#. Determine the full pathname of ``~/bin``, e.g., by using File Explorer.  Let us denote this :samp:`{full_path_to_bin}`.
#. Open Windows Settings, and click on the System category.
#. Use the search box to find ``Edit environment variables for your account``.
#. Highlight ``path`` in the user variables area and click ``Edit...``.
#. Click ``New...`` and enter :samp:`{full_path_to_bin}`.
#. Click ``OK`` and dismiss all windows.

Compile with LLVM (free)
------------------------

#. If not already done, get the turboWAVE components, see :doc:`getting-components`. You should end up with a new directory containing at least ``core`` and ``tools``.  This directory can be renamed if desired.  We refer to it generically as :samp:`{twroot}` throughout this documentation.
#. Edit :samp:`{twroot}\\core\\source\\makefile`

	* Uncomment ``PLATFORM = WIN`` and comment out all other platforms. In a makefile, comments are preceded by :samp:`#`.
	* Uncomment ``HARDWARE_ACCEL = OMP`` and comment out all other accelerators.
	* Uncomment ``COMPILER_PREF = LLVM_CLANG`` and comment out all other compilers.

#. Open an administrator PowerShell window.
#. :samp:`choco install llvm`
#. :samp:`choco install make`
#. Open a new user PowerShell window.
#. :samp:`cd` :samp:`{twroot}`:samp:`\\core\\source`
#. :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~\\bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~\\Run`, but these will not be used.

Compile with Intel (may require purchase)
-----------------------------------------

#. If not already done, get the turboWAVE components, see :doc:`getting-components`. You should end up with a new directory containing at least ``core`` and ``tools``.  This directory can be renamed if desired.  We refer to it generically as :samp:`{twroot}` throughout this documentation.
#. Edit :samp:`{twroot}\\core\\source\\win.make`

	* Uncomment ``COMPILER_PREF = INTEL`` and comment out ``COMPILER_PREF = VS``. In a makefile, comments are preceded by :samp:`#`.
	* Uncomment ``CCFLAGS = $(RELEASE_FLAGS)`` and comment out ``CCFLAGS = $(DEBUG_FLAGS)`` and ``CCFLAGS = $(PROFILE_FLAGS)``.

#. Download and install Intel Parallel Studio.

	* The Intel compiler is a commercial product, but you may be able to use it freely on a trial basis.

#. Open the special Intel compiler command prompt for the appropriate processor type.  You can find this in the start menu.

	* You cannot use the PowerShell or the usual command prompt.

#. :samp:`cd` :samp:`{twroot}`:samp:`\\core\\source`
#. :samp:`nmake /F win.make`
#. The makefile should automatically copy the executable into your :samp:`~\\bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~\\Run`, but these will not be used.
