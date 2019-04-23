Core Installation
=================

Compiler Notes
--------------

We support several compilers for desktop turboWAVE.  Explicit instructions are given for the following OS/compiler combinations:

.. csv-table:: Table I. Desktop C++ Compilers.
	:header: "Umbrella", "Compiler Command", "Operating Systems"

	"GNU", ``g++``, "Linux"
	"LLVM", ``clang++``, "Linux, MacOS"
	"Intel Parallel Studio", ``icl``, "Windows"
	"Microsoft Visual Studio", ``cl``, "Windows"

The GNU and LLVM compilers can be freely downloaded.  Intel and Microsoft compilers
are commercial products, but typically offer free trial downloads.  Other combinations
are possible, these are merely the ones that are tested and directly supported.

On HPC systems, we expect a suitable compiler to be pre-installed by the system
administrators, although there may be modules to load, unload, or swap.

Note on Packages
----------------

Linux distributions typically have a native package management system,
which we use to install the compiler and other components we need.
The package management system grabs software from internet repositories,
and importantly, checks whether the software depends on or conflicts with some other
software.  If there is a dependency, the system is supposed to grab the whole hierarchy
of software that you need in order to complete the requested installation.  MacOS and
Windows have package managers that can be installed as add-ons.

.. csv-table:: Table II. Package Managers.
	:header: "System", "Command Interface", "Operating Systems"

	"Anaconda", ``conda``, "Any"
	"Chocolatey", ``choco``, "Windows"
	"Debian", ``apt``, "Ubuntu (native)"
	"Homebrew", ``brew``, "MacOS"
	"MacPorts", ``port``, "MacOS"
	"RPM", ``yum``, "CentOS/RHEL/SL (native)"

OpenMP Option
-------------

It is recommended that turboWAVE always be compiled with OpenMP, even if shared memory parallelism is not going to be used.  There are two reasons for this:

	1. All turboWAVE SIMD support relies on OpenMP.  So even if shared memory is not used at all, there can be a significant performance penalty when OpenMP is missing.
	2. OpenMP enabled executables are tested more frequently than executables without OpenMP.

Core Install by OS
------------------

.. toctree::
	:maxdepth: 1

	core-mac
	core-RHEL
	core-ubuntu
	core-windows
	core-hpc
