Core Installation
=================

General Notes
-------------

We consider three compilers for desktop turboWAVE.  The first is the GNU C++ compiler,
typically referred to as GCC, and invoked as ``g++``.  The second is the LLVM C++
compiler, typically referred to as CLANG, and invoked as ``clang++``.
The third is the compiler Microsoft packages with Visual Studio.

UNIX operating systems typically have some kind of package management system,
which we will use to install the compiler and other components we need.
The package management system grabs software from internet repositories,
and importantly, checks whether the software depends on or conflicts with some other
software.  If there is a dependency, the system is supposed to grab the whole hierarchy
of software that you need in order to complete the requested installation.

On MacOS, the package management systems we support are MacPorts and Homebrew.
On Ubuntu, the system is Debian.  On RHEL/SL/CentOS, the system is RPM.
On Windows we do not use any package manager.

On HPC systems, we expect a suitable compiler to be pre-installed by the system
administrators, although there may be modules to load, unload, or swap.

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
