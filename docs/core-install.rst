Core Installation
=================

Core Install by OS
------------------

If you want to dive right in, select your operating system or configuration below.  If you want some background, skip to the below subsections first.

.. toctree::
	:maxdepth: 1

	core-mac
	core-RHEL
	core-ubuntu
	core-windows
	core-cluster
	core-hpc


Basic Approach
---------------

The basic installation workflow is build from source, and copy files into user space directories.  This workflow is managed using GNU ``make``.  The ``makefile`` is configured to build and copy files in one step.  This is distinct from the usual practice of separating ``make``, which builds the binaries, from ``make install`` which copies them to standard locations.  The philosophy behind this is that some runtime workflows involve frequently tweaking source code, and in this setting, combining ``make`` and ``make install`` is less error prone.  The more typical approach can be applied by using ``make tw3d_release`` followed by ``make install``.

Compiler Notes
--------------

We support several compilers for desktop turboWAVE.  Explicit instructions are given for the following OS/compiler combinations:

.. csv-table:: Table I. Desktop C++ Compilers.
	:header: "Umbrella", "Compiler Command", "Operating Systems"

	"GNU", ``g++``, "Linux, MacOS"
	"LLVM", ``clang++``, "Linux, Windows"
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

External MPI Option
-------------------

For desktop systems, the installation instructions assume the use of turboWAVE's internal MPI implementation.  This is primarily for convenience.  Better performance might be obtained by using an external MPI implementation, even on the desktop.  This can be particularly advantageous if measures for controlling thread or CPU affinity are undertaken.  For guidance on this, see :doc:`core-cluster`.

Performance Tuning Parameters
-----------------------------

There are a few parameters hard coded in source that can be used to tune performance.  These are as follows.

#. Vector length - In ``definitions.h``, you should adjust the constant ``vec_align_bytes`` to match the vector processor of the target system.  For example, systems with AVX2 should set this to 32, while systems with AVX-512 should set this to 64.

#. Particle bundle size - In ``definitions.h``, whenever ``vec_align_bytes`` is changed, you should also change ``max_bundle_size`` to be a multiple of ``vec_align_bytes`` divided by four. The optimal choice is difficult to predict, but a rule of thumb is to make it close to the typical number of particles per cell.

#. There are some constants in ``Pusher.cpp``, in the method ``Species::Push``, that it may be advantageous to adjust.  Namely, the constants ``min_particles_per_task`` and ``preferred_tasks``.  This has to do with how particles are partitioned among OpenMP threads.
