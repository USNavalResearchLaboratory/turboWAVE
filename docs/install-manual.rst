
.. toctree::
	:maxdepth: 2

Manual Installation
///////////////////

It is recommended to use the graphical installer program.  However, if you want to take full control of the configuration this page may be useful.

Getting Components
==================

General Information
-------------------

The `Git <https://git-scm.com/>`_ version control system is fundamental to the evolution, storage, and retrieval of turboWAVE components.

TurboWAVE components are maintained in a single repository.  Unless you are doing development yourself, you will typically interact only with the `master branch <https://git-scm.com/docs/git-branch/>`_. The `commits <https://git-scm.com/docs/git-commit/>`_ that make up the master branch have varying degrees of stability.  Stable releases are identified by a `tag <https://git-scm.com/docs/git-tag/>`_.  The tag may have letter suffixes ``a``, ``b``, or ``rc``, indicating alpha, beta, or release candidate.  A numerical tag with no letter suffix is the most stable.

The turboWAVE repository is at `<https://github.com/USNavalResearchLaboratory/turboWAVE>`_.

Cloning the Repository
----------------------

The first step to getting the components is to clone the repository.  This copies the entire history of the code to your local computer (storage requirement is in the tens of megabytes).  Perform the following procedure:

	#. Open a terminal window
	#. Test to see if you have Git installed by executing :samp:`git --version`
	#. Install Git if necessary.

		* Anaconda --- :samp:`conda install git`
		* CentOS/RHEL/SL --- :samp:`sudo yum install git`
		* MacOS with Homebrew --- :samp:`brew install git`
		* MacOS with MacPorts --- :samp:`sudo port install git`
		* Ubuntu --- :samp:`sudo apt install git`
		* Windows PowerShell with Chocolatey --- :samp:`choco install git` (run as administrator)

	#. Navigate to the directory where you want to install turboWAVE (you don't need to make an enclosing directory).
	#. :samp:`git clone https://github.com/USNavalResearchLaboratory/turboWAVE.git`

Switching to a Stable Version
-----------------------------

When you clone the repository the active files (the version you have checked out) will likely be the latest commit, which is not necessarily the most stable.  In order to select a stable version perform the following procedure.

	#. Open a terminal and navigate to :samp:`{turboWAVE}`.
	#. :samp:`git tag --list`
	#. Choose the latest tag without a letter suffix, :samp:`{latest_stable_tag}`.
	#. :samp:`git checkout {latest_stable_tag}`.

		* You may be in a detached state.  If you want to restore the state later you can run ``git checkout master``, or, to automatically throw out any build products or other changes, add the ``-f`` flag.

.. Note::

	if you are viewing this documentation on ``readthedocs``, you are likely viewing the ``latest`` documentation which is not necessarily the ``stable`` documentation.  There should be a hovering dropdown to select between the two.

Tools Installation
==================

General Notes
-------------

Most of the tools are built on Python 3.  It is highly recommended to operate in a virtual environment.  Virtual environments can be created using native Python tools, or Anaconda.

Tools are synchronized with core to within major and minor version numbers (the patch and build numbers may differ).  For the best results, keep the installed versions consistent.

Python Module twutils
---------------------

#. You can create a conda environment with the necessary modules and scripts using :samp:`conda create -n {NAME} -c dfxgordon twutils`, where :samp:`{NAME}` is the name of the environment (your choice).
#. You can also install from PyPi by typing ``pip install twutils``.  Note that in this case git must be installed separately.
#. You can install from the local repository using ``pip install .`` (note dot) from within the :samp:`{turboWAVE}/tools/twutils` directory.

Python DataViewer
-----------------

#. The Python DataViewer has to be run in a Jupyter Notebook.
#. Copy :samp:`{turboWAVE}/tools/DataViewer.ipynb` to :samp:`~/bin`
#. Create a directory :samp:`~/.jupyter/custom/` and copy :samp:`{turboWAVE}/tools/config-files/custom.css` to the new directory.

Input File Language Support
----------------------------

You can add language support for turboWAVE input files in various editors.

#. For ``VS Code`` and ``Atom``, language support can be installed using the editor's own GUI, search for ``turbowave``.

#. To enable turboWAVE input file syntax highlights with the :samp:`micro` editor

	* Copy :samp:`{turboWAVE}/tools/config-files/turbowave.micro.yaml` to ``%HomePath%\.config\micro\syntax\`` (Windows) or ``~/.config/micro/syntax/`` (others).

#. To enable turboWAVE input file syntax highlights with the :samp:`nano` editor

	* Copy :samp:`{turboWAVE}/tools/config-files/.turbowave.nano` to your home directory.
	* Edit (create if necessary) ``%HomePath%\nano.rc`` (windows) or ``~/.nanorc`` (others) and add the line ``include "~/.turbowave.nanorc"``.

#. To enable turboWAVE input file syntax highlights with the :samp:`vim` editor

	* Copy :samp:`{turboWAVE}/tools/config-files/filetype.vim` to ``%HomePath%\vimfiles\`` (Windows) or ``~/.vim/`` (others)
	* Copy :samp:`{turboWAVE}/tools/config-files/turbowave.vim` to ``%HomePath%\vimfiles\syntax\`` or ``~/.vim/syntax/`` (others).

#. To enable turboWAVE input file syntax highlights with the :samp:`neovim` editor

	* Following is only tested on Linux as of this writing
	* Install the ``nvim-treesitter`` package.  You can use your favorite ``neovim`` package manager, or clone ``nvim-treesitter`` into ``~/.local/share/nvim/site/pack/ts/start`` (Linux).
	* Copy :samp:`{turbowave}/tools/config-files/init.vim` to ``~/.config/nvim`` (edit the settings if you like).
	* Copy ``highlights.scm`` from ``tree-sitter-turbowave`` to ``~/.local/share/nvim/site/queries/turbowave``.
	* Create a file ``~/.config/nvim/ftdetect/turbowave.vim`` with contents ``au BufRead,BufNewFile *.tw set filetype=turbowave``.
	* Start a new ``nvim`` session and run the command ``:TSInstall turbowave``.
	* You can run ``:checkhealth nvim-treesitter`` to see if it worked and diagnose any problems.  The list of parsers should include ``turbowave``.

Core Installation
=================

Basic Approach
---------------

The basic installation workflow is build from source, and copy files into user space directories.  This workflow is managed using GNU ``make``.  The ``makefile`` is configured to build and copy files in one step.  This is distinct from the usual practice of separating ``make``, which builds the binaries, from ``make install`` which copies them to standard locations.  The philosophy behind this is that some runtime workflows involve frequently tweaking source code, and in this setting, combining ``make`` and ``make install`` is less error prone.  The more typical approach can be applied by using ``make tw3d_release`` followed by ``make install``.

Configuration
-------------

For manual installations, configuration is done by directly editing ``core/build/makefile``.  Editing the ``makefile`` is usually easy.  In most cases, you only have to adjust the constants in the input variables block, which is prominently identified by comments.

If the ``makefile`` senses an activated conda environment, it copies the executable binary into that environment.  Otherwise, it tries to copy the executable to an *existing* directory ``~/bin`` (``~\Scripts`` on Windows).  OpenCL kernel files are copied to an *existing* directory ``~/Run``.  You can change these directories by editing ``BINARY_PATH`` and ``WORK_PATH`` constants.

Compiler Notes
--------------

For desktop turboWAVE, it is recommended to use either GNU ``g++`` or LLVM ``clang++`` compilers. The Intel compiler may work, but is not being tested at present.  Support for ``nmake`` is dropped, as of this writing.

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

For desktop systems, the installation instructions assume the use of turboWAVE's internal MPI implementation.  This is primarily for convenience.  Better performance might be obtained by using an external MPI implementation, even on the desktop.  This can be particularly advantageous if measures for controlling thread or CPU affinity are undertaken.

Performance Tuning Parameters
-----------------------------

There are a few parameters hard coded in source that can be used to tune performance.  These are as follows.

#. Vector length - In the makefile, you should adjust the constant ``VBITS`` to match the vector processor of the target system.  For example, systems with AVX2 should set this to 256, while systems with AVX-512 should set this to 512.

#. Particle bundle size - In ``definitions.h``, you should change ``max_bundle_size`` to be a multiple of ``VBITS`` divided by 32. The optimal choice is difficult to predict, but a rule of thumb is to make it close to the typical number of particles per cell.

#. There are some constants in ``Pusher.cpp``, in the method ``Species::Push``, that it may be advantageous to adjust.  Namely, the constants ``min_particles_per_task`` and ``preferred_tasks``.  This has to do with how particles are partitioned among OpenMP threads.

Core Install for Local Clusters
===============================

In referring to a local cluster, we have in mind a modest scale system maintained by the user, without any job queuing system.  This document assumes that OpenMPI is used as the message passing library.  For other MPI implementations the procedure should be similar.  The makefile has a provision for Intel MPI as well.

This type of installation differs from a desktop installation in at least two ways:

	* TurboWAVE must be linked against an external MPI library

	* Provision must be made for file access across distributed nodes

At present almost all systems use little-endian binary numbers.  If your system is big-endian, you must change a flag in the makefile.

Single Node System
------------------

Even on a single shared-memory node, an external MPI library can give performance gains.  To perform the installation in such a case, carry out the following steps.

	#. First follow all the steps for the desktop installation that most closely resembles the operating system and compiler configuration of the local cluster.  After successfully compiling in this mode, issue ``make clean`` and continue as follows.

	#. Install the external MPI library.  For OpenMPI, most package managers should work, e.g., execute ``sudo apt install libopenmpi-dev`` for Debian, etc..

	#. In the makefile, uncomment :samp:`PLATFORM = OPENMPI`, and comment out the other platforms.

	#. Type :samp:`make`

	#. If there are errors, you may need to edit :samp:`makefile`.  Things to watch out for include:

		* There is typically a special compile command which creates an environment that allows the compiler to find MPI headers and libraries.  You may need to find this command for your particular MPI implementation, and set the variables :samp:`TW_Compiler` and :samp:`TW_Linker` to this command.

		* You may have to set environment variables that help the MPI enabled compiler invocation find the underlying compiler you are trying to use.

Multiple Node System
--------------------

#. Install the external MPI library on every node.

#. Choose one node as the "login node".  Set up RSA key pairs such that from the login node, you can log-in to any other node without a password.

#. Compile the executable such that it is compatible with the hardware and operating system on every node (for homogeneous clusters there is no issue).

#. Create a ``Run`` directory on each node, and copy the executable into each ``Run`` directory.

#. Create a "hosts file" listing the nodes.  This should go in the ``Run`` directory from which parallel jobs will be launched (the login node). See the MPI library's documentation for more.

#. If the MPI library is installed in user space, it may be necessary to set environment variables such as ``PATH`` and ``LD_LIBRARY_PATH``, on every node. See the MPI library's documentation for more.

Core Install for HPC
====================

In referring to High Performance Computing (HPC), we have in mind a large scale computing cluster, managed by a professional staff, and accessed remotely.

We will focus on Cray systems in this instruction.  For other vendors the below will often go through with little modification.  At present almost all systems use little-endian binary numbers.  If your system is big-endian, you must change a flag in the makefile.

HPC modules
-----------

HPC systems usually have a system that allows for easy loading, unloading, and swapping of modules.  The commands include :samp:`module load {module}`, :samp:`module unload {module}`, and :samp:`module swap {oldModule} {newModule}`.
For Cray systems, we usually want the Intel compiler.  This would be loaded with :samp:`module load PrgEnv-intel`.  If another compiler is already loaded you can use, say, :samp:`module swap PrgEnv-cray PrgEnv-intel`.
The NERSC supercomputing center automatically loads the darshan module, which doesn't play with turboWAVE.  Before compiling run :samp:`module unload darshan`.

Compiling on Cray Systems
-------------------------

  #. Make a directory on the HPC system for turboWAVE source.  We denote it :samp:`{turbowave}`

  #. Copy everything in the :samp:`{turboWAVE}/core/source/` directory to :samp:`{turbowave}`.

  #. Navigate to :samp:`{turbowave}`

  #. Edit :samp:`{turbowave}/makefile`

  #. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = CRAY`, :samp:`HARDWARE_ACCEL = OMP`, and :samp:`COMPILER_PREF = INTEL`, should be uncommented.

  #. Type :samp:`make`

  #. You must manually copy the executable to the scratch directory.  For example, at NERSC, this would be done with :samp:`cp tw3d $SCRATCH`.
