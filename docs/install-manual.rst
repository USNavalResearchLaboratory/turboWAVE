
.. toctree::
	:maxdepth: 2

Manual Installation
///////////////////

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

		* You may be in a detached state.  If you want to restore the state later you can run ``git checkout master``, or, to automatically throw out any changes, add the ``-f`` flag.

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
#. You can install from the local repository using ``pip install .`` (note dot) from within the :samp:`{turboWAVE}/tools/twutils` directory.

Python DataViewer (optional)
----------------------------

#. The Python DataViewer has to be run in a Jupyter Notebook.
#. Copy :samp:`{turboWAVE}/tools/DataViewer.ipynb` to wherever you keep your notebooks
#. Create a directory :samp:`~/.jupyter/custom/` and copy :samp:`{turboWAVE}/tools/config-files/custom.css` to the new directory.

Input File Support (optional)
-----------------------------

You can add language support for turboWAVE input files in various editors.

#. For ``VS Code``, language support can be installed using the editor's own GUI, search for ``turbowave``.

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

Configuration
-------------

#. Navigate to ``core/source`` and type ``meson setup build``

	* For MacOS you probably need to set ``CXX`` to the full path of the Homebrew or MacPorts compiler.

#. Set options

	* Options are set using :samp:`meson configure -D{key}={value} build`
	* Set vector width, e.g., if you have AVX512, set ``vbits=512``
	* To use turboWAVE's internal MPI set ``hpc=false``
	* To use external MPI, set ``hpc=true`` (suitable library must be installed)
	* To see all available options type ``meson configure build``

Performance Tuning Parameters
-----------------------------

There are a few parameters hard coded in source that can be used to tune performance.  These are as follows.

#. Particle bundle size - In ``definitions.h``, you should change ``max_bundle_size`` to be a multiple of ``vbits`` divided by 32. The optimal choice is difficult to predict, but a rule of thumb is to make it close to the typical number of particles per cell.

#. There are some constants in ``Pusher.cpp``, in the method ``Species::Push``, that it may be advantageous to adjust.  Namely, the constants ``min_particles_per_task`` and ``preferred_tasks``.  This has to do with how particles are partitioned among OpenMP threads.

Build and Install
-----------------

#. From the ``source`` directory type ``meson compile -C build``, or from the ``build`` directory type ``meson compile``

#. The executable can be installed to a per-OS standard location with ``meson install``.  You can use ``meson configure`` to change the default install location.  You can also simply move the executable.

Core Install for HPC
====================

In referring to High Performance Computing (HPC), we have in mind a large scale computing cluster, managed by a professional staff, and accessed remotely.

HPC modules
-----------

HPC systems usually have a system that allows for easy loading, unloading, and swapping of modules.  The commands include :samp:`module load {module}`, :samp:`module unload {module}`, and :samp:`module swap {oldModule} {newModule}`.

Copy Components
---------------

You can usually clone the repository and setup the ``twutils`` environment the same way as on a desktop.  What is required is at least the ``source`` directory and Meson.

Build
-----

#. :samp:`conda activate {NAME}`

#. Set ``CXX`` and ``CC`` environment variables to the compiler wrappers recommended by the HPC center, if applicable.  It is important that this be done first.

	* For example, you might have ``export CXX=$(which CC)`` and ``export CC=$(which cc)``
	* If you will be building often add these to the login environment

#. Edit ``meson.build`` and comment out ``deps += [dependency('mpi')]``

	* if the build fails you can try uncommenting this line

#. Navigate to ``core/source``

#. ``meson setup build``

#. ``cd build``

#. ``meson configure -Dhpc=true`` (and any other options that may be appropriate)

#. ``meson compile``

#. Copy the executable to the runtime directory
