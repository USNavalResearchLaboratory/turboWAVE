
.. toctree::
	:maxdepth: 2

Getting Components
==================

General Information
-------------------

The `Git <https://git-scm.com/>`_ version control system is fundamental to the evolution, storage, and retrieval of turboWAVE components.

TurboWAVE components are maintained in a single repository.  Unless you are doing development yourself, you will typically interact only with the `master branch <https://git-scm.com/docs/git-branch/>`_. The `commits <https://git-scm.com/docs/git-commit/>`_ that make up the master branch have varying degrees of stability.  Stable releases are identified by a `tag <https://git-scm.com/docs/git-tag/>`_.

The turboWAVE repository is at `<https://github.com/USNavalResearchLaboratory/turboWAVE>`_.

Cloning the Repository
----------------------

The first step to getting the components is to clone the repository.  This copies the entire history of the code to your local computer (storage requirement is in the tens of megabytes).  Perform the following procedure:

	#. Open a terminal window
	#. Test to see if you have Git installed by executing :samp:`git --version`
	#. Install Git if necessary.

		* Anaconda --- :samp:`conda install git`
		* CentOS/RHEL/SL --- :samp:`sudo yum install git`
		* Homebrew --- :samp:`brew install git`
		* MacPorts --- :samp:`sudo port install git`
		* Ubuntu --- :samp:`sudo apt install git`

	#. Navigate to the directory where you want to install turboWAVE (you don't need to make an enclosing directory).
	#. :samp:`git clone https://github.com/USNavalResearchLaboratory/turboWAVE.git`
	#. If you like you can give the turboWAVE root directory another name, we will call it :samp:`{twroot}` from now on.

Switching to a Stable Version
-----------------------------

When you clone the repository the active files (the version you have checked out) will likely be the latest commit, which is not necessarily the most stable.  In order to select a stable version perform the following procedure.

	#. Open a terminal and navigate to :samp:`{twroot}`.
	#. :samp:`git tag --list`
	#. Choose the latest tag, :samp:`{latest_tag}`.
	#. :samp:`git checkout {latest_tag}`.
