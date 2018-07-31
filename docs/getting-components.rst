
.. toctree::
	:maxdepth: 2

Getting Components
==================

General Information
-------------------

The `Git <https://git-scm.com/>`_ version control system is fundamental to the evolution, storage, and retrieval of turboWAVE components.

TurboWAVE components are maintained in two separate Git repositories, :samp:`core` and :samp:`tools`.  Unless you are doing development yourself, you will typically interact only with the `master branch <https://git-scm.com/docs/git-branch/>`_. The `commits <https://git-scm.com/docs/git-commit/>`_ that make up the master branch have varying degrees of stability.  Stable releases are identified by a `tag <https://git-scm.com/docs/git-tag/>`_.  The form of the tag is :samp:`TW.{YYYY}.{MM}.{DD}`.  The definition of a stable release is one for which all the example input files have been verified.

You will obtain the turboWAVE components by visiting a web page that contains up-to-date copies of the Git repositories.  At the time of this writing the web page is inward facing at NRL only (cannot be accessed from outside without using VPN).  The address is :samp:`predator.nrl.navy.mil`.

Cloning the Repository
----------------------

The first step to getting the components is to clone the repository.  This copies the entire history of the code to your local computer (storage requirement is in the tens of megabytes).  To clone the repository you will have to register on :samp:`predator.nrl.navy.mil` and follow the instructions for generating an SSH key (go to user settings and look in the left column).  You need an SSH key on each local computer that will access the server.  Once this is done perform the following procedure:

	#. Open a terminal window
	#. Test to see if you have Git installed by executing :samp:`git --version`
	#. Install Git if necessary.

		* Anaconda --- :samp:`conda install git`
		* CentOS/RHEL/SL --- :samp:`sudo yum install git`
		* Homebrew --- :samp:`brew install git`
		* MacPorts --- :samp:`sudo port install git`
		* Ubuntu --- :samp:`sudo apt install git`

	#. Create a directory :samp:`{twroot}` and navigate there.
	#. :samp:`git clone git@predator.nrl.navy.mil:turbowave/tools.git`
	#. :samp:`git clone git@predator.nrl.navy.mil:turbowave/core.git`

Switching to a Stable Version
-----------------------------

When you clone the repository the active files (the version you have checked out) will likely be the latest commit, which is not necessarily the most stable.  In order to select a stable version perform the following procedure.

	#. Open a terminal and navigate to :samp:`{twroot}/core`.
	#. :samp:`git tag --list`
	#. Choose the tag with the most recent date, :samp:`{latest_tag}`.
	#. :samp:`git checkout {latest_tag}`.
