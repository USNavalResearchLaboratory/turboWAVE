Core Install for Local Clusters
===============================

In referring to a local cluster, we have in mind a modest scale system maintained by the user, without any job queuing system.  This document assumes that OpenMPI is used as the message passing library.  For other MPI implementations the procedure should be similar.  The makefile has a provision for Intel MPI as well.

This type of installation differs from a desktop installation in at least two ways:

	* TurboWAVE must be linked against an external MPI library

	* Provision must be made for file access across distributed nodes

At present almost all systems use little-endian binary numbers.  If your system is big-endian, you must change the value of a boolean variable in :samp:`definitions.h`.

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
