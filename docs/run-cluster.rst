Local Cluster Runs
==================

In referring to a local cluster, we have in mind a modest scale system maintained by the user, without any job queuing system.  This document assumes that OpenMPI is used as the message passing library.  For other MPI implementations the procedure should be similar.

Running an Example
------------------

The process of running an example on a local cluster is similar to the desktop case.

#. Select an example, denoted :samp:`{example}.tw`, and copy to the login node's ``Run`` directory.
#. Edit the ``decomposition`` directive in the ``new grid`` block.

	* Make sure the product of the three integers corresponds to the number of cores you wish to test.
	* The number of cells per node along each axis must be even.

#. Copy the edited :samp:`{example}.tw` to the ``Run`` directory on every other node.
#. Let :samp:`{N}` be the number of cores you wish to test.
#. Edit the hosts file such that the number of nodes requested matches the number of cores to be tested.
#. Execute :samp:`mpirun -np {N} --hostfile {file} tw3d -c 1 --input-file {example}.tw`

	* The argument :samp:`{file}` is the hosts file.
	* At present the interactive thread is disabled for cluster runs.  To check the progress navigate to the ``Run`` directory in another terminal and execute ``cat twstat``.

#. When the run completes, the data will be on the node hosting MPI rank zero.

The above procedure uses MPI only, because the argument ``-c 1`` selects a single OpenMP thread per MPI process.  Using more than one OpenMP thread might be advantageous in some cases.

.. _args_cluster:

Command line arguments
----------------------

For cluster installations the command line reads

.. py:function:: <launcher> -np <procs> [--hostfile <nodes>] tw3d [optional arguments...]

The options for the launcher (which can depend on the MPI library) are typically

.. program:: launcher

.. option:: -np <procs>

	number of MPI processes to launch

.. option:: --hostfile <nodes>

	name of the file listing the hosts to be used in case of a multi-node run

.. option:: tw3d

	the turboWAVE executable file

Optional arguments following ``tw3d`` are as follows:

.. program:: tw3d

.. option:: -c <threads>

	number of OpenMP threads (see below for default)

.. option:: --input-file <file>, -i <file>

	name or path of the file to use as the input file (default=stdin)

.. option:: --platform <search_string>

	select an OpenCL platform with the search string in its name

.. option:: --device <search_string>

	select an OpenCL device with the search string in its name.  This can also be a comma-delimited list of device numbers.

.. option:: --output-level <level>

	an integer determining how much information to log.  If 0 only rank 0 logs, if 1 every rank logs.

.. option:: --restart

	if present, causes initial data to be loaded from a checkpoint.

.. option:: --unit-test <test>

	if present, bypasses simulation and runs units test instead.  The argument can be a single test, or ``--all``

.. option:: --no-interactive

	if present, suppresses the interactive thread.

.. option:: --version, -v

	if present, prints the version number.  If this is the only argument, no simulation is attempted.

.. option:: --help, -h

	if present, prints the command line arguments and the link to the online documentation.  If this is the only argument, no simulation is attempted.

Error Handling
--------------

It is important to pay attention to the output file if you are having problems.  If the code stops without reporting an error in the terminal window, you may still be able to get some feedback.  The procedure is as follows.

	#. Run the simulation again with the command line option ``--output-level 1``
	#. If the error is not reported on the console, try :samp:`grep ERROR *stdout*`
