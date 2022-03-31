HPC Runs
==============

In referring to High Performance Computing (HPC), we have in mind a large scale computing cluster, managed by a professional staff, and accessed remotely.

When running on an HPC system, you will be dealing with a queuing system such as PBS or SLURM.  The turboWAVE executable will usually not be run directly.  Instead, a launcher program such as :samp:`aprun` or :samp:`srun` will take :samp:`tw3d` as an argument.  The launcher program is usually invoked from a job script.

It will likely be necessary to study the documentation for the particular HPC system in question.


Running an Example
-------------------

#. Pick some example from :samp:`{turboWAVE}/core/examples`.  For this test a 3D example is appropriate.
#. For definiteness, let us use :samp:`{turboWAVE}/core/examples/pgc/beatwave-3d.tw`
#. Let us assume that you will use :samp:`scp` to copy files to the HPC system
#. :samp:`scp {turboWAVE}/core/examples/pgc/beatwave-3d.tw {user@HPC_URL:HPC_scratch_directory}/stdin`
#. This puts the input file in the HPC scratch space with the name :samp:`stdin`.  By default, turboWAVE expects the input file to have the name :samp:`stdin`, and to be in the working directory.
#. Create a batch script according to the instructions for your machine.  Note the decomposition in the example input file.  The product of the three integers is the number of MPI processes you should request.  Assuming you don't fork threads, the number of cores you request should be the same as the number of MPI processes.  If you fork threads, then you must multiply the total cores by the number of threads you fork from each MPI process.  When forming your batch script, be careful to note whether a parameter refers to the cores per node, or the total cores for the entire job.
#. Submit the script according to the instructions for your machine.
#. As the problem runs, you can monitor the progress by examining the contents of the file :samp:`twstat`.
#. When the run is finished, you should have several files with the extension :samp:`npy`.  This is a simple binary format often used with numerical Python.

The visualization workflow is identical to that described in :doc:`run-desktop`.  Remote visualization using the command line tools should work provided X-forwarding is allowed. If you are allowed to run a Jupyter notebook server on the remote machine, the interactive DataViewer should work as well.

Command line arguments
----------------------

For HPC the command line specification is technically the same as is it is for local clusters. However, the following points should be noted:

	* The command is issued from within a job script, with launcher and associated command line options determined by the particular queuing system and operating environment.
	* The hostfile is rarely used, instead the queuing system designates the nodes.
	* The number of OpenMP threads forked by each MPI process is typically determined by setting an environment variable in the job script.  On the other hand, if the ``-c <threads>`` argument is used, it takes precedence.

Error Handling
---------------

It is important to pay attention to the output file if you are having problems.  If the code stops without reporting an error in the main output file (usually named by you in the batch script), you may still be able to get some feedback.  The procedure is as follows.

	#. Run the simulation again with the command line option ``--output-level 1``
	#. If the error is not reported on the console, try :samp:`grep ERROR *stdout*`
