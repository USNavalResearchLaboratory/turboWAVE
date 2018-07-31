HPC Runs
==============

In referring to High Performance Computing (HPC), we have in mind a large scale computing cluster, managed by a professional staff, and accessed remotely.

When running on an HPC system, you will be dealing with a queuing system such as PBS or SLURM.  The turboWAVE executable will usually not be run directly.  Instead, a launcher program such as :samp:`aprun` or :samp:`srun` will take :samp:`tw3d` as an argument.  The launcher program is usually invoked from a batch script.

It will likely be necessary to study the documentation for the particular HPC system in question.

Parallel Programming Theory
----------------------------

In order to run turboWAVE, it is useful to understand a little about the parallelization methods.  TurboWAVE uses a combination of distributed memory and shared memory methods.  The distributed memory method is called "domain decomposition", and corresponds to physically paritioning the simulation region into chunks that can be worked on by different processors.  The memory that holds each chunk is distributed in the sense that it may be connected only by a network cable.  The shared memory method is called "fork-join".  In this approach independent software threads are running on different cores, but in a setting where each core has equal access to the memory (e.g., the threads could be on different cores of the same multi-core processor).

The distributed memory model is implemented using software called MPI.  The shared memory model is implemented using software called OpenMP.

As the user, you have to choose how to partition the problem.  In particular, you must choose the number of MPI processes, and the number of threads.  If you opt to have more than one MPI process, you may also want to control exactly how the physical simulation domain is partitioned into chunks.

On HPC systems in particular, the term "node" refers to a set of physical processor cores that have access to the same shared memory.  Often you must specify resources requested in terms of nodes.

Running an Example
-------------------

#. Pick some example from :samp:`{twroot}/core/examples`.  For this test a 3D example is appropriate.
#. For definiteness, let us use :samp:`{twroot}/core/examples/pgc/beatwave-3d.txt`
#. Let us assume that you will use :samp:`scp` to copy files to the HPC system
#. :samp:`scp {twroot}/core/examples/pgc/beatwave-3d.txt {user@HPC_URL:HPC_scratch_directory}/stdin`
#. This puts the input file in the HPC scratch space with the name :samp:`stdin`.  TurboWAVE always expects the input file to have the name :samp:`stdin` or :samp:`stdin.txt`, and to be in the working directory.
#. Create a batch script according to the instructions for your machine.  Note the decomposition in the example input file.  The product of the three integers is the number of MPI processes you should request.  Assuming you don't fork threads, the number of cores you request should be the same as the number of MPI processes.  If you fork threads, then you must multiply the total cores by the number of threads you fork from each MPI process.  When forming your batch script, be careful to note whether a parameter refers to the cores per node, or the total cores for the entire job.
#. Submit the script according to the instructions for your machine.
#. As the problem runs, you can monitor the progress by examining the contents of the file :samp:`twstat`.
#. When the run is finished, you should have several files with the extension :samp:`dvdat`.  This is a simple binary format.  The twutils Python package has a function to read data into numpy arrays from this type of file.  If you want to see an example of how to read this file from C++, you can look in :samp:`{twroot}/tools/twpost`.
#. Bring one or more data files over to your desktop computer (e.g., using :samp:`scp`).  For this example, retrieve :samp:`phi.dvdat`.
#. Let us plot the results using DataViewer.  If you have the native MacOS version, double-click on :samp:`phi.dvdat` and advance the "Frame" slider.  You may like to go to the "View" menu and select "Autoscale Plot" to get a better color contrast.
#. If you do not have a native DataViewer, you can run the python version.  Open a terminal window and navigate to :samp:`~/bin`, or wherever :samp:`DataViewer.ipynb` is.
#. Activate your virtual environment (see :doc:`tools-install`)
#. :samp:`jupyter notebook`
#. Click on :samp:`DataViewer.ipynb`
#. Locate the path variable in the source, and change to the directory where you downloaded the data.
#. Click on the button to run the notebook
#. Use the File dropdown to select :samp:`phi.dvdat`.
#. Advance the Frame slider to examine the frames


Error Handling
---------------

It is important to pay attention to the output file if you are having problems.  If the code stops without reporting an error in the main output file (usually named by you in the batch script), you may still be able to get some feedback.  The procedure is as follows.

	#. In the input file, add the line :samp:`stdout = full`
	#. This line can go anywhere except within a :samp:`new` block or :samp:`generate` block
	#. Run the problem again
	#. If the error is not reported in the main output, try :samp:`grep ERROR *stdout*`
