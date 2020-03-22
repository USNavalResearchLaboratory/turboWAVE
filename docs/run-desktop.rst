Desktop Runs
============

.. caution::

	This document assumes you have followed the installation instructions precisely.

.. note::

	Once installed, running turboWAVE should be very nearly the same for any desktop system. For Windows, this consistency depends on using the PowerShell as the terminal window.

Running an Example
------------------

#. Pick some example from :samp:`{twroot}/core/examples`.  For this test you should avoid the 3D examples.
#. For definiteness, let us use :samp:`{twroot}/core/examples/pgc/LWFA-coulomb.tw`
#. Open a terminal window and navigate to :samp:`~/Run`
#. :samp:`cp {twroot}/core/examples/pgc/LWFA-coulomb.tw stdin`
#. This puts the input file in :samp:`~/Run` with the name :samp:`stdin`.  By default, turboWAVE assumes the input file is in the working directory with the name :samp:`stdin`.
#. :samp:`tw3d -n 4`
#. The above command runs the problem with 4 MPI processes and 1 thread per process.  Of course this choice may not be optimal for your system, method of compiling, etc., but it should suffice for this example.
#. As the problem runs, you can press the enter key to prompt turboWAVE to report the current step.  Enter :samp:`help` to get the full list of interactive commands.
#. When the run is finished, you should have several files with the extension :samp:`dvdat`.  This is a simple binary format.  The twutils Python package has a function to read data into numpy arrays from this type of file.  If you want to see an example of how to read this file from C++, you can look in :samp:`{twroot}/tools/twpost`.
#. Let us plot the results using DataViewer.  If you have the native MacOS or Windows version, double-click on :samp:`phi.dvdat` and advance the "Frame" slider.  You may like to go to the "View" menu and select "Autoscale Plot" to get a better color contrast.
#. If you do not have a native DataViewer, you can run the python version.  Open a terminal window and navigate to :samp:`~/bin`, or wherever :samp:`DataViewer.ipynb` is.
#. Activate your virtual environment (see :doc:`tools-install`)
#. :samp:`jupyter notebook`
#. Click on :samp:`DataViewer.ipynb`
#. Locate the path variable in the source, and change to your own Run directory. Prefixing the string with ``u`` allows forward slashes to be used as directory separators irrespective of operating system.
#. Click on the button to run the notebook
#. Use the File dropdown to select :samp:`phi.dvdat`.
#. Advance the Frame slider to the last frame
#. Your window should look something like Fig. 1.

.. figure:: LWFA-coulomb.png
	:figwidth: 80%

	Fig. 1 --- Python DataViewer output of the scalar potential produced by the :file:`LWFA-coulomb.tw` example.

.. _args:

Command line arguments
----------------------

For desktop installations the command line specification is

.. py:function:: tw3d [-n <procs>] [-c <threads>] [--input-file <file>] [--no-interactive] [--restart] [--version] [--help]

	:param int procs: number of MPI processes (default=1, desktop only)
	:param int threads: number of OpenMP threads (see below for default)
	:param str file: name or path of the file to use as the input file (default=stdin)

	The :samp:`--restart` argument, if present, causes initial data to be loaded from a checkpoint.

	The :samp:`--no-interactive` argument, if present, suppresses the interactive thread.

	The :samp:`--version` argument, if present, prints the version number.  If this is the only argument, no simulation is attempted.

	The :samp:`--help` argument, if present, prints the command line arguments and the link to the online documentation.  If this is the only argument, no simulation is attempted.

If you enter only :samp:`tw3d` with no arguments, turboWAVE will use a single MPI processes, and will fork as many threads as there are logical cores on the system.  If you enter :samp:`tw3d -n {procs}`, turboWAVE will use the requested number of MPI processes, but only a single thread.  Finally, if you enter :samp:`tw3d -n {procs} -c {threads}`, turboWAVE will use the requested number for both processes and threads.

When you ran the example above, you may have noticed turboWAVE issuing a warning about the domain decomposition.  That is because if you choose to specify the domain decomposition in the input file, the product of the three integers is supposed to equal the number of processes requested.  If this is not the case, turboWAVE will try to find a suitable decomposition on its own.  There are some rules about how this can be done.  Sometimes turboWAVE will fail to find a suitable decomposition and report an error.

Finally, if you want to disable the interactive thread, add the command line argument :samp:`--no-interactive`.  This can be important for batch processing, because when the interactive thread is used, the :samp:`tw3d` process will not stop without a keystroke from the user.

Error Handling
--------------

It is important to pay attention to the output file if you are having problems.  If the code stops without reporting an error in the terminal window, you may still be able to get some feedback.  The procedure is as follows.

	#. In the input file, add the line :samp:`output level = 1`
	#. This line can go anywhere except within a :samp:`new` block or :samp:`generate` block
	#. Run the problem again
	#. If the error is not reported on the console, try :samp:`grep ERROR *stdout*`
