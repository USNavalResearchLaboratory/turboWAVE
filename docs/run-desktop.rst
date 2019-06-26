Desktop Runs
============

.. caution::

	This document assumes you have followed the installation instructions precisely.

.. note::

	Once installed, running turboWAVE should be very nearly the same for any desktop system. For Windows, this consistency depends on using the PowerShell as the terminal window.


Parallel Programming Theory
---------------------------

In order to run turboWAVE, it is useful to understand a little about the parallelization methods.  TurboWAVE uses a combination of distributed memory and shared memory methods.  The distributed memory method is called "domain decomposition", and corresponds to physically partitioning the simulation region into chunks that can be worked on by different processors.  The memory that holds each chunk is distributed in the sense that it may be connected only by a network cable.  The shared memory method is called "fork-join".  In this approach independent software threads are running on different cores, but in a setting where each core has equal access to the memory (e.g., the threads could be on different cores of the same multi-core processor).

The distributed memory model is implemented using software called MPI.  The shared memory model is implemented using software called OpenMP.

As the user, you have to choose how to partition the problem.  In particular, you must choose the number of MPI processes, and the number of threads.  If you opt to have more than one MPI process, you may also want to control exactly how the physical simulation domain is partitioned into chunks.

Running an Example
------------------

#. Pick some example from :samp:`{twroot}/core/examples`.  For this test you should avoid the 3D examples.
#. For definiteness, let us use :samp:`{twroot}/core/examples/pgc/LWFA-coulomb.txt`
#. Open a terminal window and navigate to :samp:`~/Run`
#. :samp:`cp {twroot}/core/examples/pgc/LWFA-coulomb.txt stdin`
#. This puts the input file in :samp:`~/Run` with the name :samp:`stdin`.  By default, turboWAVE assumes the input file is in the working directory with the name :samp:`stdin`.
#. :samp:`tw3d -n 4`
#. The above command runs the problem with 4 MPI proceses and 1 thread per process.  Of course this choice may not be optimal for your system, method of compiling, etc., but it should suffice for this example.
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

	Fig. 1 --- Python DataViewer output of the scalar potential produced by the :file:`LWFA-coulomb.txt` example.

.. _args:

Command line arguments
----------------------

The full command line specification is

.. py:function:: tw3d -n procs -c threads --input-file file --no-interactive --version --help

	:param int procs: number of MPI processes (default=1)
	:param int threads: number of OpenMP threads (see below for default)
	:param str file: name or path of the file to use as the input file (default=stdin)

	The :samp:`--no-interactive` argument, if present, suppresses the interactive thread.

	The :samp:`--version` argument, if present, prints the version number.  If this is the only argument, no simulation is attempted.

	The :samp:`--help` argument, if present, prints out a message pointing to the online documentation.  If this is the only argument, no simulation is attempted.

Except for :samp:`--input-file`, these arguments are only used on the desktop.  If you enter only :samp:`tw3d` with no arguments, turboWAVE will use a single MPI processes, and will fork as many threads as there are logical cores on the system.  If you enter :samp:`tw3d -n {procs}`, turboWAVE will use the requested number of MPI processes, but only a single thread.  Finally, if you enter :samp:`tw3d -n {procs} -c {threads}`, turboWAVE will use the requested number for both processes and threads.

When you ran the example above, you may have noticed turboWAVE issuing a warning about the domain decomposition.  That is because if you choose to specify the domain decomposition in the input file, the product of the three integers is supposed to equal the number of processes requested.  If this is not the case, turboWAVE will try to find a suitable decomposition on its own.  There are some rules about how this can be done.  Sometimes turboWAVE will fail to find a suitable decomposition and report an error.

Finally, if you want to disable the interactive thread, add the command line argument :samp:`--no-interactive`.  This can be important for batch processing, because when the interactive thread is used, the :samp:`tw3d` process will not stop without a keystroke from the user.

Error Handling
--------------

It is important to pay attention to the output file if you are having problems.  If the code stops without reporting an error in the terminal window, you may still be able to get some feedback.  The procedure is as follows.

	#. In the input file, add the line :samp:`stdout = full`
	#. This line can go anywhere except within a :samp:`new` block or :samp:`generate` block
	#. Run the problem again
	#. If the error is not reported on the console, try :samp:`grep ERROR *stdout*`

Test Suite
----------

The example input files comprise a test suite for the turboWAVE installation.  There is a Python script :samp:`tools/twtest/twtest.py` which automatically runs all of the example cases and generates a report.  The report contains animations produced using the ImageMagick suite's :samp:`convert` program.  Check to see if it is installed using

:kbd:`convert --version`

If you don't see the ImageMagick version number displayed, you must install it.  It should be installable with any package manager (apt, yum, homebrew, etc.).

In order to run the script navigate to :samp:`tools/twtest` and invoke

:kbd:`python twtest.py` *twroot* *args*

where *twroot* is the turboWAVE root directory path and *args* are the usual command line arguments used to specify parallelism options (you do not need to add :samp:`--no-interactive` as this is put in automatically).  Due to the large number of simulations to be run this may take several hours.  You can limit the test to specific categories by appending them to *twroot* using double colon separators.  For example,

:kbd:`python twtest.py ~/turboWAVE::hydro::pic -n 4`

would test all the examples in the ``hydro`` and ``pic`` directories using a four-way domain decomposition.

.. note::

	The :samp:`twtest.py` script will try to adjust the parallel parameters between 1D and 2D examples.  A safe choice is to use 4 MPI processes and enough threads per process to occupy all the cores on your system.  Choices that cause the script to fail are possible.

.. warning::

	The :samp:`twtest.py` script assumes the standard turboWAVE directory structure has not been disturbed.  The script freely deletes files in the :samp:`twtest` directory during cleanup operations.

There is a special comment line in most of the example files that triggers execution by :samp:`twtest.py`.  The form is

:samp:`// TWTEST matplotlib` *slicing_spec=slices* *file* *dynamic_range*

The ``TWTEST`` token tells ``twtest.py`` to process this line.  The ``matplotlib`` token indicates that the ``matplotlib`` library will be used to make the figure to include in the report.  The *slicing_spec* is replaced by a four character ordered list of axes, such as ``zxyt``.  The first two characters are the axes to plot, the last two are sliced using the indices given by *slices*.  These are separated by a comma.  An example of a full slicing specification is ``zxyt=0,-1``.  Note that the index -1 is special, indicating that an animation over all frames should be created.  Other negative indices work in the usual Python manner.  The *file* is the name of the file to plot (do not include path).  The *dynamic_range* is a floating point number, which, if zero, leads to a linear scale, if non-zero, leads to a log-plot with the requested dynamic range.

.. note::

	There are going to be speckles in the movie images due to GIF compression.

When the script completes there should be a file called :samp:`twreport.html`.  Open this in your browser to examine the results.  There should be a heading for each example subdirectory and a subheading for each example.  If the run failed any error messages are recorded.  If it succeeded an image or animation showing the data that was produced is displayed.

.. tip::

	If you would like to check on the progress of a particular run that has been executed by the script, open a separate terminal window, navigate to the :samp:`tools/twtest` directory, and type :kbd:`cat twstat`.

.. tip::

	If you would like to "comment out the comment", e.g., to skip over the longer examples, change ``TWTEST`` to lower case.
