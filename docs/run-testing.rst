Test Suite
==========

The example input files comprise a test suite for the turboWAVE installation.  There is a Python script ``twtest`` which automatically runs all of the example cases and generates a report.

In order to run the script activate your environment and navigate to a directory where you would like to generate the report.  Then execute

:kbd:`twtest <twroot> <cmd>`

where ``<twroot>`` is the turboWAVE root directory path and ``<cmd>`` is the command that the script will invoke to execute each simulation in the test suite.  The command line arguments :samp:`--no-interactive` and :samp:`--input-file` are added automatically, and should therefore be omitted.  Due to the large number of simulations to be run this may take a few hours.  You can limit the test to specific categories by appending them to ``<twroot>`` using double colon separators.  For example,

:kbd:`twtest ~/turboWAVE::hydro::pic tw3d -n 4`

would test all the examples in the ``hydro`` and ``pic`` directories using a four-way domain decomposition.  On the other hand,

:kbd:`twtest ~/turboWAVE tw3d -n 4`

would run the test cases in every directory.  If turboWAVE is compiled against an external MPI library, simply substitute the appropriate command, e.g.,

:kbd:`twtest ~/turboWAVE mpirun -np 4 tw3d -c 2`

would use OpenMPI to launch 4 processes, and OpenMP to fork 2 threads per process.

.. note::

	The :samp:`twtest` script will try to adjust the parallel parameters between 1D and 2D examples.  A safe choice is to use 4 MPI processes and enough threads per process to occupy all the cores on your system.  Choices that cause the script to fail are possible.

.. warning::

	The :samp:`twtest` script assumes the standard turboWAVE directory structure has not been disturbed.  The script freely deletes files in the working directory during cleanup operations.

There is a special comment line in most of the example files that triggers execution by :samp:`twtest`.  The form is

:samp:`// TWTEST matplotlib` *slicing_spec=slices* *file* *key=val*...

The ``TWTEST`` token tells ``twtest`` to process this line.  The ``matplotlib`` token indicates that the ``matplotlib`` library will be used to make the figure to include in the report.  The *slicing_spec* is replaced by a four character ordered list of axes, such as ``zxyt``.  The first two characters are the axes to plot, the last two are sliced using the indices given by *slices*.  These are separated by a comma.  An example of a full slicing specification is ``zxyt=0,:``.  The colon is a special slice that indicates that an animation over all frames should be created.  Negative slices are interpreted in the usual Python manner.  The *file* is the name of the file to plot (do not include path).  Trailing arguments are optional key/value pairs separated by equals signs (no intervening spaces).

.. note::

	There are going to be speckles in the movie images due to GIF compression.

When the script completes there should be a file called :samp:`twreport.html`.  Open this in your browser to examine the results.  There should be a heading for each example subdirectory and a subheading for each example.  If the run failed any error messages are recorded.  If it succeeded an image or animation showing the data that was produced is displayed.

.. tip::

	If you would like to check on the progress of a particular run that has been executed by the script, open a separate terminal window, navigate to the ``twtest`` working directory, and type :kbd:`cat twstat`.

.. tip::

	If you would like to "comment out the comment", e.g., to skip over the longer examples, change ``TWTEST`` to lower case.
