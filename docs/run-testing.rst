Testing
==========

The test runner, ``twtest``, is a python program that can invoke three types of tests, unit tests, integration tests, and "sea trials".  The ``twtest`` runner is installed with ``twutils`` and should be available from your conda environment.  To see the command line options type::

	twtest -h

Unit Testing
------------

The unit tests can be run directly from the executable using the command line option ``--unit-test``.  They can also be invoked from the test runner ``twtest``.  Example invocation::

	tw3d -n 2 --unit-test --all

Any given unit test is allowed to reject the environment (number of MPI nodes, etc.) with a warning message.

Integration Tests
-----------------

The integration tests are very small simulations that should be thought of as tests of the solution of a complex set of difference equations.  These are not tests of physical accuracy, because the resolution is typically very low.  The input files for the integration tests are found in :samp:`{turbowave}/core/test`.  Each integration test input file has a special comment line describing the test.  For example::

	// CITEST Ex.npy range=-.024,.024 tolerance=1e-3

would check that the electric field x-component minimum is -0.024, and that the maximum is 0.024, with a tolerance of 0.001.  Example invocation::

	twtest --integration --root ~/path/to/turboWAVE --command tw3d -n 4

Sea Trials
----------

The example input files comprise a test suite that verifies physics as well as differencing.  The ``twtest`` script can run all of the example cases and generate a report.  Each example input file has a special comment line describing the figure to generate for the report.  For example::

	// TWTEST matplotlib` *slicing_spec=slices* *file* *key=val*...

The ``TWTEST`` token tells ``twtest`` to process this line.  The ``matplotlib`` token indicates that the ``matplotlib`` library will be used to make the figure to include in the report.  The *slicing_spec* is replaced by a four character ordered list of axes, such as ``zxyt``.  The first two characters are the axes to plot, the last two are sliced using the indices given by *slices*.  These are separated by a comma.  An example of a full slicing specification is ``zxyt=0,:``.  The colon is a special slice that indicates that an animation over all frames should be created.  Negative slices are interpreted in the usual Python manner.  The *file* is the name of the file to plot (do not include path).  Trailing arguments are optional key/value pairs separated by equals signs (no intervening spaces).

The sea trials can be invoked as follows::

	twtest --sea-trials --root ~/path/to/turboWAVE --command tw3d -n 16

Due to the large number of simulations to be run this may take a few hours.  You can limit the test to specific categories using the ``--categories`` option::

	twtest --sea-trials --root ~/path/to/turboWAVE --categories pic,hydro --command tw3d -n 16

This would test all the examples in the ``hydro`` and ``pic`` directories.

If turboWAVE is compiled against an external MPI library, simply substitute the appropriate command::

	twtest --sea-trials --root ~/path/to/turboWAVE --categories pic,hydro --command mpirun -np 16 tw3d -c 2

.. note::

	The :samp:`twtest` script will try to adjust the parallel parameters between 1D and 2D examples.  A safe choice is to use 4 MPI processes and enough threads per process to occupy all the cores on your system.  Choices that cause the script to fail are possible.

.. warning::

	The :samp:`twtest` script assumes the standard turboWAVE directory structure has not been disturbed.  The script freely deletes files in the working directory during cleanup operations.

When the script completes there should be a file called :samp:`twreport.html`.  Open this in your browser to examine the results.  There should be a heading for each example subdirectory and a subheading for each example.  If the run failed any error messages are recorded.  If it succeeded an image or animation showing the data that was produced is displayed.

.. tip::

	If you would like to check on the progress of a particular run that has been executed by the script, open a separate terminal window, navigate to the ``twtest`` working directory, and type :kbd:`cat twstat`.

.. tip::

	If you would like to "comment out the comment", e.g., to skip over the longer examples, change ``TWTEST`` to lower case.
