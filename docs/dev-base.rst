Basic Architecture
==================

MPI Launch
-----------

At the most basic level, a turboWAVE simulation consists of multiple independently running MPI processes, which must communicate with each other.  Prior to v5 there was an internal MPI that could be used in place of an external library, but this has been discarded, since at this point, MPI libraries are easy to install on most platforms.

As a result, turboWAVE will always be launched with ``mpirun``, ``mpiexec``, or other MPI launch command.

Simulation Class
----------------

TurboWAVE will create a ``Simulation`` object for each MPI process.  The ``Simulation`` object is the master container object.  It inherits from two lower level classes, ``Task`` and ``MetricSpace``

``Task`` is a container for various MPI communicators.  It manages all kinds of communication between MPI processes.  It contains information on the structure of the domain decomposition.  ``MetricSpace`` defines the geometry of the grid cells.

It is important to remember that each MPI process has its own unique instance of the ``Task`` object, even though this object describes the entire domain.  For example, the ``Task`` object stores the ranks of neighboring domains.  This data is different on each running process.

Immediately upon launch, ``Simulation`` makes a first pass through the input file using the ``InputFileFirstPass`` method, and then executes the ``Run`` method.  The ``InputFileFirstPass`` method is used primarily to setup the grid and the domain decomposition using information in the input file.

The ``Run`` method executes the following sequence:

	#. Create status file
	#. Execute ``PrepareSimulation`` method

		#. Construct all objects using input file
		#. Verify all ``Module`` objects
		#. Initialize ``Region`` objects
		#. Initialize ``ComputeTool`` objects
		#. Exchange resources among ``Module`` objects
		#. Initialize ``Module`` objects

	#. Loop over ``FundamentalCycle`` method

		#. Write diagnostics
		#. Reset all modules
		#. Update all modules
		#. Advance all clocks
		#. Modify the grid

	#. Close status file

.. Tip::

	In turboWAVE constructing and initializing objects are two different operations.  All objects are constructed while processing the input file.  Initialization of objects happens only after all objects are constructed.  Modules also have an intermediate verification stage.

.. Warning::

	TurboWAVE's internal MPI threads use serial I/O.  Therefore message passing is illegal in an object's constructor, because it may be executed during I/O operations.  You may safely use the object's ``Initialize`` method instead.
