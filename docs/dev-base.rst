Basic Architecture
==================

MPI Launch
-----------

At the most basic level, a turboWAVE simulation consists of multiple independently running MPI processes, which must communicate with each other.  The launch procedure takes two possible forms:

	1. Internal MPI launch
	2. External MPI launch

Internal MPI launch
,,,,,,,,,,,,,,,,,,,

If the code is compiled usng the internal implementation of MPI, then the MPI processes are treated as threads.  These threads are created and manipulated using the standard C++ library.

In the internal launch mode, MPI processes are launched by a master thread from within the ``tw3d`` executable.  The master thread is encapsulated in a ``Launcher`` class defined in ``Main.cpp``.  The ``Launcher`` inherits from ``tw::Thread``, which is a wrapper for ``std::thread``.

External MPI launch
,,,,,,,,,,,,,,,,,,,

If the code is compiled using an external MPI library, then the MPI processes are under the control of that library, but are typically distinct operating system processes, as opposed to threads running in the same process.

In the external launch mode some third party program creates multiple instances of the ``tw3d`` executable.  In this case explicit threads are not used.  Implicit OpenMP threads may be forked by each process, however.

Simulation Class
----------------

Regardless of the launch mechanism, it results in the creation of a ``Simulation`` object for each MPI process.  The ``Simulation`` object is the master container object.  It inherits from two lower level classes, ``Task`` and ``MetricSpace``

``Task`` is a container for various MPI communicators.  It manages all kinds of communication between MPI processes.  It contains information on the structure of the domain decomposition.  ``MetricSpace`` defines the geometry of the grid cells.

It is important to remember that each MPI process has its own unique instance of the ``Task`` object, even though this object describes the entire domain.  For example, the ``Task`` object stores the ranks of neighboring domains.  This data is different on each running process.

Immediately upon launch, ``Simulation`` makes a first pass through the input file using the ``InputFileFirstPass`` method, and then executes the ``Run`` method.  The ``InputFileFirstPass`` method is used primarily to setup the domain decomposition using information in the input file.

The ``Run`` method executes the following sequence:

	#. Create status file
	#. Execute ``PrepareSimulation`` method

		#. Construct grid using input file
		#. Construct other objects using input file
		#. Initialize ``Region`` objects
		#. Initialize ``Wave``, ``Pulse``, and ``Conductor`` objects
		#. Initialize ``ComputeTool`` objects
		#. Initialize ``Module`` objects

	#. Loop over ``FundamentalCycle`` method

		#. Write diagnostics
		#. Reset all modules
		#. Update all modules
		#. Advance all clocks
		#. Modify the grid

	#. Close status file

.. Tip::

	In turboWAVE constructing and initializing objects are two different operations.  All objects are constructed while processing the input file.  Initialization of objects happens only after all objects are constructed.

.. Warning::

	TurboWAVE's internal MPI threads use serial I/O.  Therefore message passing is illegal in an object's constructor, because it may be executed during I/O operations.  You may safely use the object's ``Initialize`` method instead.
