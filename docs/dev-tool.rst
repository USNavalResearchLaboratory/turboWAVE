ComputeTool
===========

The ``ComputeTool`` object is used to perform calculations that may be re-used by various modules.  Unlike a ``Module``, a ``ComputeTool`` is assumed to perform its functions independently, and so does not require any data sharing paradigm beyond passing arguments into its methods.  All ``ComputeTool`` objects are owned by ``Simulation``.  When a ``Module`` wants to create or release a ``ComputeTool`` it calls a method of ``Simulation`` to do so.

When to use a ``ComputeTool``
-----------------------------

The ``ComputeTool`` provides services for free that a simple function or user defined object does not.  There are also services that are not provided, being reserved for ``Module`` objects.  The basic characteristics are:

	#. Easy management of input file and restart file interactions
	#. Easy access to OpenCL kernel functions.
	#. Management of allocation and de-allocation of resources, accounting for the possibility of multiple references to the tool.
	#. No data sharing.  A ``ComputeTool`` does not participate in publisher-consumer data movement, and one ``ComputeTool`` cannot contain another.
	#. No direct invocation within the main simulation loop.  All data processing has to be initiated by a ``Module`` which has a pointer to the ``ComputeTool``.

Example
-------

A straightforward example of a module that uses a ``ComputeTool`` inheritance tree is the ``FieldSolver`` module in ``fieldSolver.h`` and ``FieldSolver.cpp``.  This ``FieldSolver`` base class exists largely to manage various types of elliptical solvers.  The elliptical solver tools themselves are implemented in ``elliptic.h`` and ``Elliptic.cpp``.

Implementing a ``ComputeTool``
------------------------------

Declaration
,,,,,,,,,,,

When implementing a new ``ComputeTool``, first carry out the following.

	#. In ``ComputeTool.cpp``, introduce a new ``tw::tool_type`` element.  This is a label for the type of tool.
	#. In ``Factory.cpp``, add a case to the ``CreateToolFromType`` function for the new type.
	#. In an appropriate header file, derive the new type from ``ComputeTool``.
	#. In an appropriate source file, implement the ``ComputeTool``.

		* The tool can implement whatever methods are desired. A function that moves data forward by one time level is conventionally called ``Advance``.
		* Further implementation details follow.

Association with ``Module``
,,,,,,,,,,,,,,,,,,,,,,,,,,,

Carry out the following for any ``Module`` that wants to use the tool.

	#. In the ``Module`` declaration, declare a pointer to the tool.  Polymorphism may be used to whatever extent is desired.
	#. In the constructor, set the pointer to ``NULL``.
	#. In the destructor, if the pointer is not ``NULL``, call ``owner->RemoveTool`` with the pointer as the argument.

Input File Support
,,,,,,,,,,,,,,,,,,

If you want the tool to be accessible from the input file, carry out the following steps.

	#. In the tool's constructor define the input file directives. For each directive make one call to ``directives.Add(std::string&,tw::input::Directive*)``.
	#. Add an entry to the hash table returned by ``Map`` in ``ComputeTool.cpp``.  This connects the directive key with the ``tw::tool_type``.
	#. In the ``VerifyInput`` method of each ``Module`` that wants to use the tool:

		* Search the ``moduleTool`` vector for a compatible tool, and copy the dynamically typecast compatible pointer to your pointer.
		* If no compatible tool is found, either throw an error, or create a default tool using ``owner->CreateTool``.

	#. If you need to set the tool's member variables based on module data, it is safest to do this in the ``Initialize`` method, since this is called only after all modules have exchanged resources.

.. tip::

	If you want to create a tool exclusively for the use of a particular module, *and* there is no need for input file or restart file interaction, you can simply create it in the module constructor, and remove it in the module destructor.

.. tip::

	If implemented correctly, tools are ready to be used by the time ``Module::Initialize`` is called.

Restart File Support
,,,,,,,,,,,,,,,,,,,,

As of version 4.0, only time varying quantities need to be checkpointed (no need to store constants or structural information).  To support restarting a tool, carry out the following steps.

	#. Override the tool's ``ReadCheckpoint`` method.  Call the inherited ``ReadCheckpoint`` method first.  Then read any necessary data from the restart file.
	#. Override the tool's ``WriteCheckpoint`` method.  Call the inherited ``WriteCheckpoint`` method first.  Then write any necessary data to the restart file.
	#. Verify that ``ReadCheckpoint`` and ``WriteCheckpoint`` access the data in the same order.

Best Practices
--------------

#. Avoid making the tool the owner of heavyweight data.  Instead pass such data to member functions by reference.
#. Keep self-contained, i.e., avoid using references to modules or other objects in the containment hierarchy.  If this seems unavoidable consider using a ``Module``.
