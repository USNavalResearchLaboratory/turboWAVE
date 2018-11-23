ComputeTool
===========

The ``ComputeTool`` object is used to perform calculations that may be re-used by various modules.  Unlike a ``Module``, a ``ComputeTool`` is assumed to perform its functions independently, and so does not require any data sharing paradigm beyond passing arguments into its methods.  All ``ComputeTool`` objects are owned by ``Grid``.  When a ``Module`` wants to create or release a ``ComputeTool`` it calls a method of ``Grid`` to do so.

.. note::

	While ``Module`` objects have been part of turboWAVE almost from the inception, ``ComputeTool`` objects were introduced later.  As a result there are objects in the framework that over time may be transitioned to ``ComputeTool`` objects.

When to use a ``ComputeTool``
-----------------------------

The ``ComputeTool`` provides services for free that a simple function or user defined object does not.  There are also services that are not provided, being reserved for ``Module`` objects.  The basic characteristics are:

	#. Standard structure for managing input file and restart file interactions
	#. Easy access to OpenCL kernel functions.
	#. Management of allocation and de-allocation of resources, accounting for the possibility of multiple references to the tool.
	#. No data sharing.  A ``ComputeTool`` does not participate in publisher-consumer data movement, and one ``ComputeTool`` cannot contain another.
	#. No direct invocation within the main simulation loop.  All data processing has to be initiated by a ``Module`` which has a pointer to the ``ComputeTool``.

Declaration
-----------

When implementing a new ``ComputeTool``, first carry out the following.

	#. In ``computeTool.h``, introduce a new ``tw::tool_type`` element to identify the type of tool.
	#. In ``ComputeTool.cpp``, add a case to the ``CreateObjectFromType`` method for the new type.
	#. In an appropriate header file, derive the new type from ``ComputeTool``.
	#. In an appropriate source file, implement the ``ComputeTool``.

		* The tool can implement whatever methods are desired. A function that moves data forward by one time level is conventionally called ``Advance``.
		* Further implementation details follow.

Allocation
----------

Carry out the following for any ``Module`` that wants to use the tool.

	#. In the ``Module`` declaration, declare a pointer to the tool.  Use polymorphism in whatever degree is preferable.
	#. In the constructor, set the pointer to ``NULL``.
	#. In the destructor, if the pointer is not ``NULL``, call ``owner->RemoveTool`` with the pointer as the argument.
	#. In the ``Initialize`` function, you can test to see if the pointer is ``NULL``, and create a default tool, if desired.  To create the tool, use ``owner->CreateTool``.

Input File Support
------------------

There are two ways to create tools from the input file.

	1. Tools may be created at the root level and given an explicit name.  Modules then access the tool by name.
	2. Tools may be created on the fly from within a module block using a directive.  A unique name is assigned automatically.

If you want the tool to be accessible from the input file, carry out the following steps.

	#. Implement the ``ReadInputFileDirective`` method for the tool.
	#. If you want to create named tools, add a case to ``CreateTypeFromInput`` in ``ComputeTool.cpp``.
	#. If you want to create tools on the fly, add a case to ``CreateTypeFromDirective`` in ``ComputeTool.cpp``.
	#. In the ``ReadInputFileDirective`` method of each ``Module`` that wants to use the tool:

		* call ``owner->ToolFromDirective`` and assign the resulting pointer to the module's pointer to the tool.
		* Test the pointer, and if valid call the ``ReadInputFileDirective`` of the tool.

Restart File Support
--------------------

To support restarting a tool, carry out the following steps.

	#. Override the ``ReadData`` method.  Call the superclass ``ReadData`` method first.  Then read any necessary data from the restart file.
	#. Override the ``WriteData`` method.  Call the superclass ``WriteData`` method first.  Then write any necessary data to the restart file.
	#. Verify that ``ReadData`` and ``WriteData`` access the data in the same order.
	#. For any Module that uses the tool:

		* In the module's ``ReadData`` function, call ``owner->GetRestartedTool`` and save the returned pointer to a member of the module.
		* In the module's ``WriteData`` function, call ``SaveToolReference``, accessing through the pointer to the tool.
		* The two above calls must occur at the same point in the restart file.
