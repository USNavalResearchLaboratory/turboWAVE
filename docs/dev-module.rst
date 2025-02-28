Module
======

The ``Module`` class is the highest level object for managing simulations.  Its principle method is the virtual function ``Update``.  The simulation cycle invokes ``Update`` for each module.  The order of ``Update`` calls is determined by a sorting index stored with each ``Module`` instance.

Communication between modules can be accomplished in two ways.

	1. Publisher-consumer mechanism
	2. Containment hierarchy

The preferred turboWAVE style is to use the horizontal publisher-consumer model, with modules containing only lower level constructs.  Ideally these lower level constructs are formal ``ComputeTool`` objects.

In principle the publisher-consumer model can handle most of the inter-module relationships one needs, but in some cases a containment structure is irresistible.  The longest-standing such relationship is the containment of ``Species`` modules within the ``Kinetics`` module.  The SPARC hydrodynamics modules actually have two levels of containment, the ``HydroManager`` module contains ``EquilibriumGroup`` modules, and the ``EquilibriumGroup`` modules contain ``Chemical`` modules.

.. _pub-cons:

Publisher-Consumer Model
------------------------

In the publisher-consumer model, derived modules override the ``InspectResource`` and ``ExchangeResources`` methods.  Typically ``InspectResource`` copies pointers passed in by other modules:

	.. cpp:function:: bool InspectResource(void* resource,const std::string& description)

		This is the function which consumes resources provided by other modules.

		:param resource: Pointer from another module to be copied to a member of the inspecting module.
		:param description: String used to identify the resource.
		:return: Whether the resource was copied or not
		:rtype: bool

The body of the function will test several string literals against ``description``, and if a match is found copy ``resource`` to an appropriate member pointer.

The ``ExchangeResources`` method calls the ``PublishResource`` method for each pointer the module wants to share.  The ``PublishResource`` method, which does not need to be overridden, simply calls ``InspectResource`` for every module (note that nothing stops a module from consuming its own resource).

	.. cpp:function:: void ExchangeResources()

		This function is intended to be populated with one or more calls to ``PublishResource``.

	.. cpp:function:: void PublishResource(void* resource,const std::string& description)

		Shares a resource with all modules. Does not need to be overridden.

		:param resource: Pointer to share with other modules.
		:param description: String used to identify the resource, typically a literal.


Containment Model
-----------------

The turboWAVE containment hierarchy uses both a flat list and a tree structure.  The ``Simulation`` object is the root. ``Simulation`` owns every ``Module`` and stores references to them on a flat list.  The flat list is sorted to reflect the structure, with objects nearer the root listed first.

The ``Module`` objects use a simple tree structure.  Each ``Module`` has references to other containment hierarchy elements as follows.

	#. ``super`` - a single pointer to a supermodule, which is ``NULL`` if the module is directly below ``Simulation``.
	#. ``owner`` - a single pointer to ``Simulation``.
	#. ``submodule`` - a ``std::vector`` of pointers to submodules.
	#. ``moduleTool`` - a ``std::vector`` of pointers to ``ComputeTool`` objects the user associated with the module.

The containment hierarchy is created while reading the input file.

Implementing a Module
---------------------

Declaration
,,,,,,,,,,,

When implementing a new ``Module``, first carry out the following.

	#. In ``Module.cpp``, introduce a new ``tw::module_type`` element.  This is as a label for the type of module.
	#. In ``Factory.cpp``, add a case to the function ``CreateModuleFromType`` for the new type.
	#. If this is a singular module, modify the static member ``SingularType`` in ``Module.cpp`` appropriately.
	#. In an appropriate header file, derive the new type from ``Module``.
	#. In an appropriate source file, implement the ``Module``.

		* Almost every module will override the ``Update`` method.
		* Further implementation details follow.

Constructor and Initializer
,,,,,,,,,,,,,,,,,,,,,,,,,,,

The constructor is used to define input file parameters (see below), to set default parameter values, and to set up allocations that are independent of subsequent input file processing.

Values or allocations that can only be known after the whole input file is processed should be set in the ``Initialize`` method.

Input File Support
,,,,,,,,,,,,,,,,,,

With very little effort the user will be able to create the module and associate it with other objects from within the input file.  To support this, carry out the following steps.

	#. In the module's constructor define the input file directives. For each directive make one call to ``directives.Add(std::string&,tw::input::Directive*)``.
	#. Add an entry to the hash table returned by ``Map`` in ``Module.cpp``.  This connects the input file keys with the ``tw::module_type``.

Containment Support
,,,,,,,,,,,,,,,,,,,

Input file semantics automatically establish the containment tree.  However, there are some details of the relationship that have to be specified in source code.

#. If your module is intended as a supermodule:

	* If you need strongly typed pointers to submodules, use ``Module::VerifyInput`` to search the ``submodule`` vector for the desired modules.  Use ``dynamic_cast`` to identify the module type, and to create the strongly typed pointer.

#. If your module is intended as a submodule:

	* If the submodule *requires* its supermodule, add an entry to the hash table in the static member ``RequiredSupermoduleType``.

#. If you want your module to use the ``ComputeTool`` system, see :doc:`dev-tool`.

Intermodule Processing
,,,,,,,,,,,,,,,,,,,,,,

If your module needs to share data through the publisher-consumer mechanism, follow the guidance :ref:`above <pub-cons>`.  If you want to use the containment hierarchy to orchestrate more complex interactions, you may want to store explicitly typed pointers to the supermodule or certain submodules.  This should be done in the ``VerifyInput`` method.

Restart File Support
,,,,,,,,,,,,,,,,,,,,

As of version 4.0, only time varying quantities need to be checkpointed (no need to store constants or structural information).  To support restarting a module, carry out the following steps.

	#. Override the ``ReadCheckpoint`` method.  Call the superclass ``ReadCheckpoint`` method first.  Then read any necessary data from the restart file.
	#. Override the ``WriteCheckpoint`` method.  Call the superclass ``WriteCheckpoint`` method first.  Then write any necessary data to the restart file.
	#. Verify that ``ReadCheckpoint`` and ``WriteCheckpoint`` access the data in the same order.

Glossary
--------

.. glossary::

	Ownership
		When an object owns another object, it has the exclusive right to create and release that object.

	Supermodule
		Any module which has a non-empty ``submodule`` container.  The ``submodule`` container is a flat list of references to other modules.

	Submodule
		Any module which has a ``super`` pointer with a value other than ``NULL``.  The ``super`` pointer must point to a module whose ``submodule`` container includes a reference to the referencing submodule (the supermodule and submodule must point to each other).

	Singular Module
		Any module which requires that there be only one instance of it per MPI task.

Best Practices
--------------

#. Do not use upper case in defining your input file keys or directives.
#. Use descriptive English language keys and directives, but without excessive verbosity.
