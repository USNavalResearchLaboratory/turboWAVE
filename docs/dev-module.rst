Module
======

The ``Module`` class is the highest level object for managing simulations.

Communication between modules can be accomplished in two ways.

	1. Publisher-consumer mechanism
	2. Quasi-flat containment hierarchy

The preferred turboWAVE style is to use the horizontal publisher-consumer model, with modules containing only lower level constructs.  Ideally these lower level constructs are formal ``ComputeTool`` objects, but there are still elements (as of this writing) that fall outside this convention.  In order to handle these there is also a quasi-tool protocol.

In principle the publisher-consumer model can handle most of the inter-module relationships one needs, but in some cases a containment structure is irresistible.  The longest-standing such relationship is the containment of ``Species`` modules within the ``Kinetics`` module.  The SPARC hydrodynamics modules actually have two levels of containment, the ``Chemistry`` module contains ``EquilibriumGroup`` modules, and the ``EquilibriumGroup`` modules contain ``Chemical`` modules.

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

The turboWAVE containment hierarchy is quasi-flat, because the ``Grid`` object owns every module and stores references to them on a flat list. It is useful to define terms as follows.

.. glossary::

	Ownership
		When an object owns another object, it has the exclusive right to create and release that object.

	Supermodule
		Any module which has a non-empty ``submodule`` container.  The ``submodule`` container is a flat list of references to other modules.

	Submodule
		Any module which has a ``super`` pointer with a value other than ``NULL``.  The ``super`` pointer must point to a module whose ``submodule`` container includes a reference to the referencing submodule (the supermodule and submodule must point to each other).

	Singular Module
		Any module which requires that there be only one instance of it per MPI task.
