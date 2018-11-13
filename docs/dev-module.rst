Module
======

The ``Module`` class is the highest level object for managing simulations.

Communication between modules can be accomplished in two ways.

	1. Publisher-consumer mechanism
	2. Quasi-flat containment heierarchy

The preferred turboWAVE style is to use the horizontal publisher-consumer model, with modules containing only lower level constructs.  Ideally these lower level constructs are formal ``ComputeTool`` objects, but there are still elements (as of this writing) that fall outside this convention.  In order to handle these there is also a quasi-tool protocol.

In principle the publisher-consumer model can handle most of the inter-module relationships one needs, but in some cases a containment structure is irresistible.  The longest-standing such relationship is the containment of ``Species`` modules within the ``Kinetics`` module.  The SPARC hydrodynamics modules actually have two levels of containment, the ``Chemistry`` module contains ``EquilibriumGroup`` modules, and the ``EquilibriumGroup`` modules contain ``Chemical`` modules.

Publisher-Consumer Model
------------------------

Containment Model
-----------------

The turboWAVE containment heierarchy is quasi-flat, because the ``Grid`` object owns every module and stores references to them on a flat list. Explaining the nature of the containment model requires defining terms as follows.

.. glossary::

	Ownership
		When an object owns another object, it has the exclusive right to create and release that object.

	Supermodule
		Any module which has a non-empty ``submodule`` container.  The ``submodule`` container is a flat list of references to other modules.

	Submodule
		Any module which has a ``super`` pointer with a value other than ``NULL``.  The ``super`` pointer must point to a module whose ``submodule`` container includes a reference to the referencing submodule (the supermodule and submodule must point to each other).

	Singular Module
		Any module which requires that there be only one instance of it per MPI task.
