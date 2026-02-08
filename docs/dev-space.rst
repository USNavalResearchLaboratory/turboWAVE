Static Space
==============

The ``StaticSpace`` object is the basis of structured grid management in turboWAVE. It is the root of the space-time mesh inheritance tree.

* ``StaticSpace`` - basic topology and indexing that stays constant throughout a simulation
* ``DynSpace`` - adds coordinates or parameters that may evolve during a simulation
* ``MetricSpace`` - adds metric information geared toward finite volume methods

Owing to the nature of particle-in-cell, there is a tight coupling between certain concepts associated with particles and the ``StaticSpace`` object.

Boundary Conditions
-------------------

.. doxygennamespace:: tw::bc
  :members:

Grid Axes
-------------------

.. doxygennamespace:: tw::grid
  :members:

Particle Primitives
-------------------

.. doxygenstruct:: Primitive
  :members:

.. doxygenstruct:: Particle
  :members:

.. doxygenstruct:: ParticleRef
  :members:

.. doxygenstruct:: TransferParticle
  :members:

.. doxygenstruct:: weights_3D
  :members:

StaticSpace Object
--------------------

.. doxygenstruct:: StaticSpace
  :members:
  :protected-members:
