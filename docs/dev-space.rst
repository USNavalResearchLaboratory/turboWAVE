Static Space
==============

The ``StaticSpace`` object is the basis of structured grid management in turboWAVE. While ``StaticSpace`` provides topology and indexing, ``MetricSpace`` adds edge lengths, wall areas, and cell volumes, which define a specific grid geometry.

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

DynSpace Object
--------------------

.. doxygenstruct:: StaticSpace
  :members:
  :protected-members:
