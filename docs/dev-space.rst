Discrete Space
==============

The ``DiscreteSpace`` object is the basis of structured grid management in turboWAVE. While ``DiscreteSpace`` provides topology and indexing, ``MetricSpace`` adds edge lengths, wall areas, and cell volumes, which define a specific grid geometry.

There is a tight coupling between certain concepts associated with particles and the ``DiscreteSpace`` object.  As a result some of the lowest level particle objects appear in this scope.

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

DiscreteSpace Object
--------------------

.. doxygenstruct:: DiscreteSpace
  :members:
  :protected-members:
