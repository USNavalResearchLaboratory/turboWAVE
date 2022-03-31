Input File: Geometry
====================

Basic geometric shapes are created using a ``region`` block.  Named regions can be embedded inside a compound region to perform boolean operations on the volumes.  This can be done recursively to build up complex geometries.

The geometry can be viewed using ``twcad``, a Python program that reads a turboWAVE input file and displays the geometry in a CAD style window.

An important distinction in creating geometries is whether an object is defined in parameter space or real space.

.. glossary::

	Parameter Space
		Objects in parameter space are defined in a specific coordinate system, and become a different object when the coordinate system is changed.  For example, the circle is defined in Cartesian coordinates as :math:`x^2+z^2=1`.  If cylindrical coordinates are used, the equation retains its form, :math:`\varrho^2+z^2=1`, which implies the geometry of the object changes.  This can sometimes be used to advantage, e.g., one can create a torus in real space using a circle in parameter space, :math:`(\varrho-5)^2+z^2=1`.

	Real Space
		Objects in real space are the same object regardless of coordinates.  If a sphere is created in real space, it remains a sphere no matter what the grid geometry is.

At present most objects are created in parameter space.  This distinction is only important on a curvilinear grid.

When defining geometry using the input file, 3-tuples are often used to specify spatial coordinates.  The meaning of the tuple components depends on the coordinate system, as shown in Table I.

.. csv-table:: Table I. Input File Tuple Coordinate Mappings.
	:header: "System", "Tuple Order", "Comment"
	:delim: ;

	"Cartesian";  :math:`(x,y,z)`; :math:`{\bf e}_x\times{\bf e}_y = {\bf e}_z`
	"Cylindrical";  :math:`(\varrho,\varphi,z)`; :math:`\varrho^2 = x^2 + y^2`
	"Spherical";  :math:`(r,\varphi,\theta)`; "Polar angle is last"


TurboWAVE CAD Viewer
--------------------

There is a Python program for viewing turboWAVE geometries in :file:`{turboWAVE}/tools/twcad/`.  If the program is run with ``stdin`` in the working directory, a 3D viewing window is created allowing the user to inspect the geometry from any angle.  This program has to be run from a special Python environment.  For more on how to install and run this program see the documentation in :file:`{turboWAVE}/tools/twcad/`.

.. _shared-geometry:

Geometry Shared Directives
--------------------------

The following directives may be used with any type of region.

.. py:function:: translation = ( x , y , z )

	Displace region from the current position.  Each type of region has a standard initial position, generally corresponding to being centered at the origin.  This directive can be repeated to create a sequence of operations.

	:param float x: translation in the x direction
	:param float y: translation in the y direction
	:param float z: translation in the z direction

.. py:function:: rotation about axis = angle

	Rotate about the given global axis. This directive can be repeated to create a sequence of operations. Rotations about a shifted axis can be carried out by performing the sequence :math:`TRT^{-1}` with :math:`T` a translation and :math:`R` a rotation.

	:param enum axis: can be ``x``, ``y``, or ``z``
	:param float angle: the rotation angle in radians.  You can use degrees by adding a dimension, e.g., ``rotation about axis = 45 [deg]``.

.. py:function:: complement = tst

	:param bool tst: if true, transforms the region into its complement

.. py:function:: move with window = tst

	:param bool tst: if true, the region moves with the window


Basic Region Types
------------------

.. py:function:: new region rect <name> { <directives> }

	Creates a rectangular box.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: bounds = ( x0 , x1 , y0 , y1 , z0 , z1 )

			Defines the limits of the rectangular box. This implies a translation transformation, so **ordering with respect to other transformations matters**.

.. py:function:: new region circ <name> { <directives> }

	Creates a circle or sphere in parameter space.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: radius = R

			:param float R: defines the radius of the circle

.. py:function:: new region true_sphere <name> { <directives> }

	Creates a sphere in real space.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: radius = R

			:param float R: defines the radius of the sphere

.. py:function:: new region prism <name> { <directives> }

	Creates a prism in parameter space.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: bounds = ( x0 , x1 , y0 , y1 , z0 , z1 )

			Defines the limits of the bounding rectangular box. This implies a translation transformation, so **ordering with respect to other transformations matters**.
			The tip points in the +x direction.

.. py:function:: new region ellipsoid <name> { <directives> }

	Creates an ellipsoid in parameter space.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: bounds = ( x0 , x1 , y0 , y1 , z0 , z1 )

			Defines the limits of the bounding rectangular box. This implies a translation transformation, so **ordering with respect to other transformations matters**.

.. py:function:: new region cylinder <name> { <directives> }

	Creates a cylinder in parameter space. Default orientation is centered on the z-axis.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: radius = R

			:param float R: radius of the cylinder

		.. py:function:: length = L

			:param float L: length of cylinder


.. py:function:: new region rounded_cylinder <name> { <directives> }

	Creates a cylinder in parameter space, with hemispherical end-caps. Default orientation is centered on the z-axis.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: radius = R

			:param float R: radius of the cylinder

		.. py:function:: length = L

			:param float L: length of cylinder, does not count the end-caps


.. py:function:: new region cylindrical_shell <name> { <directives> }

	Creates a cylindrical shell, or "tube", in parameter space. Default orientation is centered on the z-axis.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: inner radius = R1

			:param float R1: inner radius of the tube

		.. py:function:: outer radius = R2

			:param float R2: outer radius of the tube

		.. py:function:: length = L

			:param float L: length of tube


.. py:function:: new region torus <name> { <directives> }

	Creates a torus in parameter space. Default orientation is centered on the z-axis.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: minor radius = R1

			:param float R1: radius of tube that is bent to make the torus

		.. py:function:: major radius = R2

			:param float R2: distance from the torus center to the tube center


.. py:function:: new region cone <name> { <directives> }

	Creates a cone in parameter space. Default orientation is centered on the z-axis, i.e., the origin is mid-way between the base and the tip.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: tip radius = R1

			:param float R1: radius of blunted tip

		.. py:function:: base radius = R2

			:param float R2: radius of the cone base

		.. py:function:: length = L

			:param float L: distance between base and tip

.. py:function:: new region tangent_ogive <name> { <directives> }

	Creates a spherically blunted tangent ogive in parameter space. Default orientation is centered on the z-axis, with the tip on the +z side.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: tip radius = R1

			:param float R1: radius of blunted tip

		.. py:function:: base radius = R2

			:param float R2: radius of the cone base

		.. py:function:: length = L

			:param float L: distance between base and spherical tip


.. py:function:: new region box_array <name> { <directives> }

	Creates an infinite array of box shaped regions.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: size = ( dx , dy , dz )

			size of each box

		.. py:function:: spacing = ( lx , ly , lz )

			center-to-center distance between boxes


Compound Regions
----------------

.. py:function:: new region union <name> { <directives> }

	Create the union of several other regions.  This is the *boolean or*, i.e., if the union has elements A, B, and C, then a point in the union must be in A or in B or in C.  If C is the complement of D, then the point must be in A or in B or not in D.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: elements = { R1 , R2 , R3 , ... }

			variable length list of names of regions forming the union

.. py:function:: new region intersection <name> { <directives> }

	Create the intersection of several other regions.  This is the *boolean and*, i.e., if the intersection has elements A, B, and C, then a point in the intersection must be in A and in B and in C.  If C is the complement of D, then the point must be in A and in B and not in D.

	:param str name: assigns a name to the region
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`shared-geometry`

		.. py:function:: elements = { R1 , R2 , R3 , ... }

			variable length list of names of regions forming the intersection

.. tip::

	You can create what most CAD software calls a difference using intersections.  The difference of A and B is the intersection of A and the complement of B.

Specific Example in 2D
----------------------

In this example we make a square with a rounded top in the x-z plane, with a hole in it.
First define the elements of the compound region::

	new region rect r1
	{
		bounds =  -1.0 1.0 -1.0 1.0 -1.0 1.0
	}

	new region circ c1
	{
		translation = 0.0 0.0 1.0
		radius = 1.0
	}

	new region circ !c2 // exclamation point just reminds us this is the complement of a circle
	{
		radius = 0.5
		complement = true
	}

Now make the square with rounded top::

	new region union u1
	{
		elements = { r1 , c1 }
	}

Now put a hole in it by using the union in an intersection::

	new region intersection i1
	{
		elements = { u1 , !c2 }
	}

The named region "i1" can now be used as the clipping region in :ref:`matter-loading`, in certain diagnostics, or in :ref:`conductor`.
