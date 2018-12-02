Field
=====

The ``Field`` object is essentially a 4-dimensional array, where one of the dimensions is a component, and the other three are spatial axes.  It provides the following specialized structure useful for simulation (not exhaustive):

	#. Common communication patterns between distributed compute nodes
	#. Iterators to simplify common loop ranges
	#. Flexible storage: storage pattern and index space are independent
	#. OpenMP aware: design accounts for multithreading and vectorization
	#. Finite volumes: awareness of structured grid metrics
	#. Slice operations

Index Space
-----------

No matter what the storage pattern is, the order of array indices always has the same meaning.  The first three indices are spatial, and the fourth indexes some set of elements that are known at a particular spacetime point, such as vector components.  The meaning of the spatial indices depends on the coordinate system chosen.  The ordering for the various coordinate systems is as follows.

.. csv-table:: Table I. Coordinates.
	:header: "System", "Ordering", "Comment"
	:delim: ;

	"Cartesian"; :math:`x,y,z`; :math:`{\bf e}_x\times{\bf e}_y = {\bf e}_z`
	"Cylindrical"; :math:`\varrho,\varphi,z`
	"Spherical"; :math:`r,\varphi,\theta`; "Polar angle is last"

Loopless Operations
-------------------

For operations acting uniformly on a field, there is no need for an explicit loop.  For example, to set all the components in a field to unity::

	some_field = 1.0;

Most C++ operators are supported.  For example, to divide the whole field by 3::

	some_field /= 3.0;

Index Based Loops
-----------------

You can loop through a field in legacy fashion:

.. code-block:: c++

	for (int i=0;i<5;i++)
		for (int j=0;j<5;j++)
			for (int k=0;k<5;k++)
				some_field(i,j,k,0) = 1.0;

This sets the first component to unity for the range of indices indicated.  Notice the FORTRAN style array access.

Index based loops are convenient when an unusual loop range is needed.  For commonly occurring ranges use range based loops instead.

Looping over Cells
------------------

Range based loops over cells are used when the order of iterations does not require any special control.  A typical construction is:

.. code-block:: c++

	#pragma omp parallel
	{
		for (auto cell : CellRange(*this,false))
			some_field(cell,0) *= 2;
	}

In this example the first component of ``some_field`` is doubled.  Putting the loop inside a parallel region causes the ``CellRange`` object to automatically partition the iterations among the threads. The use of ``*this`` in the range constructor assumes that the loop is formed in the method of some object that derives from ``DiscreteSpace``, such as a ``Module``.  If calling from a ``ComputeTool`` method, ``*this`` would be replaced by ``*space``.  The boolean argument determines whether ghost cells should be included in the iteration.

If an iteration index is required, one may use the more explicit notation:

.. code-block:: c++

	#pragma omp parallel
	{
		CellRange range(*this,false);
		for (auto it=range.begin();it!=range.end();++it)
		{
			tw::cell cell = *it;
			some_field(cell,1) = it.global_count();
		}
	}

Here, the iterator method ``global_count`` is used to get the global index of the iteration, which is unique across threads.  The explicit example brings out the three elements of iterating through a ``Field``: the range (specific type ``CellRange``), the iterator (automatically typed variable ``it``), and the reference (specific type ``tw::cell``).

Vectorization
-------------

In order to promote compiler vectorization, one has to commit to a particular storage pattern.  Special templated ranges and references must be used.  The template argument is an integer identifying the packed axis.  Once this type of construction is used, the storage pattern cannot be changed, unless all the code that makes use of vectorizing iterators is modified.

Suppose we have a ``Field`` with axis 3 as the packed axis.  Then an optimized loop might be constructed as follows:

.. code-block:: c++

	#pragma omp parallel
	{
		for (auto v : VectorizingRange<3>(*this,false))
		{
			#pragma omp simd
			for (tw::Int i=0;i<=Dim(3);i++)
				some_field(v,i,0) *= 2;
		}
	}

Here, we have again assumed the block is defined inside a derivative of ``DiscreteSpace``.  It is important to understand that this construction uses thread parallelism *across* strips, and vector parallelism *along* strips.  Therefore it is not effective for 1D problems.
