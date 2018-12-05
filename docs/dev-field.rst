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

Index based loops are convenient when complex ranges or differencing patterns are needed.  For commonly occurring patterns use range based loops instead.

Looping over Cells
------------------

Range based loops over cells are used when the order of iterations does not require any special control.  A typical construction is:

.. code-block:: c++

	#pragma omp parallel
	{
		for (auto cell : InteriorCellRange(*this))
			some_field(cell,0) *= 2;
	}

In this example the first component of ``some_field`` is doubled.  Putting the loop inside a parallel region causes the ``InteriorCellRange`` object to automatically partition the iterations among the threads. The use of ``*this`` in the range constructor assumes that the loop is formed in the method of some object that derives from ``DiscreteSpace``, such as a ``Module``.  If calling from a ``ComputeTool`` method, ``*this`` would be replaced by ``*space``.  The ``InteriorCellRange`` is so named because it does not include ghost cells. To include them use ``EntireCellRange``.

.. Caution::

	If for some reason you nest auto-threaded loops inside a parallel section, you will create a new team of threads for each thread in the outer team.

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

Here, the iterator method ``global_count`` is used to get the global index of the iteration, which is unique across threads.  The explicit example brings out the three elements of iterating through a ``Field``: the range (specific type ``CellRange``), the iterator (automatically typed variable ``it``), and the reference (specific type ``tw::cell``).  The ``CellRange`` range is the generalization of ``InteriorCellRange`` and ``EntireCellRange``.  The boolean argument chooses whether to include ghost cells.

.. Note::

	More elaborate ghost cell inclusion patterns are intended for future development.

Looping over Strips
-------------------

A frequent pattern is operating on strips of cells.  Often one would like to repeat the same strip-wise operations along each axis. Strip ranges make this simple.

.. code-block:: c++

	for (int ax=1;ax<=3;ax++)
	{
		#pragma omp parallel
		{
			for (auto strip : StripRange(*this,ax,strongbool::no))
				for (int s=0;s<=Dim(ax);s++)
					some_field(strip,s,0) *= 2.0;
		}
	}

The ``StripRange`` takes a new argument, an integer giving the axis parallel to the strips.  To avoid errors in the order of arguments, we require the strongly typed ``strongbool`` to indicate ghost cell inclusion.

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

Differencing
------------

The ``Field`` class provides for differencing patterns that occur often in computational physics.  For example:

.. code-block:: c++

	#pragma omp parallel
	{
		for (auto v : VectorizingRange<3>(*this,false))
		{
			#pragma omp simd
			for (int i=0;i<Dim(3);i++)
				A(v,i,0) = B.d2(v,i,0,2);
		}
	}

In mathematical notation this would be:

	:math:`A_0(x_1,x_2,x_3) = \frac{\partial^2}{\partial x_2^2}B_0(x_1,x_2,x_3)`

.. Note::

	When applying differencing operators the range must not include ghost cells.
