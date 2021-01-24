
.. toctree::
	:maxdepth: 2

Overview
========

turboWAVE Description
---------------------

TurboWAVE is a set of simulation modules built on a common framework.
The modules include PIC and cold fluid models, hydrodynamics models (SPARC),
nonlinear optics models, and quantum optics models.
All these modules can be invoked from the same executable file.
TurboWAVE consists of the following software components:

	#.	:samp:`twutils`---Python package with installer, visualization and other utilities.
	#.	:samp:`twcad`---CAD rendering program to visualize turboWAVE geometries
	#.	:samp:`DataViewer`---A Jupyter notebook for viewing and animating the data
	#.	:samp:`SPARC database`---chemical reactions for use in SPARC modules
	#.	:samp:`tw3d`---The executable file

On the desktop, turboWAVE is essentially self-contained, in that no external libraries are
needed, other than those commonly included with every C++ compiler.
TurboWAVE can take advantage of multi-core desktop systems using an internal implementation
of MPI which is transparent to the user.  TurboWAVE's internal MPI does not require any
external launch facility such as ``mpirun``, nor do any external MPI libraries need to be
installed.  However, to use distributed clusters, an external MPI library is needed.

Optional Components
-------------------

TurboWAVE can be run as a hybrid MPI-OpenMP code.
This means it can use a combination of distributed (MPI) and shared memory (OpenMP)
parallel processing models.  Unlike MPI, the OpenMP support relies on availability of
OpenMP libraries.  Fortunately, OpenMP is becoming a standard part of modern compilers.

TurboWAVE supports hardware acceleration with Graphical Processing Units (GPU) through OpenCL.
However, this support is limited to selected modules.
GPU acceleration requires that OpenCL libraries be installed, and that some hardware device
that supports OpenCL be available.
