
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

The core executable uses only three libraries: `std`, MPI, and OpenMP.  Older versions could be
compiled with only `std`, but now MPI and OpenMP are required.  You can optionally compile
with OpenCL.
