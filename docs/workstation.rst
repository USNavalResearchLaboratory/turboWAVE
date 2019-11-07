
.. toctree::
	:maxdepth: 2

.. highlight:: none

Workstation Setup
=================

Operating System
----------------

TurboWAVE can be readily installed on Linux, MacOS, or Windows.  The workflow tends to be UNIX-oriented, which is natural for Linux and MacOS, and thanks to PowerShell, fairly seamless even on Windows.

MacOS and Windows have native applications to visualize turboWAVE output files.
The native viewer for Windows is dated, however, and is not being updated.
For Linux, visualization is accomplished using either command line Python scripts,
or using interactive Python notebooks.

Desktop OS's that are explicitly supported at the time of this writing include MacOS 10.15, RHEL 8.0 and clones, Ubuntu 18.04, and Windows 10.  Supercomputer platforms that are known to work well include Cray and SGI.

Hardware
--------

At the time of this writing, the choice of CPU and memory capacity are the two major
factors affecting turboWAVE performance, while choice of GPU is not very important.
TurboWAVE is designed to effectively untilize multi-core processors.  The number of
cores is the most important performance hardware characteristic for turboWAVE. The more cores the better.
The memory requirements depend most critically on the dimensions of the problem.
A rule of thumb is to allow 1000 bytes of RAM for each grid cell.  For a serious 2D problem
this leads to about 1 GB.  For a small 3D problem, several tens of GB are required.
