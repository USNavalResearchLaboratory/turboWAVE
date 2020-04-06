GPGPU Acceleration
==================

General Notes
-------------

As of this writing, GPGPU support is only useful for Quantum Optics modules.
To offload computations to a GPGPU, turboWAVE uses OpenCL.  Typically OpenCL can be used for any modern GPU from either AMD or NVIDIA.  In the case of NVIDIA devices, the OpenCL support is packaged with CUDA.

.. warning::

	Enabling GPGPU acceleration involves manipulating display drivers, which on Linux can be treacherous.  If a Linux display manager is lost, it can be difficult to recover without a full OS reinstallation.  If you cannot afford for this to happen, you should take steps to backup your system configuration.

To run on GPGPU you must prepare a special executable.  The procedure for several operating systems is below.  Once the executable is prepared, running is almost the same as it is for a CPU.  The primary differences are as follows:

	#. Use a single MPI processes per GPGPU device, i.e., if your system has dual video cards you can use 2 MPI processes.
	#. OpenMP threads cannot be used.  The number of OpenMP threads must be one.
	#. If you want to control the particular OpenCL platform and device, use the command line arguments.  Otherwise turboWAVE will select the first available.

GPGPU Support by OS
-------------------

.. toctree::
	:maxdepth: 1

	gpgpu-mac
	gpgpu-win
	gpgpu-RHEL-amd
	gpgpu-RHEL-nvidia
	gpgpu-ubuntu-amd
	gpgpu-ubuntu-nvidia
	troubleshooting-ocl
