GPGPU Acceleration
==================

General Notes
-------------

As of this writing, GPGPU support is only useful for Quantum Optics modules.
To offload computations to a GPGPU, turboWAVE uses OpenCL.  Typically OpenCL can be used for any modern GPU from either AMD or NVIDIA.  In the case of NVIDIA devices, the OpenCL support is packaged with CUDA.

.. warning::

	Enabling GPGPU acceleration involves manipulating display drivers, which on Linux can be treacherous.  If a Linux display manager is lost, it can be difficult to recover without a full OS reinstallation.  If you cannot afford for this to happen, you should take steps to backup your system configuration.

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
