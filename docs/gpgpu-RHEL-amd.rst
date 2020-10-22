AMD GPGPU on RHEL/CentOS 8
==========================

.. caution::

	We assume core turboWAVE **and** tools have already been installed according to the documentation, with no steps omitted.

Driver
-------

#. Install AMD ROCm

	* Perform internet search to find the installation instructions and carry out.  As of this writing the simplest way appears to be to use the script ``rocminstall.py``, see `<https://github.com/srinivamd/rocminstaller>`_.
	* Be sure to test the installation per the installation instructions.
	* This may involve multiple restarts.

#. Create a symbolic link to the ROCm installation

	* :samp:`cd /opt && ls`
	* The output should include a directory in the form :samp:`rocm-{x.y.z}`.
	* If there is no symbolic link :samp:`rocm` create it using :samp:`sudo ln -s rocm-{x.y.z} rocm`

Compile with OpenCL
-------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = RADEON_PRO`, and the compiler preference, should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. :samp:`make clean`
#. :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.
