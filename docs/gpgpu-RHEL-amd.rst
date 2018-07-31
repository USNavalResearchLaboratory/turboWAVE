AMD GPGPU on RHEL 7.3
=====================

.. caution::

	We assume core turboWAVE **and** tools have already been installed according to the documentation, with no steps omitted.

Driver
-------

Graphics drivers can change rapidly, so internet searches may figure prominently into your installation effort.
As of this writing the driver to use for AMD is the radeon-pro driver.  This typically supports AMD workstation class GPU and AMD or Intel CPU.  However, AMD gaming class GPU or APU may not be supported.

#. Find and download the radeon-pro software for RHEL 7.3 (use internet search)
#. Navigate to downloaded archive
#. :samp:`tar -Jxvf {downloaded_file}`
#. Navigate into unpacked directory
#. :samp:`sudo ./amd-gpu-install`
#. Reboot the system
#. The radeon-pro files should be in :samp:`/opt/amdgpu-pro/`.  To fill in the paths below you may need to dig around in this directory.  For the headers, you are looking for the :samp:`CL` directory.  For the libraries, you are looking for :samp:`libOpenCL.so`.
#. Edit :samp:`~/.bashrc`

	- Add line :samp:`export PATH=/opt/amdgpu-pro/bin:$PATH`
	- Add line :samp:`export CPATH=/opt/amdgpu-pro/{subpath_to_headers}:$CPATH`
	- Add line :samp:`export LIBRARY_PATH=/opt/amdgpu-pro/{subpath_to_libraries}:$LIBRARY_PATH`
	- Add line :samp:`export LD_LIBRARY_PATH=/opt/amdgpu-pro/{subpath_to_libraries}:$LIBRARY_PATH`
	- Remember, no spaces around equals sign.

#. Close all terminal windows, open a new one.
#. :samp:`clinfo`
#. If the above command gives a device listing which includes your CPU and your GPU the installation likely succeeded.

Compile with OpenCL
-------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = RADEON_PRO`, and the compiler preference, should be uncommented.
#. Edit :samp:`{twroot}/core/source/definitions.h`
#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`, and :samp:`#` is **not** a comment.  For this installation, only :samp:`#define USE_DESKTOP` and :samp:`#define USE_OPENCL` should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`scl enable devtoolset-7 'make'`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.
