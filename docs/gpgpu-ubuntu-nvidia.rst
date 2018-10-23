NVIDIA GPGPU on Ubuntu 18.04
============================

.. caution::

	We assume core turboWAVE **and** tools have already been installed according to the documentation, with no steps omitted.

Driver
------

It is possible to install all the necessary packages using ``apt`` (no need to visit NVIDIA website).

	#. :samp:`sudo apt update`
	#. :samp:`sudo add-apt-repository ppa:graphics-drivers/ppa`
	#. :samp:`sudo apt install nvidia-driver-{XXX}`

		* Replace :samp:`{XXX}` with the version of your choice.  As of this writing the latest is 396.  Get a current list using :samp:`apt search nvidia-driver`.
		* As an alternative :samp:`sudo ubuntu-drivers autoinstall` is supposed to automatically select a suitable version.

	#. :samp:`sudo apt update`

Display Recovery
------------------

Installing graphics drivers in Linux can sometimes cause you to lose your display.  If this happens, try to switch to console mode by pressing :samp:`Ctrl-Alt-F2` (you may have to try different function keys).  If this succeeds you can issue the following commands to rollback the graphics driver:

	#. :samp:`sudo apt install ppa-purge`
	#. :samp:`ppa-purge ppa:graphics-drivers/ppa`
	#. Reboot using :samp:`sudo reboot`

Of course upon doing this GPU support may be lost.

Compile with OpenCL
-------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = CUDA`, and the compiler preference should be uncommented.
#. Edit :samp:`{twroot}/core/source/definitions.h`
#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`, and :samp:`#` is **not** a comment.  For this installation, only :samp:`#define USE_DESKTOP` and :samp:`#define USE_OPENCL` should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.
