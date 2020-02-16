NVIDIA GPGPU on RHEL/CentOS 8
=============================

.. Warning::

	The following is tested only for CentOS.

.. caution::

	We assume core turboWAVE **and** tools have already been installed according to the documentation, with no steps omitted.

Driver
------

#. Prepare for EPEL:

	* For CentOS, type ``sudo dnf config-manager --set-enabled PowerTools``
	* For RHEL, type ``ARCH=$( /bin/arch )`` followed by ``sudo subscription-manager repos --enable "codeready-builder-for-rhel-8-${ARCH}-rpms"``

#. Go to `EPEL <https://fedoraproject.org/wiki/EPEL>`_ and install.  As of this writing there is a link, ``epel-release-latest-8``, that allows you to use a GUI installer.
#. Go to `RPM Fusion <https://rpmfusion.org/Configuration>`_ and install the ``nonfree`` repository for RHEL 8 or compatible (there is no charge, ``nonfree`` refers to license restrictions).  The link should allow you to install via GUI.
#. Type ``sudo dnf install akmod-nvidia``

	* This automatic kernel module recompiles automatically when a new Linux kernel is installed (e.g. during a system update).  After restarting you must allow extra time for the kernel module to compile.  There could be a long delay before the login screen appears.

#. Restart the system, allow extra time for this restart.
#. ``sudo dnf install xorg-x11-drv-nvidia-cuda``
#. ``sudo dnf install opencl-headers``


Compile with OpenCL
--------------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = CUDA`, and :samp:`COMPILER_PREF = GNU` should be uncommented.
#. Open a new terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`make clean`
#. Type :samp:`make`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  The OpenCL kernel files will be copied into :samp:`~/Run`.  The OpenCL enabled code will not run without the kernel files.
