Core Install for RHEL/CentOS/SL 7.5
===================================

The GCC packaged with RHEL7 has insufficient OpenMP support for turboWAVE.  We will use the RHEL software collections mechanism to overcome this issue.

.. note::

	There is significant variability among RHEL OS installations due to the variability of options that can be selected at install time.  The below is tested for a KDE Plasma Workspaces installation.

.. note::

	CentOS and Scienfific Linux are clones of RHEL, and so should behave the same way.

Install Developer Software Collection
-------------------------------------

#. Open a terminal window
#. Point :samp:`yum` to needed repositories

	* For CentOS, :samp:`sudo yum install centos-release-scl`
	* For Scientific Linux, :samp:`sudo yum install yum-conf-softwarecollections`
	* For Red Hat, :samp:`sudo yum-config-manager --enable rhel-server-rhscl-7-rpms`

#. :samp:`sudo yum install devtoolset-7`

The developer toolset is installed in :samp:`/opt/rh/` and has to be specially enabled using the :samp:`scl` tool each time a component is invoked.  This is accounted for in the compiling step given below.

Configure File System
---------------------

#. Put the turboWAVE components (:samp:`core` and :samp:`tools`) into some directory, denoted :samp:`{twroot}`.
#. Open a terminal window
#. :samp:`cd ~`
#. :samp:`mkdir bin`
#. :samp:`mkdir Run`
#. Verify bash as the login shell

	* :samp:`echo $0` prints the currently running shell
	* :samp:`chsh -s /bin/bash` changes the default to bash

#. Edit :samp:`~/.bashrc`, adding the line :samp:`export PATH=~/bin/:$PATH`
#. When exporting variables in bash, do not put spaces around the equals sign.
#. Close all terminal windows.

Compile turboWAVE
-----------------

#. Edit :samp:`{twroot}/core/source/makefile`
#. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = LINUX`, :samp:`HARDWARE_ACCEL = OMP`, and :samp:`COMPILER_PREF = GNU`, should be uncommented.
#. Edit :samp:`{twroot}/core/source/definitions.h`
#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`.  For this installation, only :samp:`#define USE_DESKTOP` and :samp:`#define USE_OPENMP` should be uncommented.
#. Open a terminal window and navigate to :samp:`{twroot}/core/source`
#. Type :samp:`scl enable devtoolset-7 'make'`
#. The makefile should automatically copy the executable into your :samp:`~/bin` directory for later use.  OpenCL kernel files may also be copied into :samp:`~/Run`, but these will not be used.
#. If you use :samp:`vim`, you may want to copy :samp:`{twroot}/tools/config-files/filetype.vim` to :samp:`~/.vim`.  This will enable syntax highlighting while editing turboWAVE input files.

If you want to enable the developer toolset for the duration of a session, start a new shell using :samp:`scl enable devtoolset-7 'bash'`
