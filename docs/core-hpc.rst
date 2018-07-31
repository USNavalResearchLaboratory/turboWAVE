Core Install for HPC
====================

In referring to High Performance Computing (HPC), we have in mind a large scale computing cluster, managed by a professional staff, and accessed remotely.

We will focus on Cray systems in this instruction.  For other vendors the below will often go through with little modification.  At present almost all systems use little-endian binary numbers.  If your system is big-endian, you must change the value of a boolean variable in :samp:`definitions.h`.

HPC modules
-----------

HPC systems usually have a system that allows for easy loading, unloading, and swapping of modules.  The commands include :samp:`module load {module}`, :samp:`module unload {module}`, and :samp:`module swap {oldModule} {newModule}`.
For Cray systems, we usually want the Intel compiler.  This would be loaded with :samp:`module load PrgEnv-intel`.  If another compiler is already loaded you can use, say, :samp:`module swap PrgEnv-cray PrgEnv-intel`.
The NERSC supercomputing center automatically loads the darshan module, which doesn't play with turboWAVE.  Before compiling run :samp:`module unload darshan`.

Compiling on Cray Systems
-------------------------

  #. Make a directory on the HPC system for turboWAVE source.  We denote it :samp:`{turbowave}`

  #. Copy everything in the :samp:`{twroot}/core/source/` directory to :samp:`{turbowave}`.

  #. Navigate to :samp:`{turbowave}`

  #. Edit :samp:`{turbowave}/makefile`

  #. In the makefile, you must comment/uncomment lines to select platform, hardware acceleration, compiler, and package manager.  You will only be editing the lines between :samp:`BEGIN INPUT VARIABLES BLOCK` and :samp:`END INPUT VARIABLES BLOCK`.  In a makefile, comments are preceded by :samp:`#`.  For this installation, only :samp:`PLATFORM = CRAY`, :samp:`HARDWARE_ACCEL = OMP`, and :samp:`COMPILER_PREF = INTEL`, should be uncommented.

  #. Edit :samp:`{turbowave}/definitions.h`

  #. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`.  For this installation, only :samp:`#define USE_CRAY` and :samp:`#define USE_OPENMP` should be uncommented.

  #. Type :samp:`make`

  #. You must manually copy the executable to the scratch directory.  For example, at NERSC, this would be done with :samp:`cp tw3d $SCRATCH`.
