# Windows NMAKE makefile for turboWAVE
# D. Gordon, Naval Research Laboratory

# Invoke from developer/intel command prompt using nmake /F win.make
# Switch between compilers by commenting/uncommenting COMPILER_PREF defs.
# Switch between release, debug, and profile by commenting/uncommenting CCFLAGS defs.
# (nmake seems to be missing target specific variable values so we must resort to this)

#COMPILER_PREF = VS
COMPILER_PREF = INTEL

!IF "$(COMPILER_PREF)"=="VS"
TW_Compiler = cl
TW_Linker = cl
RELEASE_FLAGS = /EHsc /c /O2 /DUSE_DESKTOP /DNDEBUG
PROFILE_FLAGS = /EHsc /c /O2 /Zi /DUSE_DESKTOP
DEBUG_FLAGS = /EHsc /c /Od /Zi /DUSE_DESKTOP
LIB_PATH =
LIBS =
!ENDIF

!IF "$(COMPILER_PREF)"=="INTEL"
TW_Compiler = icl
TW_Linker = icl
RELEASE_FLAGS = /EHsc /c /O2 /Qopenmp /QxHOST /DUSE_DESKTOP /DUSE_OPENMP /DNDEBUG
PROFILE_FLAGS = /EHsc /c /O2 /Qopenmp /QxHOST /Zi /DUSE_DESKTOP /DUSE_OPENMP
DEBUG_FLAGS = /EHsc /c /Od /Qopenmp /Zi /DUSE_DESKTOP /DUSE_OPENMP
LIB_PATH =
LIBS =
!ENDIF

CCFLAGS = $(RELEASE_FLAGS)
#CCFLAGS = $(PROFILE_FLAGS)
#CCFLAGS = $(DEBUG_FLAGS)

BINARY_PATH = $(HOMEPATH)\bin
WORK_PATH = $(HOMEPATH)\Run

!IF "$(CCFLAGS)"=="$(RELEASE_FLAGS)"
LKFLAGS = /Fetw3d
!ENDIF

!IF "$(CCFLAGS)"=="$(PROFILE_FLAGS)"
LKFLAGS = /Fetw3d /Zi
!ENDIF

!IF "$(CCFLAGS)"=="$(DEBUG_FLAGS)"
LKFLAGS = /Fetw3d /Zi
!ENDIF

BASE_HEADERS = definitions.h tw_mpi.h tasks.h ctools.h 3dmath.h metricSpace.h 3dfields.h region.h tw_iterator.h
SIM_HEADERS = $(BASE_HEADERS) functions.h physics.h chemistry.h numerics.h computeTool.h elliptic.h parabolic.h hyperbolic.h fft.h injection.h input.h diagnostics.h module.h simulation.h
MODULE_HEADERS = fieldSolve.h electrostatic.h laserSolve.h fluid.h quantum.h particles.h solidState.h

all: tw3d copy_files

tw3d: 3dfields.cl particles.cl elliptic.cl hyperbolic.cl quantum.cl Main.obj TW_MPI.obj Tasks.obj MetricSpace.obj 3DFields.obj Region.obj FFT.obj Physics.obj Chemistry.obj Numerics.obj ComputeTool.obj Parabolic.obj Elliptic.obj Hyperbolic.obj Functions.obj Simulation.obj Module.obj FieldSolve.obj Electrostatic.obj LaserSolve.obj Fluid.obj Quantum.obj SolidState.obj Particles.obj Pusher.obj Diagnostics.obj Injection.obj Input.obj
	$(TW_Linker) $(LKFLAGS) Main.obj TW_MPI.obj Tasks.obj MetricSpace.obj 3DFields.obj Region.obj FFT.obj Physics.obj Chemistry.obj Numerics.obj ComputeTool.obj Parabolic.obj Elliptic.obj Hyperbolic.obj Functions.obj Simulation.obj Module.obj FieldSolve.obj Electrostatic.obj LaserSolve.obj Fluid.obj Quantum.obj SolidState.obj Particles.obj Pusher.obj Diagnostics.obj Injection.obj Input.obj $(LIB_PATH) $(LIBS)

copy_files:
	copy tw3d.exe $(BINARY_PATH)
	copy *.cl $(WORK_PATH)

clean:
	del *.obj
	del tw3d.exe

Main.obj: Main.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Main.cpp

TW_MPI.obj: TW_MPI.cpp definitions.h tw_mpi.h
	$(TW_Compiler) $(CCFLAGS) TW_MPI.cpp

Tasks.obj: Tasks.cpp definitions.h tasks.h
	$(TW_Compiler) $(CCFLAGS) Tasks.cpp

MetricSpace.obj: MetricSpace.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) MetricSpace.cpp

3DFields.obj: 3DFields.cpp $(BASE_HEADERS) fft.h
	$(TW_Compiler) $(CCFLAGS) 3DFields.cpp

Region.obj: Region.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Region.cpp

FFT.obj: FFT.cpp definitions.h fft.h
	$(TW_Compiler) $(CCFLAGS) FFT.cpp

Physics.obj: Physics.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Physics.cpp

Chemistry.obj: Chemistry.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Chemistry.cpp

Numerics.obj: Numerics.cpp $(BASE_HEADERS) numerics.h
	$(TW_Compiler) $(CCFLAGS) Numerics.cpp

ComputeTool.obj: ComputeTool.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) ComputeTool.cpp

Elliptic.obj: Elliptic.cpp $(BASE_HEADERS) numerics.h computeTool.h elliptic.h
	$(TW_Compiler) $(CCFLAGS) Elliptic.cpp

Parabolic.obj: Parabolic.cpp $(BASE_HEADERS) numerics.h computeTool.h parabolic.h
	$(TW_Compiler) $(CCFLAGS) Parabolic.cpp

Hyperbolic.obj: Hyperbolic.cpp $(BASE_HEADERS) numerics.h computeTool.h hyperbolic.h
	$(TW_Compiler) $(CCFLAGS) Hyperbolic.cpp

Functions.obj: Functions.cpp definitions.h ctools.h functions.h
	$(TW_Compiler) $(CCFLAGS) Functions.cpp

Simulation.obj: Simulation.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Simulation.cpp

Module.obj: Module.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Module.cpp

FieldSolve.obj: FieldSolve.cpp $(SIM_HEADERS) fieldSolve.h
	$(TW_Compiler) $(CCFLAGS) FieldSolve.cpp

Electrostatic.obj: Electrostatic.cpp $(SIM_HEADERS) fieldSolve.h electrostatic.h
	$(TW_Compiler) $(CCFLAGS) Electrostatic.cpp

LaserSolve.obj: LaserSolve.cpp $(SIM_HEADERS) fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) LaserSolve.cpp

Fluid.obj: Fluid.cpp $(SIM_HEADERS) fieldSolve.h fluid.h
	$(TW_Compiler) $(CCFLAGS) Fluid.cpp

Quantum.obj: Quantum.cpp $(SIM_HEADERS) fieldSolve.h quantum.h
	$(TW_Compiler) $(CCFLAGS) Quantum.cpp

SolidState.obj: SolidState.cpp $(SIM_HEADERS) solidState.h
	$(TW_Compiler) $(CCFLAGS) SolidState.cpp

Particles.obj: Particles.cpp $(SIM_HEADERS) particles.h fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) Particles.cpp

Pusher.obj: Pusher.cpp $(SIM_HEADERS) particles.h fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) Pusher.cpp

Diagnostics.obj: Diagnostics.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Diagnostics.cpp

Injection.obj: Injection.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Injection.cpp

Input.obj: Input.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Input.cpp
