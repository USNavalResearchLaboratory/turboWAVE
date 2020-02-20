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

BASE_HEADERS = definitions.h tasks.h ctools.h 3dmath.h tw_iterator.h discreteSpace.h metricSpace.h functions.h fft.h numerics.h 3dfields.h region.h input.h
TOOL_HEADERS = computeTool.h parabolic.h elliptic.h hyperbolic.h injection.h physics.h chemistry.h diagnostics.h
SIM_HEADERS = $(BASE_HEADERS) $(TOOL_HEADERS) module.h simulation.h
MODULE_HEADERS = fieldSolve.h electrostatic.h laserSolve.h fluid.h quantum.h particles.h solidState.h
OBJ_LIST = 3DFields.obj Chemistry.obj ComputeTool.obj Diagnostics.obj DiscreteSpace.obj Electrostatic.obj Elliptic.obj FFT.obj FieldSolve.obj Fluid.obj Functions.obj Hyperbolic.obj Injection.obj Input.obj LaserSolve.obj Main.obj MetricSpace.obj Module.obj Numerics.obj Parabolic.obj Particles.obj Physics.obj Pusher.obj Quantum.obj Region.obj Simulation.obj SolidState.obj Tasks.obj TW_MPI.obj

all: tw3d copy_files

tw3d: 3dfields.cl particles.cl elliptic.cl hyperbolic.cl quantum.cl $(OBJ_LIST)
	$(TW_Linker) $(LKFLAGS) $(OBJ_LIST) $(LIB_PATH) $(LIBS)

copy_files:
	copy tw3d.exe $(BINARY_PATH)
	copy *.cl $(WORK_PATH)

clean:
	del *.obj
	del tw3d.exe

3DFields.obj: 3DFields.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) 3DFields.cpp

Chemistry.obj: Chemistry.cpp $(BASE_HEADERS) computeTool.h physics.h chemistry.h
	$(TW_Compiler) $(CCFLAGS) Chemistry.cpp

ComputeTool.obj: ComputeTool.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) ComputeTool.cpp

Diagnostics.obj: Diagnostics.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Diagnostics.cpp

DiscreteSpace.obj: DiscreteSpace.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) DiscreteSpace.cpp

Electrostatic.obj: Electrostatic.cpp $(SIM_HEADERS) fieldSolve.h electrostatic.h
	$(TW_Compiler) $(CCFLAGS) Electrostatic.cpp

Elliptic.obj: Elliptic.cpp $(BASE_HEADERS) computeTool.h elliptic.h
	$(TW_Compiler) $(CCFLAGS) Elliptic.cpp

FFT.obj: FFT.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) FFT.cpp

FieldSolve.obj: FieldSolve.cpp $(SIM_HEADERS) fieldSolve.h
	$(TW_Compiler) $(CCFLAGS) FieldSolve.cpp

Fluid.obj: Fluid.cpp $(SIM_HEADERS) fieldSolve.h fluid.h
	$(TW_Compiler) $(CCFLAGS) Fluid.cpp

Functions.obj: Functions.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Functions.cpp

Hyperbolic.obj: Hyperbolic.cpp $(BASE_HEADERS) computeTool.h hyperbolic.h
	$(TW_Compiler) $(CCFLAGS) Hyperbolic.cpp

Injection.obj: Injection.cpp $(BASE_HEADERS) computeTool.h injection.h
	$(TW_Compiler) $(CCFLAGS) Injection.cpp

Input.obj: Input.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Input.cpp

LaserSolve.obj: LaserSolve.cpp $(SIM_HEADERS) fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) LaserSolve.cpp

Main.obj: Main.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Main.cpp

MetricSpace.obj: MetricSpace.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) MetricSpace.cpp

Module.obj: Module.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Module.cpp

Numerics.obj: Numerics.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Numerics.cpp

Parabolic.obj: Parabolic.cpp $(BASE_HEADERS) computeTool.h parabolic.h
	$(TW_Compiler) $(CCFLAGS) Parabolic.cpp

Particles.obj: Particles.cpp $(SIM_HEADERS) particles.h fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) Particles.cpp

Physics.obj: Physics.cpp $(BASE_HEADERS) computeTool.h physics.h
	$(TW_Compiler) $(CCFLAGS) Physics.cpp

Pusher.obj: Pusher.cpp $(SIM_HEADERS) particles.h fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) Pusher.cpp

Quantum.obj: Quantum.cpp $(SIM_HEADERS) fieldSolve.h quantum.h
	$(TW_Compiler) $(CCFLAGS) Quantum.cpp

Region.obj: Region.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Region.cpp

Simulation.obj: Simulation.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Simulation.cpp

SolidState.obj: SolidState.cpp $(SIM_HEADERS) solidState.h
	$(TW_Compiler) $(CCFLAGS) SolidState.cpp

Tasks.obj: Tasks.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Tasks.cpp

TW_MPI.obj: TW_MPI.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) TW_MPI.cpp
