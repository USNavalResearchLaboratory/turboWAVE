# makefile for turboWAVE
# D. Gordon, Naval Research Laboratory

# 'make' produces optimized executable and copies files
# 'make debug' produces debugging executable and copies files
# 'make tw3d' produces optimized executable (no copy files)
# 'make tw3d_debug' produces debugging executable (no copy files)
# 'make copy_files' copies files only
# For specific platforms and hardware accelerators comment/uncomment variable definitions below

all: tw3d copy_files
debug: tw3d_debug copy_files

#################################################
# BEGIN INPUT VARIABLES BLOCK

# Uncomment exactly one platform

#PLATFORM = OSX
PLATFORM = LINUX
#PLATFORM = CRAY
#PLATFORM = SGI

# Uncomment one, and at most one, of the following to use hardware acceleration
# N.b. specific accelerator devices can be specified in the turboWAVE input file

HARDWARE_ACCEL = OMP
#HARDWARE_ACCEL = APPLE_CL
#HARDWARE_ACCEL = CUDA
#HARDWARE_ACCEL = RADEON_PRO

# Specify a preference for compiler (ignored in some cases)

COMPILER_PREF = GNU
#COMPILER_PREF = LLVM_CLANG
#COMPILER_PREF = INTEL
#COMPILER_PREF = CRAY

# Specify a package manager preference (ignored in some cases)

#PACKAGE_PREF = HOMEBREW
PACKAGE_PREF = MACPORTS

# END INPUT VARIABLES BLOCK
#################################################

ifndef HARDWARE_ACCEL
	HARDWARE_ACCEL = NONE
endif

ifeq ($(PLATFORM),OSX)
	ifeq ($(HARDWARE_ACCEL),NONE)
		TW_Compiler = clang++
		TW_Linker = clang++
		CCFLAGS = -c -std=c++11 -O3
		DBFLAGS = -c -g -std=c++11 -O0
		LIB_PATH =
		LIBS = -lc++
	endif
	ifeq ($(HARDWARE_ACCEL),APPLE_CL)
		TW_Compiler = clang++
		TW_Linker = clang++
		CCFLAGS = -c -std=c++11 -O3
		DBFLAGS = -c -g -std=c++11 -O0
		LIB_PATH =
		LIBS = -lc++ -framework OpenCL
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		ifeq ($(PACKAGE_PREF),MACPORTS)
			# Tested for Mac OS 10.13
			# Remember to run port select clang mp-clang-x.x
			TW_Compiler = /opt/local/bin/clang++
			TW_Linker = /opt/local/bin/clang++
			CCFLAGS = -c -O3 -std=c++11 -march=native -fopenmp -I/opt/local/include
			#CCFLAGS = -c -O3 -std=c++11 -fopenmp -Rpass-analysis=loop-vectorize -I/opt/local/include
			DBFLAGS = -c -g -O0 -std=c++11 -fopenmp -I/opt/local/include
			LIB_PATH = -L/opt/local/lib
			LIBS = -lc++ -fopenmp
		else
			# DEFAULT TO HOMEBREW
			# Tested for Mac OS 10.13
			TW_Compiler = /usr/local/opt/llvm/bin/clang++
			TW_Linker = /usr/local/opt/llvm/bin/clang++
			CCFLAGS = -c -O3 -std=c++11 -march=native -fopenmp -I/usr/local/opt/llvm/include/
			#CCFLAGS = -c -O3 -std=c++11 -fopenmp -Rpass-analysis=loop-vectorize -I/usr/local/opt/llvm/include
			DBFLAGS = -c -g -O0 -std=c++11 -fopenmp -I/usr/local/opt/llvm/include/
			LIB_PATH = -L/usr/local/opt/llvm/lib
			LIBS = -lc++ -fopenmp
		endif
	endif
	BINARY_PATH = ~/bin
	WORK_PATH = ~/Run
	TW_MPI = yes
endif

ifeq ($(PLATFORM),LINUX)
	ifeq ($(COMPILER_PREF),LLVM_CLANG)
		TW_Compiler = clang++
		TW_Linker = clang++
		LIBS0 = -lc++
	else
		# DEFAULT TO GNU
		TW_Compiler = g++
		TW_Linker = g++
		LIBS0 = -lstdc++
	endif
	ifeq ($(HARDWARE_ACCEL),NONE)
		CCFLAGS = -c -std=c++11 -O3
		DBFLAGS = -c -g -std=c++11 -O0
		LIB_PATH =
		LIBS = $(LIBS0) -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		CCFLAGS = -c -std=c++11 -O3 -march=native -fopenmp# -Rpass-analysis=loop-vectorize
		DBFLAGS = -c -g -std=c++11 -O0 -fopenmp
		LIB_PATH =
		LIBS = $(LIBS0) -fopenmp
	endif
	ifeq ($(HARDWARE_ACCEL),CUDA)
		CCFLAGS = -c -std=c++11 -O3
		DBFLAGS = -c -g -std=c++11 -O0
		LIB_PATH =
		LIBS = $(LIBS0) -lOpenCL -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),RADEON_PRO)
		CCFLAGS = -c -std=c++11 -O3
		DBFLAGS = -c -g -std=c++11 -O0
		LIB_PATH =
		LIBS = $(LIBS0) -lOpenCL -pthread
	endif
	BINARY_PATH = ~/bin
	WORK_PATH = ~/Run
	TW_MPI = yes
endif

ifeq ($(PLATFORM),CRAY)
	ifeq ($(COMPILER_PREF),CRAY)
		ifeq ($(HARDWARE_ACCEL),NONE)
			TW_Compiler = CC
			TW_Linker = CC
			CCFLAGS = -c -O2 -h std=c++11 -h noomp
			DBFLAGS = -c -g -O0 -h std=c++11 -h noomp
			LIB_PATH =
			LIBS =
		endif
		ifeq ($(HARDWARE_ACCEL),OMP)
			TW_Compiler = CC
			TW_Linker = CC
			CCFLAGS = -c -O2 -h std=c++11
			#CCFLAGS = -c -O2 -hlist=m -h std=c++11
			DBFLAGS = -c -g -O0 -h std=c++11
			LIB_PATH =
			LIBS =
		endif
	else
		# DEFAULT TO INTEL
		ifeq ($(HARDWARE_ACCEL),NONE)
			TW_Compiler = CC
			TW_Linker = CC
			CCFLAGS = -c -std=c++11 -fast
			DBFLAGS = -c -g -std=c++11
			LIB_PATH =
			LIBS =
		endif
		ifeq ($(HARDWARE_ACCEL),OMP)
			TW_Compiler = CC
			TW_Linker = CC
			CCFLAGS = -c -std=c++11 -qopenmp -qopenmp-simd -fast
			#CCFLAGS = -c -std=c++11 -qopenmp -qopenmp-simd -qopt-report -qopt-report-phase=vec -fast
			DBFLAGS = -c -g -std=c++11 -qopenmp -qopenmp-simd
			LIB_PATH =
			LIBS = -qopenmp -qopenmp-simd
			#LIBS = -qopenmp -qopenmp-simd -qopt-report -qopt-report-phase=vec
		endif
	endif
endif

ifeq ($(PLATFORM),SGI)
	ifeq ($(HARDWARE_ACCEL),NONE)
		TW_Compiler = icpc
		TW_Linker = icpc
		CCFLAGS = -c -std=c++11 -O3 -finline-functions
		DBFLAGS = -c -g -std=c++11 -O0
		LIB_PATH =
		LIBS = -lmpi
	endif
	ifeq ($(HARDWARE_ACCEL),CUDA)
		CUDA_ROOT = /p/home/apps/cuda/7.5
		TW_Compiler = icpc
		TW_Linker = icpc
		CCFLAGS = -c -std=c++11 -O3 -finline-functions -I$(CUDA_ROOT)/include
		DBFLAGS = -c -g -std=c++11 -O0 -I$(CUDA_ROOT)/include
		LIB_PATH = -L$(CUDA_ROOT)/lib64
		LIBS = -lmpi -lOpenCL
	endif
endif

ifdef TW_MPI
BASE_HEADERS = definitions.h tw_mpi.h tasks.h ctools.h 3dmath.h metricSpace.h 3dfields.h region.h tw_iterator.h
SIM_HEADERS = $(BASE_HEADERS) functions.h physics.h numerics.h computeTool.h elliptic.h parabolic.h hyperbolic.h fft.h injection.h input.h diagnostics.h module.h sim.h
MODULE_HEADERS = fieldSolve.h electrostatic.h laserSolve.h fluid.h quantum.h particles.h solidState.h
tw3d: 3dfields.cl particles.cl elliptic.cl hyperbolic.cl quantum.cl Main.o TW_MPI.o Tasks.o MetricSpace.o 3DFields.o Region.o FFT.o Physics.o Numerics.o ComputeTool.o Parabolic.o Elliptic.o Hyperbolic.o Functions.o Grid.o Module.o FieldSolve.o Electrostatic.o LaserSolve.o Fluid.o Quantum.o SolidState.o Particles.o Pusher.o Diagnostics.o Injection.o Input.o
	$(TW_Linker) -o tw3d Main.o TW_MPI.o Tasks.o MetricSpace.o 3DFields.o Region.o FFT.o Physics.o Numerics.o ComputeTool.o Parabolic.o Elliptic.o Hyperbolic.o Functions.o Grid.o Module.o FieldSolve.o Electrostatic.o LaserSolve.o Fluid.o Quantum.o SolidState.o Particles.o Pusher.o Diagnostics.o Injection.o Input.o $(LIB_PATH) $(LIBS)
tw3d_debug: 3dfields.cl particles.cl elliptic.cl hyperbolic.cl quantum.cl Main.d TW_MPI.d Tasks.d MetricSpace.d 3DFields.d Region.d FFT.d Physics.d Numerics.d ComputeTool.d Parabolic.d Elliptic.d Hyperbolic.d Functions.d Grid.d Module.d FieldSolve.d Electrostatic.d LaserSolve.d Fluid.d Quantum.d SolidState.d Particles.d Pusher.d Diagnostics.d Injection.d Input.d
	$(TW_Linker) -o tw3d Main.o TW_MPI.o Tasks.o MetricSpace.o 3DFields.o Region.o FFT.o Physics.o Numerics.o ComputeTool.o Parabolic.o Elliptic.o Hyperbolic.o Functions.o Grid.o Module.o FieldSolve.o Electrostatic.o LaserSolve.o Fluid.o Quantum.o SolidState.o Particles.o Pusher.o Diagnostics.o Injection.o Input.o $(LIB_PATH) $(LIBS)
else
BASE_HEADERS = definitions.h tasks.h ctools.h 3dmath.h metricSpace.h 3dfields.h region.h tw_iterator.h
SIM_HEADERS = $(BASE_HEADERS) functions.h physics.h numerics.h computeTool.h elliptic.h parabolic.h hyperbolic.h fft.h injection.h input.h diagnostics.h module.h sim.h
MODULE_HEADERS = fieldSolve.h electrostatic.h laserSolve.h fluid.h quantum.h particles.h solidState.h
tw3d: 3dfields.cl particles.cl elliptic.cl hyperbolic.cl quantum.cl Main.o Tasks.o MetricSpace.o 3DFields.o Region.o FFT.o Physics.o Numerics.o ComputeTool.o Parabolic.o Elliptic.o Hyperbolic.o Functions.o Grid.o Module.o FieldSolve.o Electrostatic.o LaserSolve.o Fluid.o Quantum.o SolidState.o Particles.o Pusher.o Diagnostics.o Injection.o Input.o
	$(TW_Linker) -o tw3d Main.o Tasks.o MetricSpace.o 3DFields.o Region.o FFT.o Physics.o Numerics.o ComputeTool.o Parabolic.o Elliptic.o Hyperbolic.o Functions.o Grid.o Module.o FieldSolve.o Electrostatic.o LaserSolve.o Fluid.o Quantum.o SolidState.o Particles.o Pusher.o Diagnostics.o Injection.o Input.o $(LIB_PATH) $(LIBS)
tw3d_debug: 3dfields.cl particles.cl elliptic.cl hyperbolic.cl quantum.cl Main.d Tasks.d MetricSpace.d 3DFields.d Region.d FFT.d Physics.d Numerics.d ComputeTool.d Parabolic.d Elliptic.d Hyperbolic.d Functions.d Grid.d Module.d FieldSolve.d Electrostatic.d LaserSolve.d Fluid.d Quantum.d SolidState.d Particles.d Pusher.d Diagnostics.d Injection.d Input.d
	$(TW_Linker) -o tw3d Main.o Tasks.o MetricSpace.o 3DFields.o Region.o FFT.o Physics.o Numerics.o ComputeTool.o Parabolic.o Elliptic.o Hyperbolic.o Functions.o Grid.o Module.o FieldSolve.o Electrostatic.o LaserSolve.o Fluid.o Quantum.o SolidState.o Particles.o Pusher.o Diagnostics.o Injection.o Input.o $(LIB_PATH) $(LIBS)
endif

copy_files:
ifdef BINARY_PATH
	cp tw3d $(BINARY_PATH)
endif
ifdef WORK_PATH
	cp *.cl $(WORK_PATH)
endif

clean:
	rm -f *.o
	rm -f tw3d

Main.o: Main.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Main.cpp

ifdef TW_MPI
TW_MPI.o: TW_MPI.cpp definitions.h tw_mpi.h
	$(TW_Compiler) $(CCFLAGS) TW_MPI.cpp
endif

Tasks.o: Tasks.cpp definitions.h tasks.h
	$(TW_Compiler) $(CCFLAGS) Tasks.cpp

MetricSpace.o: MetricSpace.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) MetricSpace.cpp

3DFields.o: 3DFields.cpp $(BASE_HEADERS) fft.h
	$(TW_Compiler) $(CCFLAGS) 3DFields.cpp

Region.o: Region.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Region.cpp

FFT.o: FFT.cpp definitions.h fft.h
	$(TW_Compiler) $(CCFLAGS) FFT.cpp

Physics.o: Physics.cpp $(BASE_HEADERS) physics.h functions.h
	$(TW_Compiler) $(CCFLAGS) Physics.cpp

Numerics.o: Numerics.cpp $(BASE_HEADERS) numerics.h
	$(TW_Compiler) $(CCFLAGS) Numerics.cpp

ComputeTool.o: ComputeTool.cpp $(BASE_HEADERS) numerics.h computeTool.h
	$(TW_Compiler) $(CCFLAGS) ComputeTool.cpp

Elliptic.o: Elliptic.cpp $(BASE_HEADERS) numerics.h computeTool.h elliptic.h
	$(TW_Compiler) $(CCFLAGS) Elliptic.cpp

Parabolic.o: Parabolic.cpp $(BASE_HEADERS) numerics.h computeTool.h parabolic.h
	$(TW_Compiler) $(CCFLAGS) Parabolic.cpp

Hyperbolic.o: Hyperbolic.cpp $(BASE_HEADERS) numerics.h computeTool.h hyperbolic.h
	$(TW_Compiler) $(CCFLAGS) Hyperbolic.cpp

Functions.o: Functions.cpp definitions.h ctools.h functions.h
	$(TW_Compiler) $(CCFLAGS) Functions.cpp

Grid.o: Grid.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Grid.cpp

Module.o: Module.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Module.cpp

FieldSolve.o: FieldSolve.cpp $(SIM_HEADERS) fieldSolve.h
	$(TW_Compiler) $(CCFLAGS) FieldSolve.cpp

Electrostatic.o: Electrostatic.cpp $(SIM_HEADERS) fieldSolve.h electrostatic.h
	$(TW_Compiler) $(CCFLAGS) Electrostatic.cpp

LaserSolve.o: LaserSolve.cpp $(SIM_HEADERS) fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) LaserSolve.cpp

Fluid.o: Fluid.cpp $(SIM_HEADERS) fieldSolve.h fluid.h
	$(TW_Compiler) $(CCFLAGS) Fluid.cpp

Quantum.o: Quantum.cpp $(SIM_HEADERS) fieldSolve.h quantum.h
	$(TW_Compiler) $(CCFLAGS) Quantum.cpp

SolidState.o: SolidState.cpp $(SIM_HEADERS) solidState.h
	$(TW_Compiler) $(CCFLAGS) SolidState.cpp

Particles.o: Particles.cpp $(SIM_HEADERS) particles.h fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) Particles.cpp

Pusher.o: Pusher.cpp $(SIM_HEADERS) particles.h fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) Pusher.cpp

Diagnostics.o: Diagnostics.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Diagnostics.cpp

Injection.o: Injection.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Injection.cpp

Input.o: Input.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Input.cpp

################################################
#              DEBUG DEPENDENCIES              #
################################################

Main.d: Main.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(DBFLAGS) Main.cpp

ifdef TW_MPI
TW_MPI.d: TW_MPI.cpp definitions.h tw_mpi.h
	$(TW_Compiler) $(DBFLAGS) TW_MPI.cpp
endif

Tasks.d: Tasks.cpp definitions.h tasks.h
	$(TW_Compiler) $(DBFLAGS) Tasks.cpp

MetricSpace.d: MetricSpace.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(DBFLAGS) MetricSpace.cpp

3DFields.d: 3DFields.cpp $(BASE_HEADERS) fft.h
	$(TW_Compiler) $(DBFLAGS) 3DFields.cpp

Region.d: Region.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(DBFLAGS) Region.cpp

FFT.d: FFT.cpp definitions.h fft.h
	$(TW_Compiler) $(DBFLAGS) FFT.cpp

Physics.d: Physics.cpp $(BASE_HEADERS) physics.h functions.h
	$(TW_Compiler) $(DBFLAGS) Physics.cpp

Numerics.d: Numerics.cpp $(BASE_HEADERS) numerics.h
	$(TW_Compiler) $(DBFLAGS) Numerics.cpp

ComputeTool.d: ComputeTool.cpp $(BASE_HEADERS) numerics.h computeTool.h
	$(TW_Compiler) $(DBFLAGS) ComputeTool.cpp

Elliptic.d: Elliptic.cpp $(BASE_HEADERS) numerics.h computeTool.h elliptic.h
	$(TW_Compiler) $(DBFLAGS) Elliptic.cpp

Parabolic.d: Parabolic.cpp $(BASE_HEADERS) numerics.h computeTool.h parabolic.h
	$(TW_Compiler) $(DBFLAGS) Parabolic.cpp

Hyperbolic.d: Hyperbolic.cpp $(BASE_HEADERS) numerics.h computeTool.h hyperbolic.h
	$(TW_Compiler) $(DBFLAGS) Hyperbolic.cpp

Functions.d: Functions.cpp definitions.h ctools.h functions.h
	$(TW_Compiler) $(DBFLAGS) Functions.cpp

Grid.d: Grid.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(DBFLAGS) Grid.cpp

Module.d: Module.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(DBFLAGS) Module.cpp

FieldSolve.d: FieldSolve.cpp $(SIM_HEADERS) fieldSolve.h
	$(TW_Compiler) $(DBFLAGS) FieldSolve.cpp

Electrostatic.d: Electrostatic.cpp $(SIM_HEADERS) fieldSolve.h electrostatic.h
	$(TW_Compiler) $(DBFLAGS) Electrostatic.cpp

LaserSolve.d: LaserSolve.cpp $(SIM_HEADERS) fieldSolve.h laserSolve.h
	$(TW_Compiler) $(DBFLAGS) LaserSolve.cpp

Fluid.d: Fluid.cpp $(SIM_HEADERS) fieldSolve.h fluid.h
	$(TW_Compiler) $(DBFLAGS) Fluid.cpp

Quantum.d: Quantum.cpp $(SIM_HEADERS) fieldSolve.h quantum.h
	$(TW_Compiler) $(DBFLAGS) Quantum.cpp

SolidState.d: SolidState.cpp $(SIM_HEADERS) solidState.h
	$(TW_Compiler) $(DBFLAGS) SolidState.cpp

Particles.d: Particles.cpp $(SIM_HEADERS) particles.h fieldSolve.h laserSolve.h
	$(TW_Compiler) $(DBFLAGS) Particles.cpp

Pusher.d: Pusher.cpp $(SIM_HEADERS) particles.h fieldSolve.h laserSolve.h
	$(TW_Compiler) $(DBFLAGS) Pusher.cpp

Diagnostics.d: Diagnostics.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(DBFLAGS) Diagnostics.cpp

Injection.d: Injection.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(DBFLAGS) Injection.cpp

Input.d: Input.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(DBFLAGS) Input.cpp
