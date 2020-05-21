# makefile for turboWAVE
# D. Gordon, Naval Research Laboratory

# 'make' produces optimized executable and copies files
# 'make debug' produces debugging executable and copies files
# 'make profile' produces profiling executable and copies files
# 'make tw3d_release' produces optimized executable (no copy files)
# 'make tw3d_debug' produces debugging executable (no copy files)
# 'make tw3d_profile' produces profiling executable (no copy files)
# 'make copy_files' copies files only
# 'make clean' removes all object files and the tw3d executable
# For specific platforms and hardware accelerators comment/uncomment variable definitions below

# WARNING: release,debug, and profile objects are in the same directory.
# This means you must 'make clean' before switching build modes.

all: tw3d_release copy_files
debug: tw3d_debug copy_files
profile: tw3d_profile copy_files

#################################################
# BEGIN INPUT VARIABLES BLOCK

# Uncomment exactly one platform

PLATFORM = OSX
#PLATFORM = WIN
#PLATFORM = LINUX
#PLATFORM = OPENMPI
#PLATFORM = INTELMPI
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

PACKAGE_PREF = HOMEBREW
#PACKAGE_PREF = MACPORTS

# END INPUT VARIABLES BLOCK
#################################################

ifndef HARDWARE_ACCEL
	HARDWARE_ACCEL = NONE
endif

ifndef PACKAGE_PREF
	PACKAGE_PREF = NONE
endif

ifeq ($(COMPILER_PREF),GNU)
	OMP_FLAG = -fopenmp
endif
ifeq ($(COMPILER_PREF),LLVM_CLANG)
	OMP_FLAG = -fopenmp
endif
ifeq ($(COMPILER_PREF),INTEL)
	OMP_FLAG = -qopenmp
endif

ifeq ($(PLATFORM),OSX)
	# Tested for Mac OS 10.15
	BASE_FLAGS = -c -std=c++17 -stdlib=libc++ -DUSE_DESKTOP
	LIBS = -lc++
	ifeq ($(HARDWARE_ACCEL),NONE)
		TW_Compiler = clang++
		TW_Linker = clang++
	endif
	ifeq ($(HARDWARE_ACCEL),APPLE_CL)
		TW_Compiler = clang++
		TW_Linker = clang++
		BASE_FLAGS += -DUSE_OPENCL
		LIBS += -framework OpenCL
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		BASE_FLAGS += -DUSE_OPENMP -fopenmp
		LIBS += -fopenmp
		ifeq ($(PACKAGE_PREF),MACPORTS)
			TW_Compiler = /opt/local/bin/clang++
			TW_Linker = /opt/local/bin/clang++
			BASE_FLAGS += -I/opt/local/include
			LIBS += -L/opt/local/lib -Wl,-rpath,/opt/local/lib
		endif
		ifeq ($(PACKAGE_PREF),HOMEBREW)
			TW_Compiler = /usr/local/opt/llvm/bin/clang++
			TW_Linker = /usr/local/opt/llvm/bin/clang++
			BASE_FLAGS += -isysroot /Applications/XCode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk
			BASE_FLAGS += -I/usr/local/opt/llvm/include/
			LIBS += -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib
		endif
	endif
	RELEASE_FLAGS = $(BASE_FLAGS) -O3 -march=native
	PROFILE_FLAGS = $(BASE_FLAGS) -O3 -march=native -Rpass-analysis=loop-vectorize
	DEBUG_FLAGS = $(BASE_FLAGS) -g -O0
	BINARY_PATH = ~/bin/
	WORK_PATH = ~/Run/
	TW_MPI = yes
endif

ifeq ($(PLATFORM),WIN)
	BASE_FLAGS = -c -std=c++17 -DUSE_DESKTOP -DSF_STD
	TW_Compiler = clang++
	TW_Linker = clang++
	LIBS =
	ifeq ($(HARDWARE_ACCEL),OMP)
		BASE_FLAGS += -DUSE_OPENMP -fopenmp
		LIBS += -fopenmp
	endif
	ifeq ($(HARDWARE_ACCEL),CUDA)
		CL_INCLUDE = /Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.0/include
		CL_LIB = /Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.0/lib/x64
		BASE_FLAGS += -DUSE_OPENCL -I'$(CL_INCLUDE)'
		LIBS +=  -L'$(CL_LIB)' -lOpenCL
	endif
	RELEASE_FLAGS = $(BASE_FLAGS) -O3 -march=native
	PROFILE_FLAGS = $(BASE_FLAGS) -O3 -march=native -Rpass-analysis=loop-vectorize
	DEBUG_FLAGS = $(BASE_FLAGS) -g -O0
	BINARY_PATH = $(HOMEPATH)\bin\$(strip)
	WORK_PATH = $(HOMEPATH)\Run\$(strip)
	TW_MPI = yes
endif

ifeq ($(PLATFORM),LINUX)
	BASE_FLAGS = -c -std=c++17 -DUSE_DESKTOP
	ifeq ($(COMPILER_PREF),LLVM_CLANG)
		TW_Compiler = clang++
		TW_Linker = clang++
		LIBS = -lc++
		BASE_FLAGS += -stdlib=libc++
	endif
	ifeq ($(COMPILER_PREF),GNU)
		TW_Compiler = g++
		TW_Linker = g++
		LIBS = -lstdc++
		BASE_FLAGS += -DSF_STD
	endif
	ifeq ($(COMPILER_PREF),INTEL)
		TW_Compiler = icpc
		TW_Linker = icpc
		LIBS = -lstdc++
		BASE_FLAGS += -DSF_STD
	endif
	ifeq ($(HARDWARE_ACCEL),NONE)
		LIBS += -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		BASE_FLAGS += -DUSE_OPENMP $(OMP_FLAG)
		LIBS += $(OMP_FLAG)
	endif
	ifeq ($(HARDWARE_ACCEL),CUDA)
		BASE_FLAGS += -DUSE_OPENCL
		LIBS += -lOpenCL -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),RADEON_PRO)
		BASE_FLAGS += -DUSE_OPENCL
		LIBS += -lOpenCL -pthread
	endif
	RELEASE_FLAGS = $(BASE_FLAGS) -O3 -march=native
	PROFILE_FLAGS = $(BASE_FLAGS) -O3 -march=native -Rpass-analysis=loop-vectorize
	DEBUG_FLAGS = $(BASE_FLAGS) -g -O0
	BINARY_PATH = ~/bin/
	WORK_PATH = ~/Run/
	TW_MPI = yes
endif

ifeq ($(PLATFORM),OPENMPI)
	BASE_FLAGS = -c -std=c++17 -DUSE_HPC
	TW_Compiler = mpicxx
	TW_Linker = mpicxx
	ifeq ($(COMPILER_PREF),LLVM_CLANG)
		export OMPI_CXX=clang++
		LIBS = -lc++
		BASE_FLAGS += -stdlib=libc++
	endif
	ifeq ($(COMPILER_PREF),INTEL)
		export OMPI_CXX=icpc
		LIBS = -lstdc++
	endif
	ifeq ($(COMPILER_PREF),GNU)
		export OMPI_CXX=g++
		LIBS = -lstdc++
	endif
	ifeq ($(HARDWARE_ACCEL),NONE)
		LIBS += -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		BASE_FLAGS += -DUSE_OPENMP $(OMP_FLAG)
		LIBS += $(OMP_FLAG)
	endif
	ifeq ($(HARDWARE_ACCEL),CUDA)
		BASE_FLAGS += -DUSE_OPENCL
		LIBS += -lOpenCL -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),RADEON_PRO)
		BASE_FLAGS += -DUSE_OPENCL
		LIBS += -lOpenCL -pthread
	endif
	RELEASE_FLAGS = $(BASE_FLAGS) -O3
	PROFILE_FLAGS = $(BASE_FLAGS) -O3 -Rpass-analysis=loop-vectorize
	DEBUG_FLAGS = $(BASE_FLAGS) -g -O0
endif

ifeq ($(PLATFORM),INTELMPI)
	BASE_FLAGS = -c -std=c++17 -DUSE_HPC
	ifeq ($(COMPILER_PREF),INTEL)
		TW_Compiler = mpiicpc
		TW_Linker = mpiicpc
		LIBS = -lstdc++
	endif
	ifeq ($(COMPILER_PREF),GNU)
		TW_Compiler = mpicxx
		TW_Linker = mpicxx
		LIBS = -lstdc++
	endif
	ifeq ($(HARDWARE_ACCEL),NONE)
		LIBS += -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		BASE_FLAGS += -DUSE_OPENMP $(OMP_FLAG)
		LIBS += $(OMP_FLAG)
	endif
	ifeq ($(HARDWARE_ACCEL),CUDA)
		BASE_FLAGS += -DUSE_OPENCL
		LIBS += -lOpenCL -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),RADEON_PRO)
		BASE_FLAGS += -DUSE_OPENCL
		LIBS += -lOpenCL -pthread
	endif
	RELEASE_FLAGS = $(BASE_FLAGS) -O3
	PROFILE_FLAGS = $(BASE_FLAGS) -O3 -Rpass-analysis=loop-vectorize
	DEBUG_FLAGS = $(BASE_FLAGS) -g -O0
endif

ifeq ($(PLATFORM),CRAY)
	BASE_FLAGS = -c -DUSE_HPC
	ifeq ($(COMPILER_PREF),CRAY)
		BASE_FLAGS += -h std=c++17
		TW_Compiler = CC
		TW_Linker = CC
		ifeq ($(HARDWARE_ACCEL),NONE)
			BASE_FLAGS += -h noomp
		endif
		ifeq ($(HARDWARE_ACCEL),OMP)
			BASE_FLAGS += -DUSE_OPENMP
		endif
		RELEASE_FLAGS = $(BASE_FLAGS) -O2
		PROFILE_FLAGS = $(BASE_FLAGS) -O2 -hlist=m
		DEBUG_FLAGS = $(BASE_FLAGS) -g -O0
		LIBS =
	else
		# DEFAULT TO INTEL
		BASE_FLAGS += -std=c++17 -DSF_STD
		TW_Compiler = CC
		TW_Linker = CC
		ifeq ($(HARDWARE_ACCEL),NONE)
			LIBS =
		endif
		ifeq ($(HARDWARE_ACCEL),OMP)
			BASE_FLAGS += -qopenmp -qopenmp-simd -DUSE_OPENMP
			LIBS = -qopenmp -qopenmp-simd
		endif
		RELEASE_FLAGS = $(BASE_FLAGS) -fast
		PROFILE_FLAGS = $(BASE_FLAGS) -qopt-report -qopt-report-phase=vec -fast
		DEBUG_FLAGS = $(BASE_FLAGS) -g
	endif
endif

ifeq ($(PLATFORM),SGI)
	BASE_FLAGS = -c -std=c++17 -DUSE_HPC
	TW_Compiler = icpc
	TW_Linker = icpc
	ifeq ($(HARDWARE_ACCEL),NONE)
		LIBS = -lmpi -lmpi++
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		BASE_FLAGS += -qopenmp -DUSE_OPENMP
		LIBS = -lmpi -lmpi++ -qopenmp
	endif
	ifeq ($(HARDWARE_ACCEL),CUDA)
		CUDA_ROOT = /p/home/apps/cuda/7.5
		BASE_FLAGS += -I$(CUDA_ROOT)/include -DUSE_OPENCL
		LIBS = -L$(CUDA_ROOT)/lib64
		LIBS += -lmpi -lmpi++ -lOpenCL
	endif
	RELEASE_FLAGS = $(BASE_FLAGS) -fast
	PROFILE_FLAGS = $(BASE_FLAGS) -qopt-report -qopt-report-phase=vec -fast
	DEBUG_FLAGS = $(BASE_FLAGS) -g
endif

ifndef PROFILE_FLAGS
	PROFILE_FLAGS = $(RELEASE_FLAGS)
endif

tw3d_release: CCFLAGS = $(RELEASE_FLAGS) -DNDEBUG
tw3d_release: tw3d

tw3d_debug: CCFLAGS = $(DEBUG_FLAGS)
tw3d_debug: tw3d

tw3d_profile: CCFLAGS = $(PROFILE_FLAGS)
tw3d_profile: tw3d

BASE_HEADERS = definitions.h tasks.h ctools.h 3dmath.h tw_iterator.h discreteSpace.h metricSpace.h functions.h fft.h numerics.h 3dfields.h region.h input.h units.h
TOOL_HEADERS = computeTool.h parabolic.h elliptic.h hyperbolic.h injection.h qstate.h physics.h chemistry.h diagnostics.h
SIM_HEADERS = $(BASE_HEADERS) $(TOOL_HEADERS) module.h simulation.h
MODULE_HEADERS = fieldSolve.h electrostatic.h laserSolve.h fluid.h quantum.h particles.h solidState.h
OBJ_LIST = 3DFields.o Chemistry.o ComputeTool.o Diagnostics.o DiscreteSpace.o Electrostatic.o Elliptic.o FFT.o FieldSolve.o Fluid.o Functions.o Hyperbolic.o Injection.o Input.o LaserSolve.o Main.o MetricSpace.o Module.o Numerics.o Parabolic.o Particles.o Physics.o Pusher.o Qstate.o Quantum.o Region.o Simulation.o SolidState.o Tasks.o Units.o
ifdef TW_MPI
	BASE_HEADERS += tw_mpi.h
	OBJ_LIST += TW_MPI.o
endif

tw3d: 3dfields.cl particles.cl elliptic.cl hyperbolic.cl quantum.cl $(OBJ_LIST)
	$(TW_Linker) -o tw3d $(OBJ_LIST) $(LIBS)

copy_files:
ifdef BINARY_PATH
ifeq ($(PLATFORM),WIN)
	copy tw3d $(BINARY_PATH)tw3d.exe
else
	cp tw3d $(BINARY_PATH)
endif
endif
ifdef WORK_PATH
ifeq ($(PLATFORM),WIN)
	copy *.cl $(WORK_PATH)
else
	cp *.cl $(WORK_PATH)
endif
endif

clean:
ifeq ($(PLATFORM),WIN)
	del *.o
	del tw3d
else
	rm -f *.o
	rm -f tw3d
endif

3DFields.o: 3DFields.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) 3DFields.cpp

Chemistry.o: Chemistry.cpp $(BASE_HEADERS) computeTool.h physics.h chemistry.h
	$(TW_Compiler) $(CCFLAGS) Chemistry.cpp

ComputeTool.o: ComputeTool.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) ComputeTool.cpp

Diagnostics.o: Diagnostics.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Diagnostics.cpp

DiscreteSpace.o: DiscreteSpace.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) DiscreteSpace.cpp

Electrostatic.o: Electrostatic.cpp $(SIM_HEADERS) fieldSolve.h electrostatic.h
	$(TW_Compiler) $(CCFLAGS) Electrostatic.cpp

Elliptic.o: Elliptic.cpp $(BASE_HEADERS) computeTool.h elliptic.h
	$(TW_Compiler) $(CCFLAGS) Elliptic.cpp

FFT.o: FFT.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) FFT.cpp

FieldSolve.o: FieldSolve.cpp $(SIM_HEADERS) fieldSolve.h
	$(TW_Compiler) $(CCFLAGS) FieldSolve.cpp

Fluid.o: Fluid.cpp $(SIM_HEADERS) fieldSolve.h fluid.h
	$(TW_Compiler) $(CCFLAGS) Fluid.cpp

Functions.o: Functions.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Functions.cpp

Hyperbolic.o: Hyperbolic.cpp $(BASE_HEADERS) computeTool.h hyperbolic.h
	$(TW_Compiler) $(CCFLAGS) Hyperbolic.cpp

Injection.o: Injection.cpp $(BASE_HEADERS) computeTool.h injection.h
	$(TW_Compiler) $(CCFLAGS) Injection.cpp

Input.o: Input.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Input.cpp

LaserSolve.o: LaserSolve.cpp $(SIM_HEADERS) fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) LaserSolve.cpp

Main.o: Main.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Main.cpp

MetricSpace.o: MetricSpace.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) MetricSpace.cpp

Module.o: Module.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Module.cpp

Numerics.o: Numerics.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Numerics.cpp

Parabolic.o: Parabolic.cpp $(BASE_HEADERS) computeTool.h parabolic.h
	$(TW_Compiler) $(CCFLAGS) Parabolic.cpp

Particles.o: Particles.cpp $(SIM_HEADERS) particles.h fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) Particles.cpp

Physics.o: Physics.cpp $(BASE_HEADERS) computeTool.h physics.h
	$(TW_Compiler) $(CCFLAGS) Physics.cpp

Pusher.o: Pusher.cpp $(SIM_HEADERS) particles.h fieldSolve.h laserSolve.h
	$(TW_Compiler) $(CCFLAGS) Pusher.cpp

Qstate.o: Qstate.cpp $(BASE_HEADERS) computeTool.h qstate.h
	$(TW_Compiler) $(CCFLAGS) Qstate.cpp

Quantum.o: Quantum.cpp $(SIM_HEADERS) quantum.h
	$(TW_Compiler) $(CCFLAGS) Quantum.cpp

Region.o: Region.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Region.cpp

Simulation.o: Simulation.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Simulation.cpp

SolidState.o: SolidState.cpp $(SIM_HEADERS) solidState.h
	$(TW_Compiler) $(CCFLAGS) SolidState.cpp

Tasks.o: Tasks.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Tasks.cpp

ifdef TW_MPI
TW_MPI.o: TW_MPI.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) TW_MPI.cpp
endif

Units.o: Units.cpp $(BASE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Units.cpp
