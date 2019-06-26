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

#PLATFORM = OSX
PLATFORM = LINUX
#PLATFORM = OPENMPI
#PLATFORM = CRAY
#PLATFORM = SGI

# Uncomment one, and at most one, of the following to use hardware acceleration
# N.b. specific accelerator devices can be specified in the turboWAVE input file

HARDWARE_ACCEL = OMP
#HARDWARE_ACCEL = APPLE_CL
#HARDWARE_ACCEL = CUDA
#HARDWARE_ACCEL = RADEON_PRO

# Specify a preference for compiler (ignored in some cases)

#COMPILER_PREF = GNU
COMPILER_PREF = LLVM_CLANG
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

ifndef PACKAGE_PREF
	PACKAGE_PREF = NONE
endif

ifeq ($(PLATFORM),OSX)
	ifeq ($(HARDWARE_ACCEL),NONE)
		TW_Compiler = clang++ -stdlib=libc++
		TW_Linker = clang++
		RELEASE_FLAGS = -c -std=c++11 -O3 -DUSE_DESKTOP
		DEBUG_FLAGS = -c -g -std=c++11 -O0 -DUSE_DESKTOP
		LIB_PATH =
		LIBS = -lc++
	endif
	ifeq ($(HARDWARE_ACCEL),APPLE_CL)
		TW_Compiler = clang++
		TW_Linker = clang++
		RELEASE_FLAGS = -c -std=c++11 -O3 -DUSE_DESKTOP -DUSE_OPENCL
		DEBUG_FLAGS = -c -g -std=c++11 -O0 -DUSE_DESKTOP -DUSE_OPENCL
		LIB_PATH =
		LIBS = -lc++ -framework OpenCL
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		ifeq ($(PACKAGE_PREF),NONE)
			TW_Compiler = "ERROR: bad settings."
		endif
		ifeq ($(PACKAGE_PREF),MACPORTS)
			# Tested for Mac OS 10.14
			# Remember to run port select clang mp-clang-x.x
			TW_Compiler = /opt/local/bin/clang++ -stdlib=libc++
			TW_Linker = /opt/local/bin/clang++
			RELEASE_FLAGS = -c -O3 -std=c++11 -march=native -fopenmp -I/opt/local/include -DUSE_DESKTOP -DUSE_OPENMP
			PROFILE_FLAGS = -c -O3 -std=c++11 -fopenmp -Rpass-analysis=loop-vectorize -I/opt/local/include -DUSE_DESKTOP -DUSE_OPENMP
			DEBUG_FLAGS = -c -g -O0 -std=c++11 -fopenmp -I/opt/local/include -DUSE_DESKTOP -DUSE_OPENMP
			LIB_PATH = -L/opt/local/lib -Wl,-rpath,/opt/local/lib
			LIBS = -lc++ -fopenmp
		endif
		ifeq ($(PACKAGE_PREF),HOMEBREW)
			# Tested for Mac OS 10.13
			TW_Compiler = /usr/local/opt/llvm/bin/clang++ -stdlib=libc++
			TW_Linker = /usr/local/opt/llvm/bin/clang++
			RELEASE_FLAGS = -c -O3 -std=c++11 -march=native -fopenmp -I/usr/local/opt/llvm/include/ -DUSE_DESKTOP -DUSE_OPENMP
			PROFILE_FLAGS = -c -O3 -std=c++11 -fopenmp -Rpass-analysis=loop-vectorize -I/usr/local/opt/llvm/include -DUSE_DESKTOP -DUSE_OPENMP
			DEBUG_FLAGS = -c -g -O0 -std=c++11 -fopenmp -I/usr/local/opt/llvm/include/ -DUSE_DESKTOP -DUSE_OPENMP
			LIB_PATH = -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib
			LIBS = -lc++ -fopenmp
		endif
	endif
	BINARY_PATH = ~/bin
	WORK_PATH = ~/Run
	TW_MPI = yes
endif

ifeq ($(PLATFORM),LINUX)
	ifeq ($(COMPILER_PREF),LLVM_CLANG)
		TW_Compiler = clang++ -stdlib=libc++
		TW_Linker = clang++
		LIBS0 = -lc++
	endif
	ifeq ($(COMPILER_PREF),GNU)
		TW_Compiler = g++
		TW_Linker = g++
		LIBS0 = -lstdc++
	endif
	ifeq ($(COMPILER_PREF),INTEL)
		TW_Compiler = icpc
		TW_Linker = icpc
		LIBS0 = -lstdc++
	endif
	ifeq ($(HARDWARE_ACCEL),NONE)
		RELEASE_FLAGS = -c -std=c++11 -O3 -DUSE_DESKTOP
		DEBUG_FLAGS = -c -g -std=c++11 -O0 -DUSE_DESKTOP
		LIB_PATH =
		LIBS = $(LIBS0) -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		RELEASE_FLAGS = -c -std=c++11 -O3 -march=native -fopenmp -DUSE_DESKTOP -DUSE_OPENMP
		PROFILE_FLAGS = -c -std=c++11 -O3 -march=native -fopenmp -Rpass-analysis=loop-vectorize -DUSE_DESKTOP -DUSE_OPENMP
		DEBUG_FLAGS = -c -g -std=c++11 -O0 -fopenmp -DUSE_DESKTOP -DUSE_OPENMP
		LIB_PATH =
		LIBS = $(LIBS0) -fopenmp
	endif
	ifeq ($(HARDWARE_ACCEL),CUDA)
		RELEASE_FLAGS = -c -std=c++11 -O3 -DUSE_DESKTOP -DUSE_OPENCL
		DEBUG_FLAGS = -c -g -std=c++11 -O0 -DUSE_DESKTOP -DUSE_OPENCL
		LIB_PATH =
		LIBS = $(LIBS0) -lOpenCL -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),RADEON_PRO)
		RELEASE_FLAGS = -c -std=c++11 -O3 -DUSE_DESKTOP -DUSE_OPENCL
		DEBUG_FLAGS = -c -g -std=c++11 -O0 -DUSE_DESKTOP -DUSE_OPENCL
		LIB_PATH =
		LIBS = $(LIBS0) -lOpenCL -pthread
	endif
	BINARY_PATH = ~/bin
	WORK_PATH = ~/Run
	TW_MPI = yes
endif

ifeq ($(PLATFORM),OPENMPI)
	TW_Compiler = mpicxx
	TW_Linker = mpicxx
	ifeq ($(COMPILER_PREF),LLVM_CLANG)
		export OMPI_CXX=clang++ -stdlib=libc++
		LIBS0 = -lc++
	endif
	ifeq ($(COMPILER_PREF),INTEL)
		export OMPI_CXX=icpc
		LIBS0 = -lstdc++
	endif
	ifeq ($(COMPILER_PREF),GNU)
		export OMPI_CXX=g++
		LIBS0 = -lstdc++
	endif
	ifeq ($(HARDWARE_ACCEL),NONE)
		RELEASE_FLAGS = -c -std=c++11 -O3 -DUSE_HPC
		DEBUG_FLAGS = -c -g -std=c++11 -O0 -DUSE_HPC
		LIB_PATH =
		LIBS = $(LIBS0) -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		RELEASE_FLAGS = -c -std=c++11 -O3 -march=native -fopenmp -DUSE_HPC -DUSE_OPENMP
		PROFILE_FLAGS = -c -std=c++11 -O3 -march=native -fopenmp -Rpass-analysis=loop-vectorize -DUSE_HPC -DUSE_OPENMP
		DEBUG_FLAGS = -c -g -std=c++11 -O0 -fopenmp -DUSE_HPC -DUSE_OPENMP
		LIB_PATH =
		LIBS = $(LIBS0) -fopenmp
	endif
	ifeq ($(HARDWARE_ACCEL),CUDA)
		RELEASE_FLAGS = -c -std=c++11 -O3 -DUSE_HPC -DUSE_OPENCL
		DEBUG_FLAGS = -c -g -std=c++11 -O0 -DUSE_HPC -DUSE_OPENCL
		LIB_PATH =
		LIBS = $(LIBS0) -lOpenCL -pthread
	endif
	ifeq ($(HARDWARE_ACCEL),RADEON_PRO)
		RELEASE_FLAGS = -c -std=c++11 -O3 -DUSE_HPC -DUSE_OPENCL
		DEBUG_FLAGS = -c -g -std=c++11 -O0 -DUSE_HPC -DUSE_OPENCL
		LIB_PATH =
		LIBS = $(LIBS0) -lOpenCL -pthread
	endif
	BINARY_PATH = ~/bin
	WORK_PATH = ~/Run
endif

ifeq ($(PLATFORM),CRAY)
	ifeq ($(COMPILER_PREF),CRAY)
		ifeq ($(HARDWARE_ACCEL),NONE)
			TW_Compiler = CC
			TW_Linker = CC
			RELEASE_FLAGS = -c -O2 -h std=c++11 -h noomp -DUSE_HPC
			DEBUG_FLAGS = -c -g -O0 -h std=c++11 -h noomp -DUSE_HPC
			LIB_PATH =
			LIBS =
		endif
		ifeq ($(HARDWARE_ACCEL),OMP)
			TW_Compiler = CC
			TW_Linker = CC
			RELEASE_FLAGS = -c -O2 -h std=c++11 -DUSE_HPC -DUSE_OPENMP
			PROFILE_FLAGS = -c -O2 -hlist=m -h std=c++11 -DUSE_HPC -DUSE_OPENMP
			DEBUG_FLAGS = -c -g -O0 -h std=c++11 -DUSE_HPC -DUSE_OPENMP
			LIB_PATH =
			LIBS =
		endif
	else
		# DEFAULT TO INTEL
		ifeq ($(HARDWARE_ACCEL),NONE)
			TW_Compiler = CC
			TW_Linker = CC
			RELEASE_FLAGS = -c -std=c++11 -fast -DUSE_HPC
			DEBUG_FLAGS = -c -g -std=c++11 -DUSE_HPC
			LIB_PATH =
			LIBS =
		endif
		ifeq ($(HARDWARE_ACCEL),OMP)
			TW_Compiler = CC
			TW_Linker = CC
			RELEASE_FLAGS = -c -std=c++11 -qopenmp -qopenmp-simd -fast -DUSE_HPC -DUSE_OPENMP
			PROFILE_FLAGS = -c -std=c++11 -qopenmp -qopenmp-simd -qopt-report -qopt-report-phase=vec -fast -DUSE_HPC -DUSE_OPENMP
			DEBUG_FLAGS = -c -g -std=c++11 -qopenmp -qopenmp-simd -DUSE_HPC -DUSE_OPENMP
			LIB_PATH =
			LIBS = -qopenmp -qopenmp-simd
			LD_PROFILE_FLAGS = -qopt-report -qopt-report-phase=vec
		endif
	endif
endif

ifeq ($(PLATFORM),SGI)
	ifeq ($(HARDWARE_ACCEL),NONE)
		TW_Compiler = icpc
		TW_Linker = icpc
		RELEASE_FLAGS = -c -std=c++11 -fast -DUSE_HPC
		DEBUG_FLAGS = -c -g -std=c++11 -DUSE_HPC
		LIB_PATH =
		LIBS = -lmpi -lmpi++
	endif
	ifeq ($(HARDWARE_ACCEL),OMP)
		TW_Compiler = icpc
		TW_Linker = icpc
		RELEASE_FLAGS = -c -std=c++11 -qopenmp -fast -DUSE_HPC -DUSE_OPENMP
		DEBUG_FLAGS = -c -g -std=c++11 -qopenmp -DUSE_HPC -DUSE_OPENMP
		LIB_PATH =
		LIBS = -lmpi -lmpi++ -qopenmp
	endif
	ifeq ($(HARDWARE_ACCEL),CUDA)
		CUDA_ROOT = /p/home/apps/cuda/7.5
		TW_Compiler = icpc
		TW_Linker = icpc
		RELEASE_FLAGS = -c -std=c++11 -fast -I$(CUDA_ROOT)/include -DUSE_HPC -DUSE_OPENCL
		DEBUG_FLAGS = -c -g -std=c++11 -I$(CUDA_ROOT)/include -DUSE_HPC -DUSE_OPENCL
		LIB_PATH = -L$(CUDA_ROOT)/lib64
		LIBS = -lmpi -lmpi++ -lOpenCL
	endif
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

BASE_HEADERS = definitions.h tasks.h ctools.h 3dmath.h metricSpace.h 3dfields.h region.h tw_iterator.h
SIM_HEADERS = $(BASE_HEADERS) functions.h physics.h chemistry.h numerics.h computeTool.h elliptic.h parabolic.h hyperbolic.h fft.h injection.h input.h diagnostics.h module.h simulation.h
MODULE_HEADERS = fieldSolve.h electrostatic.h laserSolve.h fluid.h quantum.h particles.h solidState.h
OBJ_LIST = Main.o Tasks.o MetricSpace.o 3DFields.o Region.o FFT.o Physics.o Chemistry.o Numerics.o ComputeTool.o Parabolic.o Elliptic.o Hyperbolic.o Functions.o Simulation.o Module.o FieldSolve.o Electrostatic.o LaserSolve.o Fluid.o Quantum.o SolidState.o Particles.o Pusher.o Diagnostics.o Injection.o Input.o
ifdef TW_MPI
	BASE_HEADERS += tw_mpi.h
	OBJ_LIST += TW_MPI.o
endif

tw3d: 3dfields.cl particles.cl elliptic.cl hyperbolic.cl quantum.cl $(OBJ_LIST)
	$(TW_Linker) -o tw3d $(OBJ_LIST) $(LIB_PATH) $(LIBS)

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

Physics.o: Physics.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Physics.cpp

Chemistry.o: Chemistry.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Chemistry.cpp

Numerics.o: Numerics.cpp $(BASE_HEADERS) numerics.h
	$(TW_Compiler) $(CCFLAGS) Numerics.cpp

ComputeTool.o: ComputeTool.cpp $(SIM_HEADERS)
	$(TW_Compiler) $(CCFLAGS) ComputeTool.cpp

Elliptic.o: Elliptic.cpp $(BASE_HEADERS) numerics.h computeTool.h elliptic.h
	$(TW_Compiler) $(CCFLAGS) Elliptic.cpp

Parabolic.o: Parabolic.cpp $(BASE_HEADERS) numerics.h computeTool.h parabolic.h
	$(TW_Compiler) $(CCFLAGS) Parabolic.cpp

Hyperbolic.o: Hyperbolic.cpp $(BASE_HEADERS) numerics.h computeTool.h hyperbolic.h
	$(TW_Compiler) $(CCFLAGS) Hyperbolic.cpp

Functions.o: Functions.cpp definitions.h ctools.h functions.h
	$(TW_Compiler) $(CCFLAGS) Functions.cpp

Simulation.o: Simulation.cpp $(SIM_HEADERS) $(MODULE_HEADERS)
	$(TW_Compiler) $(CCFLAGS) Simulation.cpp

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
