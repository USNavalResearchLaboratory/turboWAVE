cmake_minimum_required(VERSION 4.0.3)

# Build Reminders
#
# For MSVC, we have to specify release at build time `cmake --build . --config Release`
# For non-MSVC, we specify release at config time `cmake [options] -DCMAKE_BUILD_TYPE=Release ..`
# For clangd support, add -DCMAKE_EXPORT_COMPILE_COMMANDS=1 during config
#
# Runtime Reminders
#
# For Mac/Linux the OpenMP shared library may not be found.
# This is fixed using an environment variable, e.g.:
# export LD_LIBRARY_PATH=$HOMEBREW_PREFIX/opt/llvm/lib/x86_64-unknown-linux-gnu
#
# About `import std`
#
# To tell TW to use the `std` module instead of includes, define USE_STD_MODULE in base/config.h.in.
# As of this writing CMake is not very much help in finding/generating the PCM.
# To get it working manually:
# Compile the CPPM, for Apple/Homebrew it is at $HOMEBREW_PREFIX/opt/llvm/share/libc++/v1/.
# The compile option must be consistent with the rest of the project.
#    example: clang++ -std=gnu++23 -stdlib=libc++ -arch arm64 -fopenmp=libomp -Wno-reserved-module-identifier --precompile -o std.pcm std.cppm
# Once we have the PCM we point to it by adding to CMakeLists the line:
#    target_compile_options(tw3d PRIVATE -fmodule-file=std=/full/path/to/std.pcm)
# It is hoped someday CMake will take care of all that.  As part of this hope, compilers
# provide libc++.modules.json to help the build system automate things.
# But it can happen that this file itself cannot be located.  So Kitware provided
# CMAKE_CXX_STDLIB_MODULES_JSON as a way for users to point to this explicitly.

if(WIN32)
    # Windows requires installing some things with winget:
    # winget install LLVM.LLVM
    # winget install Kitware.CMake
    # winget install Microsoft.msmpisdk
    # winget install Microsoft.msmpi
    # winget install Ninja-build.Ninja
    # set(CMAKE_CXX_STDLIB_MODULES_JSON "$ENV{ProgramFiles}/Microsoft Visual Studio/18/Community/VC/Tools/MSVC/14.50/35717/modules/modules.json")
    set(CMAKE_CXX_COMPILER "$ENV{ProgramFiles}/llvm/bin/clang++.exe")
    set(CMAKE_C_COMPILER "$ENV{ProgramFiles}/llvm/bin/clang.exe")
    set(CMAKE_SYSTEM_LIBRARY_PATH "$ENV{ProgramFiles}/LLVM/lib")
    set(CMAKE_SYSTEM_INCLUDE_PATH "$ENV{ProgramFiles}/LLVM/include")
elseif(APPLE)
    # Apple requires installing some things with homebrew:
    # brew install llvm
    # brew install cmake
    # brew install libomp (maybe packed with llvm)
    # brew install open-mpi (or substitute favorite MPI)
    # brew install ninja
    # set(CMAKE_CXX_STDLIB_MODULES_JSON "$ENV{HOMEBREW_PREFIX}/Cellar/llvm/21.1.6/lib/c++/libc++.modules.json")
	set(CMAKE_CXX_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/llvm/bin/clang++")
	set(CMAKE_C_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/llvm/bin/clang")
	# set(CMAKE_SYSTEM_LIBRARY_PATH "$ENV{HOMEBREW_PREFIX}/opt/llvm/lib")
	# set(CMAKE_SYSTEM_INCLUDE_PATH "$ENV{HOMEBREW_PREFIX}/opt/llvm/include")
	# set(CMAKE_CXX_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/gcc/bin/g++-15")
	# set(CMAKE_C_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/gcc/bin/gcc-15")
	# set(CMAKE_SYSTEM_LIBRARY_PATH "$ENV{HOMEBREW_PREFIX}/opt/gcc/lib")
	# set(CMAKE_SYSTEM_INCLUDE_PATH "$ENV{HOMEBREW_PREFIX}/opt/gcc/include")
    # For Apple we may be forced to use libc++
    # set(CMAKE_CXX_FLAGS "-stdlib=libstdc++")
    set(CMAKE_CXX_FLAGS "-stdlib=libc++")
    # It seems that libstdc++ handles c++17 special functions correctly and libc++ does not.
    # So if we are using libc++ we have to use our own wrappers.
    add_definitions(-DUSE_TWSF)
elseif(UNIX)
    # Ubuntu requires the following installations:
    #   maybe `ln -s /var/lib/snapd/snap snap`
    #   install snapd and homebrew
    #   snap install cmake --classic
    #   apt install ninja-build
    #   apt install gcc-14
    #   apt install libstdc++-14-dev
    #   apt install libopenmpi-dev
    #   brew install llvm
    # Fedora requires the following installations:
    #   maybe `ln -s /var/lib/snapd/snap snap`
    #   install snapd and homebrew
    #   snap install cmake --classic
    #   dnf install ninja-build
    #   dnf install gcc
    #   dnf install g++
    #   dnf install openmpi-devel
    #   brew install llvm
    #   module load mpi
	set(CMAKE_CXX_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/llvm/bin/clang++")
	set(CMAKE_C_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/llvm/bin/clang")
	# set(CMAKE_SYSTEM_LIBRARY_PATH "$ENV{HOMEBREW_PREFIX}/opt/llvm/lib")
	# set(CMAKE_SYSTEM_INCLUDE_PATH "$ENV{HOMEBREW_PREFIX}/opt/llvm/include")
endif()

set(CMAKE_CXX_STANDARD 23)
# set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")

project(turbowave
        VERSION "5.0.0"
        LANGUAGES C CXX
)
add_executable(tw3d)

# set(CMAKE_CXX_MODULE_STD ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS ON)

configure_file(base/config.h.in config.h)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# Regular sources or module implementations go here
target_sources(tw3d PUBLIC
    base/Main.cpp
    io/parser.c
    io/tree-sitter/lib/src/lib.c
    io/ReadSim.cpp
    io/ReadDriver.cpp
    particles/Mover_Bohmian.cpp
    particles/Mover_Boris.cpp
    particles/Mover_HC.cpp
    particles/Mover_PGC.cpp
    particles/Mover_Photon.cpp
    particles/Mover_Unitary.cpp
    test/mover_test.cpp
    test/particles_test.cpp
    test/pusher_test.cpp
    test/elliptic_test.cpp
    test/fluid_test.cpp
    test/physics_test.cpp
    test/injection_test.cpp
)

target_sources(tw3d PUBLIC
    FILE_SET HEADERS
    BASE_DIRS base io particles fields solver sparc test
)

# Module interfaces go here
target_sources(tw3d PUBLIC
    FILE_SET CXX_MODULES
    BASE_DIRS base driver io particles fields solver sparc test
    FILES
        base/base.cpp
        base/logger.cpp
        base/pic_primitives.cpp
        base/tensor.cpp
        base/StaticSpace.cpp
        base/DynSpace.cpp
        base/MetricSpace.cpp
        base/tw_iterator.cpp
        base/Factory.cpp
        base/FFT.cpp
        base/Functions.cpp
        base/Tasks.cpp
        driver/Driver.cpp
        driver/Diagnostic.cpp
        driver/Engine.cpp
        driver/Profile.cpp
        driver/Region.cpp
        driver/Warp.cpp
        driver/Tool.cpp
        driver/Bounded.cpp
        driver/Simulation.cpp
        fields/Fields.cpp
        fields/primitives.cpp
        fields/base.cpp
        fields/vector_ops.cpp
        fields/transform.cpp
        fields/aggregates.cpp
        io/Navigate.cpp
        io/Preprocess.cpp
        io/Assignment.cpp
        io/Input.cpp
        io/Diagnostics.cpp
        io/Units.cpp
        io/ReadGrid.cpp
        particles/Bundle.cpp
        particles/Tiler.cpp
        particles/Pusher.cpp
        particles/Mover.cpp
        particles/Particles.cpp
        particles/Qed.cpp
        sparc/Physics.cpp
        sparc/fct.cpp
        sparc/Chemistry.cpp
        sparc/Fluid.cpp
        solver/numerics_serial.cpp
        solver/numerics_parallel.cpp
        solver/Numerics.cpp
        solver/Injection.cpp
        solver/SolidState.cpp
        solver/Parabolic.cpp
        solver/Qstate.cpp
        solver/Quantum.cpp
        solver/FieldSolve.cpp
        solver/Hyperbolic.cpp
        solver/LaserSolve.cpp
        solver/Electrostatic.cpp
        solver/Elliptic.cpp
        test/iterator_test.cpp
        test/metric_space_test.cpp
        test/fft_test.cpp
)

if (MSVC)
    target_compile_options(tw3d PUBLIC -openmp:experimental -openmp:llvm)
endif()

if (APPLE)
    target_compile_options(tw3d PRIVATE -fmodule-file=std=/opt/homebrew/opt/llvm/share/libc++/v1/std.pcm)
endif()

target_include_directories(tw3d PUBLIC
    io/tree-sitter/lib/include
    io/tree-sitter/lib/src
    io/tree-sitter/lib/src/unicode
    ${MPI_CXX_INCLUDE_PATH}
    ${PROJECT_BINARY_DIR} # for config.h
)

target_link_libraries(tw3d PUBLIC
    OpenMP::OpenMP_CXX
    ${MPI_CXX_LIBRARIES}
)
