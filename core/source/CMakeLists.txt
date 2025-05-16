cmake_minimum_required(VERSION 3.30)
if(WIN32)
    # This requires installing some things with winget:
    # winget install LLVM.LLVM
    # winget install Kitware.CMake
    # winget install Microsoft.msmpisdk
    # winget install Microsoft.msmpi
    # winget install Ninja-build.Ninja
    set(CMAKE_CXX_COMPILER "$ENV{ProgramFiles}/llvm/bin/clang++.exe")
    set(CMAKE_C_COMPILER "$ENV{ProgramFiles}/llvm/bin/clang.exe")
    set(CMAKE_SYSTEM_LIBRARY_PATH "$ENV{ProgramFiles}/LLVM/lib")
endif()
if(APPLE OR UNIX)
    # Supposing homebrew/linuxbrew is used for LLVM:
    # brew install llvm
    # brew install cmake
    # brew install libomp (maybe packed with llvm)
    # brew install open-mpi (or substitute favorite MPI)
    # brew install ninja
    # maybe necessary to set LD_LIBRARY_PATH to find libomp.so
	set(CMAKE_CXX_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/llvm/bin/clang++")
	set(CMAKE_C_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/llvm/bin/clang")
	set(CMAKE_SYSTEM_LIBRARY_PATH "$ENV{HOMEBREW_PREFIX}/opt/llvm/lib")
	set(CMAKE_SYSTEM_INCLUDE_PATH "$ENV{HOMEBREW_PREFIX}/opt/llvm/include")
    # set(CMAKE_CXX_FLAGS "-stdlib=libc++")
endif()
if(APPLE)
    # Shocking we still need this in 2025:
    # Use if errors are encounted pertaining to `j0`, `j1`, or `std::some_special_function`.
    # It appears the LLVM libc++ creates this situation.
    add_definitions(-DUSE_TWSF)
endif()

# set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD ON)
# set(CMAKE_CXX_MODULE_STD ON)

project(turbowave
        VERSION "5.0.0"
        LANGUAGES C CXX
)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS ON)

add_executable(tw3d)

configure_file(base/config.h.in config.h)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

set(src_base
    base/Main.cpp
    base/Simulation.cpp
    base/Injection_Test.cpp
    base/MetricSpace_Test.cpp
)
set(src_io
    io/parser.c
    io/tree-sitter/lib/src/lib.c
    io/MainReader.cpp
)
set(src_par
    particles/Mover_Test.cpp
    particles/Particles_Test.cpp
    particles/Pusher_Test.cpp
    particles/Mover_Bohmian.cpp
    particles/Mover_Boris.cpp
    particles/Mover_HC.cpp
    particles/Mover_PGC.cpp
    particles/Mover_Photon.cpp
    particles/Mover_Unitary.cpp
)
set(src_solv
    solver/Elliptic_Test.cpp
)
set(src_sparc
    sparc/Fluid_Test.cpp
    sparc/Physics_Test.cpp
)

target_sources(tw3d PUBLIC
    ${src_base}
    ${src_io}
    ${src_par}
    ${src_solv}
    ${src_sparc}
)

target_sources(tw3d PUBLIC
    FILE_SET HEADERS
    BASE_DIRS base io particles solver sparc
)

target_sources(tw3d PUBLIC
    FILE_SET CXX_MODULES
    BASE_DIRS base io particles solver sparc
    FILES
        base/base.cpp
        base/pic_primitives.cpp
        base/tensor.cpp
        base/DiscreteSpace.cpp
        base/tw_iterator.cpp
        base/MetricSpace.cpp
        base/Region.cpp
        base/ComputeTool.cpp
        base/Module.cpp
        base/Injection.cpp
        base/Factory.cpp
        base/FFT.cpp
        base/Fields.cpp
        base/Functions.cpp
        base/Tasks.cpp
        io/Navigate.cpp
        io/Preprocess.cpp
        io/Assignment.cpp
        io/Input.cpp
        io/Diagnostics.cpp
        io/Units.cpp
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
        solver/Numerics.cpp
        solver/SolidState.cpp
        solver/Parabolic.cpp
        solver/Qstate.cpp
        solver/Quantum.cpp
        solver/FieldSolve.cpp
        solver/Hyperbolic.cpp
        solver/LaserSolve.cpp
        solver/Electrostatic.cpp
        solver/Elliptic.cpp)

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
