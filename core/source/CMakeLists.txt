cmake_minimum_required(VERSION 3.30)
if(WIN32)
    # LLVM and Ninja can be installed using winget
    set(CMAKE_CXX_COMPILER "$ENV{ProgramFiles}/llvm/bin/clang++.exe")
    set(CMAKE_C_COMPILER "$ENV{ProgramFiles}/llvm/bin/clang.exe")
    set(CMAKE_SYSTEM_LIBRARY_PATH "$ENV{ProgramFiles}/LLVM/lib")
endif()
if(APPLE)
    # for this we will need to use
    # brew install llvm
    # brew install libomp (maybe packed with llvm)
    # brew install cmake
    # brew install ninja
    set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
    set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
    set(CMAKE_SYSTEM_LIBRARY_PATH "/opt/homebrew/opt/llvm/lib")
    set(CMAKE_SYSTEM_INCLUDE_PATH "/opt/homebrew/opt/llvm/include")
endif()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(turbowave
        VERSION "5.0.0"
        LANGUAGES C CXX)

find_package(OpenMP REQUIRED)

set(src_base
    base/Simulation.cpp
    base/FFT.cpp
    base/Functions.cpp
    base/Injection_Test.cpp
    base/Main.cpp
    base/MetricSpace_Test.cpp
)
set(src_con
    concurrency/Tasks.cpp
)
set(src_io
    io/Units.cpp io/parser.c io/tree-sitter/lib/src/lib.c
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

add_executable(tw3d)

target_sources(tw3d PUBLIC
    ${src_base}
    ${src_con}
    ${src_io}
    ${src_par}
    ${src_solv}
    ${src_sparc}
    concurrency/TW_MPI.cpp
)

target_sources(tw3d PUBLIC
    FILE_SET HEADERS
    BASE_DIRS base concurrency io particles solver sparc
)

target_sources(tw3d PUBLIC
    FILE_SET CXX_MODULES
    BASE_DIRS base concurrency io particles solver sparc
    FILES
        base/DiscreteSpace.cpp
        base/tw_iterator.cpp
        base/MetricSpace.cpp
        base/Region.cpp
        base/ComputeTool.cpp
        base/Module.cpp
        base/Injection.cpp
        base/Factory.cpp
        base/3DFields.cpp
        io/Assignment.cpp
        io/Input.cpp
        io/Diagnostics.cpp
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
)

target_link_libraries(tw3d PUBLIC OpenMP::OpenMP_CXX)
