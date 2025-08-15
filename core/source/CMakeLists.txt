cmake_minimum_required(VERSION 4.0.3)

# REMINDERS
#
# For MSVC, we have to specify release at build time `cmake --build . --config Release`
# For non-MSVC, we specify release at config time `cmake [options] -DCMAKE_BUILD_TYPE=Release ..`
# For clangd support, add -DCMAKE_EXPORT_COMPILE_COMMANDS=1 during config
# it may be necessary to build std by hand
# 
# To build module std for LLVM libc++:
# git clone https://github.com/llvm/llvm-project.git
# cd llvm-project
# mkdir build
# cmake -G Ninja -S runtimes -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind" -B build -DCMAKE_CXX_COMPILER=/home/linuxbrew/.linuxbrew/opt/llvm/bin/clang++ -DCMAKE_C_COMPILER=/home/linuxbrew/.linuxbrew/opt/llvm/bin/clang
# ninja -C build
# this puts std.cppm into llvm-project/build/modules/c++/v1

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
    set(CMAKE_SYSTEM_INCLUDE_PATH "$ENV{ProgramFiles}/LLVM/include")
elseif(APPLE)
    # brew install llvm
    # brew install cmake
    # brew install libomp (maybe packed with llvm)
    # brew install open-mpi (or substitute favorite MPI)
    # brew install ninja
    # maybe necessary to set LD_LIBRARY_PATH to find libomp.so
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
    # set(CMAKE_CXX_FLAGS "-stdlib=libc++")
    # It seems that libstdc++ handles c++17 special functions correctly and libc++ does not.
    # So if we are using libc++ we have to use our own wrappers.
    # add_definitions(-DUSE_TWSF)
elseif(UNIX)
    # n.b. we need snap, apt, and linuxbrew.
    # snap install cmake
    # apt install ninja-build
    # apt install gcc-14
    # apt install libstdc++-14-dev
    # apt install libopenmpi-dev
    # brew install llvm
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
    base/Simulation.cpp
    io/parser.c
    io/tree-sitter/lib/src/lib.c
    io/MainReader.cpp
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
    BASE_DIRS base io particles fields solver sparc test
    FILES
        base/base.cpp
        base/logger.cpp
        base/pic_primitives.cpp
        base/tensor.cpp
        base/StaticSpace.cpp
        base/DynSpace.cpp
        base/MetricSpace.cpp
        base/tw_iterator.cpp
        base/Region.cpp
        base/ComputeTool.cpp
        base/Module.cpp
        base/Injection.cpp
        base/Factory.cpp
        base/FFT.cpp
        base/Functions.cpp
        base/Tasks.cpp
        fields/Fields.cpp
        fields/primitives.cpp
        fields/base.cpp
        fields/vector_ops.cpp
        fields/transform.cpp
        fields/aggregates.cpp
        fields/tools.cpp
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
        solver/numerics_serial.cpp
        solver/numerics_parallel.cpp
        solver/Numerics.cpp
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
