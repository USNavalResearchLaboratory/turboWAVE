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
# Standard Library
#
# We are compiling the standard library as a user module as a workaround for the fact
# that the build system is not yet able to handle `import std`.  This leads to warnings
# since the compiler thinks we are doing unusual things with a user module.
# It is possible to revert to includes by editing `config.h.in` and a few lines here.

# Uncomment exactly one of the following
set(TW_COMPILER clang)
# set(TW_COMPILER gcc)
# set(TW_COMPILER msvc)

if(WIN32)
    # Windows requires installing some things with winget:
    # winget install LLVM.LLVM
    # winget install Kitware.CMake
    # winget install Microsoft.msmpisdk
    # winget install Microsoft.msmpi
    # winget install Ninja-build.Ninja
    if(TW_COMPILER STREQUAL clang)
        set(CMAKE_CXX_COMPILER "$ENV{ProgramFiles}/llvm/bin/clang++.exe")
        set(CMAKE_C_COMPILER "$ENV{ProgramFiles}/llvm/bin/clang.exe")
        set(CMAKE_SYSTEM_LIBRARY_PATH "$ENV{ProgramFiles}/LLVM/lib")
        set(CMAKE_SYSTEM_INCLUDE_PATH "$ENV{ProgramFiles}/LLVM/include")
    endif()
    set(TW_STD_MODULE_PATH "C:/Program Files/Microsoft Visual Studio/18/Community/VC/Tools/MSVC/14.50.35717/modules")
    set(TW_STD_MODULE_FILE "std.ixx")
elseif(APPLE)
    # Apple requires installing some things with homebrew:
    # brew install llvm
    # brew install cmake
    # brew install libomp (maybe packed with llvm)
    # brew install open-mpi (or substitute favorite MPI)
    # brew install ninja
    if(TW_COMPILER STREQUAL clang)
    	set(CMAKE_CXX_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/llvm/bin/clang++")
	    set(CMAKE_C_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/llvm/bin/clang")
        set(CMAKE_CXX_FLAGS "-stdlib=libc++")
        # It seems that libstdc++ handles c++17 special functions correctly and libc++ does not.
        # So if we are using libc++ we have to use our own wrappers, and incur some limitations.
        add_definitions(-DUSE_TWSF)
        set(TW_STD_MODULE_PATH "$ENV{HOMEBREW_PREFIX}/opt/llvm/share/libc++/v1")
        set(TW_STD_MODULE_FILE "std.cppm")
    endif()
    if(TW_COMPILER STREQUAL gcc)
	    set(CMAKE_CXX_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/gcc/bin/g++-15")
	    set(CMAKE_C_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/gcc/bin/gcc-15")
        set(TW_STD_MODULE_PATH "$ENV{HOMEBREW_PREFIX}/include/c++/15/bits")
        set(TW_STD_MODULE_FILE "std.cc")
    endif()
elseif(UNIX)
    # Ubuntu requires the following installations:
    #   maybe `ln -s /var/lib/snapd/snap snap`
    #   install snapd and homebrew
    #   snap install cmake --classic
    #   apt install ninja-build
    #   apt install libopenmpi-dev
    #   apt install libc++-20-dev
    #   brew install llvm (>=21)
    #   brew install gcc (>=15)
    # Fedora 43 requires the following installations:
    #   maybe `ln -s /var/lib/snapd/snap snap`
    #   install snapd and homebrew
    #   snap install cmake --classic
    #   dnf install ninja-build
    #   dnf install openmpi-devel
    #   dnf install libcxx-devel
    #   brew install llvm (>=21)
    #   brew install gcc (>=15)
    #   module load mpi (issue with every new shell!)
    if(TW_COMPILER STREQUAL clang)
   	    set(CMAKE_CXX_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/llvm/bin/clang++")
	    set(CMAKE_C_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/llvm/bin/clang")
        set(CMAKE_CXX_FLAGS "-stdlib=libc++")
        # It seems that libstdc++ handles c++17 special functions correctly and libc++ does not.
        # So if we are using libc++ we have to use our own wrappers, and incur some limitations.
        add_definitions(-DUSE_TWSF)
        set(TW_STD_MODULE_PATH "$ENV{HOMEBREW_PREFIX}/opt/llvm/share/libc++/v1")
        set(TW_STD_MODULE_FILE "std.cppm")
    endif()
    if(TW_COMPILER STREQUAL gcc)
	    set(CMAKE_CXX_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/gcc/bin/g++-15")
	    set(CMAKE_C_COMPILER "$ENV{HOMEBREW_PREFIX}/opt/gcc/bin/gcc-15")
        set(TW_STD_MODULE_PATH "$ENV{HOMEBREW_PREFIX}/include/c++/15/bits")
        set(TW_STD_MODULE_FILE "std.cc")
    endif()
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
    main/Main.cpp
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
    BASE_DIRS base io
)

# Module interfaces go here
target_sources(tw3d PUBLIC
    FILE_SET CXX_MODULES
    BASE_DIRS main base driver io particles fields solver sparc test ${TW_STD_MODULE_PATH}
    FILES
        ${TW_STD_MODULE_PATH}/${TW_STD_MODULE_FILE}
        main/Simulation.cpp
        main/Factory.cpp
        base/base.cpp
        base/logger.cpp
        base/pic_primitives.cpp
        base/tensor.cpp
        base/StaticSpace.cpp
        base/DynSpace.cpp
        base/MetricSpace.cpp
        base/tw_iterator.cpp
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
        test/regions_test.cpp
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
