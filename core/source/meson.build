project('turbowave','cpp','c',
    meson_version: '>=1.1.0',  
    version: '"4.8.0"',
    default_options: ['cpp_std=c++17','buildtype=release']
)

# Meson philosophy is to give the source list explicitly.
# We can always copy and paste from a glob elsewhere.
src_base = ['base/3DFields.cpp', 'base/ComputeTool.cpp', 'base/DiscreteSpace.cpp', 'base/FFT.cpp', 'base/Functions.cpp', 'base/Injection.cpp',
    'base/Injection_Test.cpp', 'base/Main.cpp', 'base/MetricSpace.cpp', 'base/MetricSpace_Test.cpp', 'base/Module.cpp', 'base/Region.cpp', 'base/Simulation.cpp']
src_con = ['concurrency/Tasks.cpp']
src_io = ['io/Diagnostics.cpp', 'io/Input.cpp', 'io/Units.cpp']
src_par = ['particles/Bundle.cpp', 'particles/Mover.cpp', 'particles/Mover_Test.cpp', 'particles/Particles.cpp', 'particles/Particles_Test.cpp',
    'particles/Pusher.cpp', 'particles/Pusher_Test.cpp', 'particles/Slicer.cpp', 'particles/Tiler.cpp', 'particles/Qed.cpp']
src_solv = ['solver/Electrostatic.cpp', 'solver/Elliptic.cpp', 'solver/Elliptic_Test.cpp', 'solver/FieldSolve.cpp', 'solver/Hyperbolic.cpp',
    'solver/LaserSolve.cpp', 'solver/Numerics.cpp', 'solver/Parabolic.cpp', 'solver/Qstate.cpp', 'solver/Quantum.cpp', 'solver/SolidState.cpp']
src_sparc = ['sparc/Chemistry.cpp', 'sparc/Fluid.cpp', 'sparc/Fluid_Test.cpp', 'sparc/Physics.cpp', 'sparc/Physics_Test.cpp']

conf_data = configuration_data()
conf_data.set('TW_VERSION_STRING',meson.project_version())
conf_data.set('BIGENDIAN',get_option('bigendian'))
conf_data.set('VBITS',get_option('vbits'))
conf_data.set('USE_DESKTOP',not get_option('hpc'))
conf_data.set('USE_HPC',get_option('hpc'))
conf_data.set('USE_OPENMP',get_option('omp'))
conf_data.set('USE_OPENCL',get_option('ocl'))
conf_data.set('USE_TWSF',get_option('twsf'))

if get_option('headers')!=''
    incdir = include_directories('base','concurrency','io','particles','solver','sparc',get_option('headers'))
else
    incdir = include_directories('base','concurrency','io','particles','solver','sparc')
endif

deps = []
if get_option('omp')
    deps += [dependency('openmp', language: 'cpp')]
endif
if get_option('ocl')
    deps += [dependency('opencl')]
    conf_data.set('CL_TARGET_OPENCL_VERSION','120')
endif
if get_option('hpc')
    deps += [dependency('mpi')]
else
    src_con += ['concurrency/TW_MPI.cpp']
endif

configure_file(output: 'config.h', configuration: conf_data)
executable(
    'tw3d',
    src_base + src_con + src_io + src_par + src_solv + src_sparc,
    include_directories : incdir,
    dependencies: deps,
    install: true
)
