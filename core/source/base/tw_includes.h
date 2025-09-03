// ubiquitous things that are not provided as modules

#include <config.h>

#ifdef USE_STD_MODULE
	#include <assert.h>
#else
	#include <cassert>
	#include <cmath>
	#include <complex>
	#include <iostream>
	#include <vector>
	#include <valarray>
	#include <map>
	#include <string>
	#include <format>
	#include <print>
	#include <memory>
#endif

#include <omp.h>
#include <mpi.h>
#ifdef USE_OPENCL
	#ifdef __APPLE__
		#include <OpenCL/cl.h>
	#else
		#include <CL/cl.h>
	#endif
#endif
