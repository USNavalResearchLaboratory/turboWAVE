// This is used to bring in whatever we cannot bring in with the module system.

#include <config.h>

#ifdef USE_STD_MODULE
	#include <malloc.h>
	#include <assert.h>
	#include <memory.h>
#else
	#include <ranges>
	#include <bit>
	#include <cassert>
	#include <cmath>
	#include <complex>
	#include <iostream>
	#include <fstream>
	#include <sstream>
	#include <ctime>
	#include <vector>
	#include <valarray>
	#include <map>
	#include <limits>
	#include <exception>
	#include <string>
	#include <cctype>
	#include <format>
	#include <print>
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
