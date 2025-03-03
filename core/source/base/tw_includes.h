// This header must be separate from the `base` module because the following
// libraries are not available as modules with sufficient reliability:
// * C++ standard library
// * OpenMP
// * OpenCL
// * MPI

#include <config.h>
#include <memory.h>
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
#include <omp.h>
#include <mpi.h>
#ifdef USE_OPENCL
	#ifdef __APPLE__
		#include <OpenCL/cl.h>
	#else
		#include <CL/cl.h>
	#endif
#endif
#ifdef BIGENDIAN
	constexpr bool LittleEndian() { return false; }
#else
	constexpr bool LittleEndian() { return true; }
#endif

namespace tw
{
	typedef double Float;
	typedef int64_t Int;
	typedef uint64_t Uint;
	typedef std::complex<tw::Float> Complex;
	class FatalError : public std::exception
	{
		char messg[256];
		public:
		FatalError(const std::string& s)
		{
			s.copy(messg,sizeof(messg));
			if (s.size()<sizeof(messg))
				messg[s.size()] = 0;
		}
		virtual const char* what() const throw()
		{
			return messg;
		}
	};
}

// ANSI terminal codes
namespace term
{
	const std::string ok("\u2713");
	const std::string err("\u2717");
	const std::string red("\u001b[31m");
	const std::string green("\u001b[32m");
	const std::string blue("\u001b[94m");
	const std::string yellow("\u001b[33m");
	const std::string cyan("\u001b[96m");
	const std::string reset_color("\u001b[39;49m");
	const std::string reset_all("\u001b[0m");
	const std::string bold("\u001b[1m");
	const std::string warning("\u001b[33mWARNING\u001b[39;49m");
}
