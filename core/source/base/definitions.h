#ifndef _did_definitions
#define TW_VERSION_STRING "4.7.0a2"

/////////////////////////////////////////////////////////////////////
//                                                                 //
// This file contains definitions of basic data types and provides //
// a generalized interface for platform specific features.         //
// It also contains constants affecting optimization               //
//                                                                 //
// Following are now defined on the compile line in the makefile:  //
// USE_DESKTOP : if defined compile for a desktop system           //
// USE_HPC : if defined compile for a generic MPI cluster          //
// USE_OPENCL : if defined executable uses OpenCL kernels          //
// USE_OPENMP : if defined executable uses OpenMP threads          //
// VBITS : used to align data for vector units                     //
//                                                                 //
// Coordinate system is fully controllable from input file.        //
//                                                                 //
/////////////////////////////////////////////////////////////////////

//#include <stdint.h>
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
#include <algorithm>
#include <limits>
#include <exception>
#include <string>
#include <cctype>
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace tw
{
	typedef double Float;
	typedef int32_t Int;
	typedef uint32_t Uint;
	typedef std::complex<tw::Float> Complex;
	static const tw::Int cache_align_bytes = 64;
	static const tw::Int vec_align_bytes = VBITS/8; // if not matched to hardware can lead to failures
	static const tw::Int max_bundle_size = 16; // must be multiple of vec_align_bytes / sizeof(float)
	static const tw::Float small_neg = -1e9*std::numeric_limits<tw::Float>::min();
	static const tw::Float small_pos = 1e9*std::numeric_limits<tw::Float>::min();
	static const tw::Float big_neg = -1e-9*std::numeric_limits<tw::Float>::max();
	static const tw::Float big_pos = 1e-9*std::numeric_limits<tw::Float>::max();
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
// Define a strongly typed boolean for better type checking in some situations
enum class strongbool { yes , no };

static const tw::Float pi = 3.1415926535897932385;
static const tw::Complex ii = tw::Complex(0,1);
// Define some trivial constants useful as binary operands with complex numbers.
// This comes about because std::complex<T> binary operators do not perform automatic conversions.
static const tw::Float one = tw::Float(1.0);
static const tw::Float two = tw::Float(2.0);
static const tw::Float half = tw::Float(0.5);
static const tw::Float root2 = std::sqrt(2);

// ANSI terminal codes
namespace term
{
	static const std::string ok("\u2713");
	static const std::string err("\u2717");
	static const std::string red("\u001b[31m");
	static const std::string green("\u001b[32m");
	static const std::string blue("\u001b[94m");
	static const std::string yellow("\u001b[33m");
	static const std::string cyan("\u001b[96m");
	static const std::string reset_color("\u001b[39;49m");
	static const std::string reset_all("\u001b[0m");
	static const std::string bold("\u001b[1m");
}

/////////////////////////////////
//                             //
//  ITEMS FOR OpenMP PROGRAMS  //
//                             //
/////////////////////////////////

#ifdef USE_OPENMP
namespace tw
{
	inline tw::Int GetOMPThreadNum() { return omp_get_thread_num(); }
	inline tw::Int GetOMPNumThreads() { return omp_get_num_threads(); }
	inline tw::Int GetOMPMaxThreads() { return omp_get_max_threads(); }
	inline void GetOMPTaskLoopRange(tw::Int task_id,tw::Int num,tw::Int num_tasks,tw::Int *first,tw::Int *last)
	{
		tw::Int locNum = num / num_tasks;
		*first = task_id * locNum;
		*last = *first + locNum - 1;
		*last += task_id==num_tasks-1 ? num % num_tasks : 0;
	}
}
#else
namespace tw
{
	inline tw::Int GetOMPThreadNum() { return 0; }
	inline tw::Int GetOMPNumThreads() { return 1; }
	inline tw::Int GetOMPMaxThreads() { return 1; }
	inline void GetOMPTaskLoopRange(tw::Int task_id,tw::Int num,tw::Int num_tasks,tw::Int *first,tw::Int *last)
	{
		// If no OpenMP, still split tasks for serial execution, if requested.
		// To suppress loop splitting in serial code, must pass in task_id=0 and num_tasks=1.
		tw::Int locNum = num / num_tasks;
		*first = task_id * locNum;
		*last = *first + locNum - 1;
		*last += task_id==num_tasks-1 ? num % num_tasks : 0;
	}
}
#endif

/////////////////////////////////
//                             //
//  ITEMS FOR OpenCL PROGRAMS  //
//                             //
/////////////////////////////////

#ifdef USE_OPENCL

	#ifdef __APPLE__
		#include <OpenCL/cl.h>
	#else
		#include <CL/cl.h>
	#endif

	inline std::string CLDoublePragma(cl_device_id theDevice)
	{
		char buff[4096];
		size_t buffSize;
		std::string extensionList;
		if (sizeof(tw::Float)==4)
			return "";
		else
		{
			clGetDeviceInfo(theDevice,CL_DEVICE_EXTENSIONS,sizeof(buff),buff,&buffSize);
			extensionList = std::string(buff,buffSize-1);
			extensionList = buff;
			if (extensionList.find("cl_amd_fp64")!=std::string::npos)
				return "#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n";
			if (extensionList.find("cl_khr_fp64")!=std::string::npos)
				return "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
			return "";
		}
	}
	inline std::string CLBasicTypes()
	{
		std::string types_string;
		if (sizeof(tw::Float)==4)
			types_string = "typedef float tw_Float;typedef float2 tw_Complex;typedef float4 tw_vec4; ";
		else
			types_string = "typedef double tw_Float;typedef double2 tw_Complex;typedef double4 tw_vec4; ";
		if (sizeof(tw::Int)==4)
			types_string += "typedef int tw_Int;typedef int4 tw_cell_id;";
		else
			types_string += "typedef long tw_Int;typedef int4 tw_cell_id;";
		return types_string;
	}
	inline std::string CLDiscreteSpace()
	{
		return "typedef struct { tw_Int xDim,yDim,zDim,xN0,xN1,yN0,yN1,zN0,zN1; } tw_DiscreteSpace;";
	}
	inline std::string CLMetrics()
	{
		return "typedef struct { tw_Float dt,dx,dy,dz; tw_Float to,xo,yo,zo; tw_Float car,cyl,sph,par; } tw_Metrics;";
	}
	inline std::string CLStrip()
	{
		return "typedef struct { tw_Int axis,di,dj,dk; tw_Int stride,dim; } tw_Strip;";
	}
	inline std::map<std::string,std::string> CLProtocols()
	{
		return {
			{"LOCAL_PROTOCOL_MACRO","const tw_Int ig = get_global_id(0); \
			const tw_Int jg = get_global_id(1); \
			const tw_Int kg = get_global_id(2); \
			const tw_Int nx = get_global_size(0) + 2*get_global_offset(0); \
			const tw_Int ny = get_global_size(1) + 2*get_global_offset(1); \
			const tw_Int nz = get_global_size(2) + 2*get_global_offset(2); \
			const tw_Int xs = (nx!=1); \
			const tw_Int ys = (ny!=1)*nx; \
			const tw_Int zs = (nz!=1)*nx*ny; \
			const tw_Int cs = nx*ny*nz; \
			const tw_Int n = ig*xs + jg*ys + kg*zs; \
			tw_cell_id cell = (tw_cell_id)(0,ig,jg,kg); \
			tw_Metrics met; \
			GetMetrics((tw_Float*)&met,met_g);"},

			{"POINT_PROTOCOL_MACRO","const tw_Int ig = get_global_id(0); \
			const tw_Int jg = get_global_id(1); \
			const tw_Int kg = get_global_id(2); \
			const tw_Int nx = get_global_size(0); \
			const tw_Int ny = get_global_size(1); \
			const tw_Int nz = get_global_size(2); \
			const tw_Int xs = 1; \
			const tw_Int ys = nx; \
			const tw_Int zs = nx*ny; \
			const tw_Int cs = nx*ny*nz; \
			const tw_Int n = ig*xs + jg*ys + kg*zs; \
			tw_cell_id cell = (tw_cell_id)(0,ig,jg,kg);"},

			{"POINT_PROTOCOL_METRIC_MACRO","const tw_Int ig = get_global_id(0); \
			const tw_Int jg = get_global_id(1); \
			const tw_Int kg = get_global_id(2); \
			const tw_Int nx = get_global_size(0); \
			const tw_Int ny = get_global_size(1); \
			const tw_Int nz = get_global_size(2); \
			const tw_Int xs = 1; \
			const tw_Int ys = nx; \
			const tw_Int zs = nx*ny; \
			const tw_Int cs = nx*ny*nz; \
			const tw_Int n = ig*xs + jg*ys + kg*zs; \
			tw_cell_id cell = (tw_cell_id)(0,ig,jg,kg); \
			tw_Metrics met; \
			GetMetrics((tw_Float*)&met,met_g);"}
		};
	}
	inline std::string CLDefinitions(cl_device_id theDevice)
	{
		std::string definitions_string;
		definitions_string = CLDoublePragma(theDevice) + '\n';
		definitions_string += CLBasicTypes() + '\n';
		definitions_string += CLDiscreteSpace() + '\n';
		definitions_string += CLMetrics() + '\n';
		definitions_string += CLStrip() + '\n';
		return definitions_string;
	}

#endif

#ifdef BIGENDIAN
	inline bool LittleEndian() { return false; }
#else
	inline bool LittleEndian() { return true; }
#endif


//////////////////////////////
//                          //
//    DESKTOP DEFINITIONS   //
//   (Mac, Windows, Linux)  //
//                          //
//////////////////////////////


#ifdef USE_DESKTOP

#define USE_TW_MPI

inline tw::Int GetSeconds()
{
	return tw::Int(time(NULL));
}

#include <thread>
#include <mutex>
#include <condition_variable>
#include "tw_thread_stl.h"
#include "tw_mpi.h"

static const bool parallelFileSystem = false;

#endif



//////////////////////
//                  //
// CRAY DEFINITIONS //
//                  //
//////////////////////


#ifdef USE_HPC

#include "mpi.h"

inline tw::Int GetSeconds()
{
	return tw::Int(MPI_Wtime());
}

static const bool parallelFileSystem = true;

#endif


#define _did_definitions

#endif
