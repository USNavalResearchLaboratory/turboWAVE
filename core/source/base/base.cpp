module;

#include "tw_includes.h"

export module base;
#ifdef USE_STD_MODULE
	import module std;
#endif

constexpr bool LittleEndian() {
	return std::endian::native == std::endian::little;
}

export namespace tw
{
	typedef double Float;
	typedef int64_t Int;
	typedef uint64_t Uint;
	using Complex = std::complex<tw::Float>;
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
export namespace term
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
	const std::string error("\u001b[31mERROR\u001b[39;49m");
}

export namespace tw
{
	const tw::Int cache_align_bytes = 64;
	const tw::Int vec_align_bytes = VBITS/8; // if not matched to hardware can lead to failures
	const tw::Int max_bundle_size = 16; // must be multiple of vec_align_bytes / sizeof(float)
	const tw::Float small_neg = -1e9*std::numeric_limits<tw::Float>::min();
	const tw::Float small_pos = 1e9*std::numeric_limits<tw::Float>::min();
	const tw::Float big_neg = -1e-9*std::numeric_limits<tw::Float>::max();
	const tw::Float big_pos = 1e-9*std::numeric_limits<tw::Float>::max();
	const tw::Float tiny = std::numeric_limits<tw::Float>::epsilon();
	#ifdef _WIN32
	tw::Float * alloc_aligned_floats(size_t count,size_t alignment) {
		return (tw::Float*)_aligned_malloc(count*sizeof(tw::Float),alignment);
	}
	void free_aligned_floats(tw::Float *ptr) {
		_aligned_free(ptr);
	}
	#else
	tw::Float * alloc_aligned_floats(size_t count,size_t alignment) {
		return (tw::Float*)std::aligned_alloc(alignment,count*sizeof(tw::Float));
	}
	void free_aligned_floats(tw::Float *ptr) {
		free(ptr);
	}
	#endif
}
// Define a strongly typed boolean for better type checking in some situations
export enum class strongbool { yes , no };

export const tw::Float pi = 3.1415926535897932385;
export const tw::Complex ii = tw::Complex(0,1);
// Define some trivial constants useful as binary operands with complex numbers.
// This comes about because std::complex<T> binary operators do not perform automatic conversions.
export const tw::Float one = tw::Float(1.0);
export const tw::Float two = tw::Float(2.0);
export const tw::Float half = tw::Float(0.5);
export const tw::Float root2 = std::sqrt(2);

/////////////////////////////////
//                             //
//      ITEMS FOR TESTING      //
//    (macros in tw_test.h)    //
//                             //
/////////////////////////////////

export struct Testable {
	std::string testName;
};

export void assertFailed(tw::Float actual,tw::Float expected,const std::string& expr,const std::string& file,int line,const std::string& func)
{
	std::ostringstream mess;
	mess << std::endl << term::err << " function " << term::red << func << term::reset_color << std::endl;
	mess << "  Assertion " << actual << " (actual) " << expr << " " << expected << " (expected) failed." << std::endl;
	mess << "  File: " << file <<  " , Line: " << line << std::endl;
	throw tw::FatalError(mess.str());
}

export void assertClose(tw::Float actual,tw::Float expected,tw::Float tol,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (std::fabs(actual-expected)>tol)
		assertFailed(actual,expected,"~",file,line,func);
}

export void assertEqualInt(tw::Int actual,tw::Int expected,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (actual!=expected)
		assertFailed(actual,expected,"==",file,line,func);
}

export void assertGtrEq(tw::Int actual,tw::Int expected,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (actual<expected)
		assertFailed(actual,expected,">=",file,line,func);
}

export void assertLessEq(tw::Int actual,tw::Int expected,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (actual>expected)
		assertFailed(actual,expected,"<=",file,line,func);
}

/////////////////////////////////
//                             //
//  ITEMS FOR OpenMP PROGRAMS  //
//                             //
/////////////////////////////////

export namespace tw
{
	tw::Int GetOMPThreadNum() { return omp_get_thread_num(); }
	tw::Int GetOMPNumThreads() { return omp_get_num_threads(); }
	tw::Int GetOMPMaxThreads() { return omp_get_max_threads(); }
	void GetOMPTaskLoopRange(tw::Int task_id,tw::Int num,tw::Int num_tasks,tw::Int *first,tw::Int *last)
	{
		tw::Int locNum = num / num_tasks;
		*first = task_id * locNum;
		*last = *first + locNum - 1;
		*last += task_id==num_tasks-1 ? num % num_tasks : 0;
	}
}

/////////////////////////////////
//                             //
//  ITEMS FOR OpenCL PROGRAMS  //
//                             //
/////////////////////////////////

#ifdef USE_OPENCL

	export std::string CLBasicTypes()
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
	export std::string CLDiscreteSpace()
	{
		return "typedef struct { tw_Int xDim,yDim,zDim,xN0,xN1,yN0,yN1,zN0,zN1; } tw_DiscreteSpace;";
	}
	export std::string CLMetrics()
	{
		return "typedef struct { tw_Float dt,dx,dy,dz; tw_Float to,xo,yo,zo; tw_Float car,cyl,sph,par; } tw_Metrics;";
	}
	export std::string CLStrip()
	{
		return "typedef struct { tw_Int axis,di,dj,dk; tw_Int stride,dim; } tw_Strip;";
	}
	export std::map<std::string,std::string> CLProtocols()
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
	export std::string CLDefinitions(cl_device_id theDevice)
	{
		std::string definitions_string;
		definitions_string = CLBasicTypes() + '\n';
		definitions_string += CLDiscreteSpace() + '\n';
		definitions_string += CLMetrics() + '\n';
		definitions_string += CLStrip() + '\n';
		return definitions_string;
	}

#endif

////////////////////////
//                    //
//  HELPER FUNCTIONS  //
//                    //
////////////////////////

export long MaxSeed() { return 2147483647; }
export const long ia = 16807;
export const long im = 2147483647;
export const long iq = 127773;
export const long ir = 2836;
export const long ntab = 32;
export const long ndiv = (1+(im-1)/ntab);

export tw::Int GetSeconds() {
    return MPI_Wtime();
}

export bool isPowerOfTwo(tw::Uint x)
{
	return ((x != 0) && !(x & (x-1)));
}

export void complex_multiply_assign(tw::Float& z_re,tw::Float& z_im,const tw::Float& c_re,const tw::Float& c_im)
{
	tw::Float temp = z_re;
	z_re = z_re*c_re - z_im*c_im;
	z_im = temp*c_im + z_im*c_re;
}

export tw::Float sqr(const tw::Float& a)
{
	return a*a;
}

export tw::Float cub(const tw::Float& a)
{
	return a*a*a;
}

export tw::Float quad(const tw::Float& a)
{
	return a*a*a*a;
}

export tw::Float SafeDiv(const tw::Float& numerator,const tw::Float& denominator)
{
	return numerator / (denominator + tw::small_pos);
	//return denominator==0.0?0.0:numerator/denominator;
}

export tw::Float arcsinh(const tw::Float& a)
{
	return std::log(a + sqrt(a*a+1.0));
}

export tw::Float pythag(tw::Float a,tw::Float b)
{
	tw::Float absa,absb;
	absa = std::fabs(a);
	absb = std::fabs(b);
	if (absa > absb) return absa*sqrt(1.0 + sqr(absb/absa));
	else return (absb==0.0 ? 0.0 : absb*sqrt(1.0 + sqr(absa/absb)));
}

export tw::Float QuadraticRoot1(const tw::Float& a,const tw::Float& b,const tw::Float& c)
{
	tw::Float sgn_b = b < 0.0 ? -1.0 : 1.0;
	return (-0.5*(b + sgn_b*sqrt(b*b-4.0*a*c)))/a;
}

export tw::Float QuadraticRoot2(const tw::Float& a,const tw::Float& b,const tw::Float& c)
{
	tw::Float sgn_b = b < 0.0 ? -1.0 : 1.0;
	return c/(-0.5*(b + sgn_b*sqrt(b*b-4.0*a*c)));
}

export int NearestIdx(const std::vector<tw::Float>& vector, tw::Float value)
{
	auto const iter = std::lower_bound(vector.begin(), vector.end(), value);
	if (iter == vector.end())
		return -1;
	return iter - vector.begin();
}

export void ReverseBytes(char *bytes,tw::Int n)
{
	tw::Int i;
	char temp;
	for (i=0;i<n/2;i++)
	{
		temp = bytes[i];
		bytes[i] = bytes[n-1-i];
		bytes[n-1-i] = temp;
	}
}

export void WriteBigEndian(const char *bytes,tw::Int n,tw::Int blockSize,std::ofstream& outFile)
{
	if (!blockSize) blockSize = n;
	char *temp = new char[blockSize];
	for (tw::Int i=0;i<n/blockSize;i++)
	{
		memcpy(temp,&bytes[i*blockSize],blockSize);
		if (LittleEndian())
			ReverseBytes(temp,blockSize);
		outFile.write(temp,blockSize);
	}
	delete [] temp;
}

export void ReadLittleEndian(char *bytes,tw::Int n,tw::Int blockSize,std::ifstream& inFile)
{
	if (!blockSize) blockSize = n;
	char *temp = new char[blockSize];
	for (tw::Int i=0;i<n/blockSize;i++)
	{
		inFile.read(temp,blockSize);
		if (!LittleEndian())
			ReverseBytes(temp,blockSize);
		memcpy(&bytes[i*blockSize],temp,blockSize);
	}
	delete [] temp;
}

export void WriteLittleEndian(const char *bytes,tw::Int n,tw::Int blockSize,std::fstream& outFile)
{
	if (!blockSize) blockSize = n;
	char *temp = new char[blockSize];
	for (tw::Int i=0;i<n/blockSize;i++)
	{
		memcpy(temp,&bytes[i*blockSize],blockSize);
		if (!LittleEndian())
			ReverseBytes(temp,blockSize);
		outFile.write(temp,blockSize);
	}
	delete [] temp;
}

export void BufferLittleEndian(char *bytes,tw::Int n)
{
	char *temp = new char[n];
	memcpy(temp,bytes,n);
	if (!LittleEndian())
		ReverseBytes(temp,n);
	memcpy(bytes,temp,n);
	delete [] temp;
}

export tw::Int MyFloor(tw::Float a)
{
	// floor() on OS X seems to be broken in some cases
	// use different return type than usual
	tw::Int ans = a;
	if (a<0.0) ans--;
	return ans;
}

export tw::Int MyCeil(tw::Float a)
{
	// floor() on OS X seems to be broken in some cases
	// use different return type than usual
	tw::Int ans = a;
	if (a>0.0) ans++;
	return ans;
}

export tw::Int TruncateUnlessClose(tw::Float x)
{
	tw::Int ans;
	ans = tw::Int(x);
	if (x - tw::Float(ans) > 0.99999)
		ans++;
	return ans;
}

export struct UniformDeviate
{
	// From numerical recipes "ran1" routine
	// Multiplicative congruential algorithm with shuffle

	tw::Float am,eps,rnmx;
	long valueNow;
	long iy;
	long iv[ntab];

	UniformDeviate(long seed)
	{
		am = tw::Float(1.0)/tw::Float(im);
		eps = 1.2e-7;
		rnmx = 1.0 - eps;

		long j,k;

		valueNow = seed<0 ? -seed : seed;
		if (valueNow == 0) valueNow = 1;
		for (j=ntab+7;j>=0;j--)
		{
			k = valueNow/iq;
			valueNow = ia*(valueNow-k*iq)-ir*k;
			if (valueNow<0) valueNow += im;
			if (j<ntab) iv[j] = valueNow;
		}
		iy = iv[0];
	}
	tw::Float Next()
	{
		long j,k;
		tw::Float temp;

		k = valueNow/iq;
		valueNow = ia*(valueNow-k*iq)-ir*k;
		if (valueNow<0) valueNow += im;
		j = iy/ndiv;
		iy = iv[j];
		iv[j] = valueNow;
		if ((temp=am*iy)>rnmx)
			return rnmx;
		else
			return temp;
	}
	void WriteCheckpoint(std::ofstream& outFile)
	{
		outFile.write((char *)this,sizeof(UniformDeviate));
	}
	void ReadCheckpoint(std::ifstream& inFile)
	{
		inFile.read((char *)this,sizeof(UniformDeviate));
	}
};

export struct GaussianDeviate
{
	// produce distribution with std dev of 1
	// this gives a density of exp(-x^2/2) (note factor of 2)

	UniformDeviate* ud;
	tw::Float x1,x2;

	GaussianDeviate(long seed)
	{
		ud = new UniformDeviate(seed);
	}
	~GaussianDeviate()
	{
		delete ud;
	}
	tw::Float Next()
	{
		x1 = ud->Next();
		x2 = ud->Next();
		return sqrt(std::fabs(2.0*std::log(0.0001 + x1)))*std::cos(6.28*x2);
	}
	void WriteCheckpoint(std::ofstream& outFile)
	{
		ud->WriteCheckpoint(outFile);
		outFile.write((char *)&x1,sizeof(tw::Float));
		outFile.write((char *)&x2,sizeof(tw::Float));
	}
	void ReadCheckpoint(std::ifstream& inFile)
	{
		ud->ReadCheckpoint(inFile);
		inFile.read((char *)&x1,sizeof(tw::Float));
		inFile.read((char *)&x2,sizeof(tw::Float));
	}
};
