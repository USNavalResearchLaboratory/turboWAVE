static long MaxSeed() { return 2147483647; }
static const long ia = 16807;
static const long im = 2147483647;
static const long iq = 127773;
static const long ir = 2836;
static const long ntab = 32;
static const long ndiv = (1+(im-1)/ntab);

#define ASSERT_EQ(actual,expected) assertEqualInt(actual,expected,__FILE__,__LINE__,__func__,testName)
#define ASSERT_NEAR(actual,expected,tol) assertClose(actual,expected,tol,__FILE__,__LINE__,__func__,testName)
#define ASSERT_GTREQ(actual,expected) assertGtrEq(actual,expected,__FILE__,__LINE__,__func__,testName)
#define ASSERT_LESSEQ(actual,expected) assertLessEq(actual,expected,__FILE__,__LINE__,__func__,testName)
#define REGISTER_TEST() testName = __func__

inline void assertFailed(tw::Float actual,tw::Float expected,const std::string& expr,const std::string& file,int line,const std::string& func)
{
	std::ostringstream mess;
	mess << std::endl << term::err << " function " << term::red << func << term::reset_color << std::endl;
	mess << "  Assertion " << actual << " (actual) " << expr << " " << expected << " (expected) failed." << std::endl;
	mess << "  File: " << file <<  " , Line: " << line << std::endl;
	throw tw::FatalError(mess.str());
}

inline void assertClose(tw::Float actual,tw::Float expected,tw::Float tol,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (fabs(actual-expected)>tol)
		assertFailed(actual,expected,"~",file,line,func);
}

inline void assertEqualInt(tw::Int actual,tw::Int expected,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (actual!=expected)
		assertFailed(actual,expected,"==",file,line,func);
}

inline void assertGtrEq(tw::Int actual,tw::Int expected,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (actual<expected)
		assertFailed(actual,expected,">=",file,line,func);
}

inline void assertLessEq(tw::Int actual,tw::Int expected,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (actual>expected)
		assertFailed(actual,expected,"<=",file,line,func);
}

inline bool isPowerOfTwo(tw::Uint x)
{
	return ((x != 0) && !(x & (x-1)));
}

inline void complex_multiply_assign(tw::Float& z_re,tw::Float& z_im,const tw::Float& c_re,const tw::Float& c_im)
{
	tw::Float temp = z_re;
	z_re = z_re*c_re - z_im*c_im;
	z_im = temp*c_im + z_im*c_re;
}

inline tw::Float sqr(const tw::Float& a)
{
	return a*a;
}

inline tw::Float cub(const tw::Float& a)
{
	return a*a*a;
}

inline tw::Float SafeDiv(const tw::Float& numerator,const tw::Float& denominator)
{
	return numerator / (denominator + tw::small_pos);
	//return denominator==0.0?0.0:numerator/denominator;
}

inline tw::Float arcsinh(const tw::Float& a)
{
	return log(a + sqrt(a*a+1.0));
}

inline tw::Float pythag(tw::Float a,tw::Float b)
{
	tw::Float absa,absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa*sqrt(1.0 + sqr(absb/absa));
	else return (absb==0.0 ? 0.0 : absb*sqrt(1.0 + sqr(absa/absb)));
}

inline tw::Float QuadraticRoot1(const tw::Float& a,const tw::Float& b,const tw::Float& c)
{
	tw::Float sgn_b = b < 0.0 ? -1.0 : 1.0;
	return (-0.5*(b + sgn_b*sqrt(b*b-4.0*a*c)))/a;
}

inline tw::Float QuadraticRoot2(const tw::Float& a,const tw::Float& b,const tw::Float& c)
{
	tw::Float sgn_b = b < 0.0 ? -1.0 : 1.0;
	return c/(-0.5*(b + sgn_b*sqrt(b*b-4.0*a*c)));
}

inline int NearestIdx(const std::vector<tw::Float>& vector, tw::Float value)
{
	auto const iter = std::lower_bound(vector.begin(), vector.end(), value);
    if (iter == vector.end())
    	return -1;
    return iter - vector.begin();
}

inline void ReverseBytes(char *bytes,tw::Int n)
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

inline void WriteBigEndian(const char *bytes,tw::Int n,tw::Int blockSize,std::ofstream& outFile)
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

inline void ReadLittleEndian(char *bytes,tw::Int n,tw::Int blockSize,std::ifstream& inFile)
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

inline void WriteLittleEndian(const char *bytes,tw::Int n,tw::Int blockSize,std::fstream& outFile)
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

inline void BufferLittleEndian(char *bytes,tw::Int n)
{
	char *temp = new char[n];
	memcpy(temp,bytes,n);
	if (!LittleEndian())
		ReverseBytes(temp,n);
	memcpy(bytes,temp,n);
	delete [] temp;
}

inline tw::Int MyFloor(tw::Float a)
{
	// floor() on OS X seems to be broken in some cases
	// use different return type than usual
	tw::Int ans = a;
	if (a<0.0) ans--;
	return ans;
}

inline tw::Int MyCeil(tw::Float a)
{
	// floor() on OS X seems to be broken in some cases
	// use different return type than usual
	tw::Int ans = a;
	if (a>0.0) ans++;
	return ans;
}

inline tw::Int TruncateUnlessClose(tw::Float x)
{
	tw::Int ans;
	ans = tw::Int(x);
	if (x - tw::Float(ans) > 0.99999)
		ans++;
	return ans;
}

struct UniformDeviate
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

struct GaussianDeviate
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
		return sqrt(fabs(2.0*log(0.0001 + x1)))*cos(6.28*x2);
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

// namespace tw
// {
// 	template <class T>
// 	class valarray
// 	{
// 		std::size_t num;
// 		T* xbuff;
// 		T* data;
//
// 		public:
// 		valarray()
// 		{
// 			num = 0;
// 			xbuff = NULL;
// 			data = NULL;
// 		}
// 		valarray(std::size_t n)
// 		{
// 			num = n;
// 			xbuff = (T*)std::malloc(num*sizeof(T));
// 			data = xbuff;
// 		}
// 		valarray(tw::Int AB,std::size_t n)
// 		{
// 			// Dangerous? see notes on resize below.
// 			num = n;
// 			std::size_t ext_bytes = num*sizeof(T)+AB;
// 			const std::size_t req_bytes = num*sizeof(T);
// 			xbuff = (T*)std::malloc(ext_bytes);
// 			void* temp = (void*)xbuff;
// 			data = (T*)std::align(AB,req_bytes,temp,ext_bytes);
// 		}
// 		~valarray()
// 		{
// 			if (xbuff!=NULL)
// 				std::free(xbuff);
// 		}
// 		void resize(std::size_t n)
// 		{
// 			if (xbuff!=NULL)
// 				std::free(xbuff);
// 			num = n;
// 			xbuff = (T*)std::malloc(num*sizeof(T));
// 			data = xbuff;
// 		}
// 		void resize(tw::Int AB,std::size_t n)
// 		{
// 			// Dangerous way to do this.
// 			// std::valarray.resize can be called the same way with
// 			// different interpretation of arguments.
// 			if (xbuff!=NULL)
// 				std::free(xbuff);
// 			num = n;
// 			std::size_t ext_bytes = num*sizeof(T)+AB;
// 			const std::size_t req_bytes = num*sizeof(T);
// 			xbuff = (T*)std::malloc(ext_bytes);
// 			void* temp = (void*)xbuff;
// 			data = (T*)std::align(AB,req_bytes,temp,ext_bytes);
// 		}
// 		std::size_t size() const
// 		{
// 			return num;
// 		}
// 		T& operator [] (const tw::Int& x)
// 		{
// 			return data[x];
// 		}
// 		T operator [] (const tw::Int& x) const
// 		{
// 			return data[x];
// 		}
// 		tw::valarray<T>& operator = (const tw::valarray<T>& A)
// 		{
// 			if (num!=A.size())
// 				resize(A.size());
// 			for (tw::Int i=0;i<num;i++)
// 				data[i] = A[i];
// 			return *this;
// 		}
// 		tw::valarray<T>& operator += (const tw::valarray<T>& A)
// 		{
// 			for (tw::Int i=0;i<num;i++)
// 				data[i] += A[i];
// 			return *this;
// 		}
// 		tw::valarray<T>& operator -= (const tw::valarray<T>& A)
// 		{
// 			for (tw::Int i=0;i<num;i++)
// 				data[i] -= A[i];
// 			return *this;
// 		}
// 		tw::valarray<T>& operator *= (const tw::valarray<T>& A)
// 		{
// 			for (tw::Int i=0;i<num;i++)
// 				data[i] *= A[i];
// 			return *this;
// 		}
// 		tw::valarray<T>& operator /= (const tw::valarray<T>& A)
// 		{
// 			for (tw::Int i=0;i<num;i++)
// 				data[i] /= A[i];
// 			return *this;
// 		}
// 		tw::valarray<T>& operator = (const T& A)
// 		{
// 			for (tw::Int i=0;i<num;i++)
// 				data[i] = A;
// 			return *this;
// 		}
// 		tw::valarray<T>& operator += (const T& A)
// 		{
// 			for (tw::Int i=0;i<num;i++)
// 				data[i] += A;
// 			return *this;
// 		}
// 		tw::valarray<T>& operator -= (const T& A)
// 		{
// 			for (tw::Int i=0;i<num;i++)
// 				data[i] -= A;
// 			return *this;
// 		}
// 		tw::valarray<T>& operator *= (const T& A)
// 		{
// 			for (tw::Int i=0;i<num;i++)
// 				data[i] *= A;
// 			return *this;
// 		}
// 		tw::valarray<T>& operator /= (const T& A)
// 		{
// 			for (tw::Int i=0;i<num;i++)
// 				data[i] /= A;
// 			return *this;
// 		}
// 	};
// }
