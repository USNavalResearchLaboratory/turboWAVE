int64 CountTrueBits(int64 theBits);

void ReverseBits(int64& theBits,int64 numBits);

void ComplexFFT(float *realPart,float *imagPart,int64 numPoints,int64 interval,float inversion);

void RealFFT(float *array,int64 num,int64 interval,int64 inversion);

void CosineTransform(float *array,int64 num,int64 interval,int64 inversion);

void SineTransform(float *array,int64 num,int64 interval,int64 inversion);

void ShiftForSineTransform(float *array,int64 num,int64 interval);

void ShiftForInverseSineTransform(float *array,int64 num,int64 interval);


inline float frequency_RFFT(int64 i,int64 num,float dxi)
{
	return 6.28*dxi*(i/2)/float(num);
}

// eigenvalues of the laplacian

inline float eigenvalue_CFFT(int64 i,int64 num,float dxi)
{
	if (num>1)
		return 2.0*dxi*dxi*(cos(float(2.0*pi*i)/float(num)) - 1.0);
	else
		return 0.0;
}

inline float eigenvalue_RFFT(int64 i,int64 num,float dxi)
{
	if (num>1)
	{
		if (i!=1)
			return 2.0*dxi*dxi*(cos(float(2.0*pi*(i/2))/float(num)) - 1.0);
		else
			return -4.0*dxi*dxi;
	}
	else
		return 0.0;
}

inline float eigenvalue_FST(int64 i,int64 num,float dxi)
{
	if (num>1)
		return 2.0*dxi*dxi*(cos(float(pi*i)/float(num)) - 1.0);
	else
		return 0.0;
}

inline float eigenvalue_FCT(int64 i,int64 num,float dxi)
{
	if (num>1)
		return 2.0*dxi*dxi*(cos(float(pi*i)/float(num)) - 1.0);
	else
		return 0.0;
}

// utilities

void ToggleOrdering(complx *array,int64 num,int64 interval);

inline void Swap(float& a,float& b)
{
	float temp;
	temp = a;
	a = b;
	b = temp;
}

inline void Swap(complx& a,complx& b)
{
	complx temp;
	temp = a;
	a = b;
	b = temp;
}
