
tw::Int CountTrueBits(tw::Int theBits);

void ReverseBits(tw::Int& theBits,tw::Int numBits);

void ComplexFFT(tw::Float *realPart,tw::Float *imagPart,tw::Int numPoints,tw::Int interval,tw::Float inversion);

void RealFFT(tw::Float *array,tw::Int num,tw::Int interval,tw::Int inversion);

void CosineTransform(tw::Float *array,tw::Int num,tw::Int interval,tw::Int inversion);

void SineTransform(tw::Float *array,tw::Int num,tw::Int interval,tw::Int inversion);

void ShiftForSineTransform(tw::Float *array,tw::Int num,tw::Int interval);

void ShiftForInverseSineTransform(tw::Float *array,tw::Int num,tw::Int interval);

// eigenvalues of the laplacian

inline tw::Float eigenvalue_CFFT(tw::Int i,tw::Int num,tw::Float dxi)
{
	if (num>1)
		return 2.0*dxi*dxi*(cos(tw::Float(2.0*pi*i)/tw::Float(num)) - 1.0);
	else
		return 0.0;
}

inline tw::Float eigenvalue_RFFT(tw::Int i,tw::Int num,tw::Float dxi)
{
	if (num>1)
	{
		if (i!=1)
			return 2.0*dxi*dxi*(cos(tw::Float(2.0*pi*(i/2))/tw::Float(num)) - 1.0);
		else
			return -4.0*dxi*dxi;
	}
	else
		return 0.0;
}

inline tw::Float eigenvalue_FST(tw::Int i,tw::Int num,tw::Float dxi)
{
	if (num>1)
		return 2.0*dxi*dxi*(cos(tw::Float(pi*i)/tw::Float(num)) - 1.0);
	else
		return 0.0;
}

inline tw::Float eigenvalue_FCT(tw::Int i,tw::Int num,tw::Float dxi)
{
	if (num>1)
		return 2.0*dxi*dxi*(cos(tw::Float(pi*i)/tw::Float(num)) - 1.0);
	else
		return 0.0;
}

// utilities

void ToggleOrdering(tw::Complex *array,tw::Int num,tw::Int interval);

inline void Swap(tw::Float& a,tw::Float& b)
{
	tw::Float temp;
	temp = a;
	a = b;
	b = temp;
}

inline void Swap(tw::Complex& a,tw::Complex& b)
{
	tw::Complex temp;
	temp = a;
	a = b;
	b = temp;
}
