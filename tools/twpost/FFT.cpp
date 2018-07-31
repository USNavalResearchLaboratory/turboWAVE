#include "definitions.h"
#include "fft.h"

int64 CountTrueBits(int64 theBits)
{
	int64 i,numTrue = 0;
	int64 numBits = sizeof(int)*8;
	int64 test = 1;
	
	for (i=0;i<numBits;i++)
	{
		if (test & theBits)
			numTrue++;
		theBits >>= 1;
	}
	return numTrue;
}

void ReverseBits(int64& theBits,int64 numBits)
{
	int64 i;
	int64 ans = 0;
	int64 test = 1;
	int64 currVal = 1;
	currVal <<= numBits - 1;
	
	for (i=0;i<numBits;i++)
	{
		if (test & theBits)
			ans += currVal;
		currVal >>= 1;
		test <<= 1;
	}
	theBits = ans;
}

void ComplexFFT(float *realPart,float *imagPart,int64 numPoints,int64 interval,float inversion)
{
	int64 i,j;
	int64 m,n,numBits;
	
	/*** ARRANGE DATA IN BIT REVERSED ORDER ***/
	
	n = numPoints;
	numBits = 0;
	while (n>1)
	{
		n >>= 1;
		numBits++;
	}
	for (n=0;n<numPoints;n++)
	{
		m = n;
		ReverseBits(m,numBits);
		if (m>n)
		{
			Swap(realPart[m*interval],realPart[n*interval]);
			Swap(imagPart[m*interval],imagPart[n*interval]);
		}
	}
	
	/*** CONSTRUCT TRANSFORM ***/
	
	int64 transforms,points,a,b;
	complx W,FW,recurrenceFactor;
	
	transforms = numPoints;
	points = 1;
	
	while (transforms > 1)
	{		
		recurrenceFactor = exp(ii*pi*inversion/float(points));
		W = complx(1,0);
		for (i=0;i<points;i++)
		{
			for (j=0;j<transforms;j+=2)
			{
				a = interval*(j*points + i);
				b = interval*((j+1)*points + i);
				FW = complx(realPart[b],imagPart[b])*W;
				realPart[b] = realPart[a] - real(FW);
				imagPart[b] = imagPart[a] - imag(FW);
				realPart[a] = realPart[a] + real(FW);
				imagPart[a] = imagPart[a] + imag(FW);
			}
			W *= recurrenceFactor;
		}
		points *= 2;
		transforms /= 2;
	}
	
	if (inversion==-1.0)
	{
		for (i=0;i<numPoints;i++)
		{
			realPart[i*interval] /= float(numPoints);
			imagPart[i*interval] /= float(numPoints);
		}
	}
}

void RealFFT(float *array,int64 num,int64 interval,int64 inversion)
{
	// Computes N/2-1 complex numbers and 2 real numbers
	// The two reals are the lowest and highest frequencies
	// the lowest freq is stored in array[0] and the highest in array[1]
	// Then complex(array[2],array[3]) is the second frequency, etc.
	// The negative frequencies are not calculated but can be obtained from symmetry
	
	int64 i,i1,i2,i3,i4;
	float c1=0.5;
	float c2,wtemp,theta;
	float h1r,h1i,h2r,h2i,wr,wi,wpr,wpi;
	
	theta = 3.141592653589793/float(num/2);
	if (inversion==1)
	{
		c2 = -0.5;
		ComplexFFT(array,&array[interval],num/2,2*interval,1.0);
	}
	else
	{
		c2 = 0.5;
		theta = -theta;
	}
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0 + wpr;
	wi = wpi;
	for (i=2;i<=(num/4);i++)
	{
		i1 = interval*(2*i - 2);
		i2 = interval*(2*i - 1);
		i3 = interval*(num + 2 - 2*i);
		i4 = interval*(num + 3 - 2*i);
		h1r = c1*(array[i1] + array[i3]);
		h1i = c1*(array[i2] - array[i4]);
		h2r = -c2*(array[i2] + array[i4]);
		h2i = c2*(array[i1] - array[i3]);
		array[i1] = h1r + wr*h2r - wi*h2i;
		array[i2] = h1i + wr*h2i + wi*h2r;
		array[i3] = h1r - wr*h2r + wi*h2i;
		array[i4] = -h1i + wr*h2i + wi*h2r;
		wr = (wtemp=wr)*wpr - wi*wpi + wr;
		wi = wi*wpr + wtemp*wpi + wi;
	}
	if (inversion==1)
	{
		array[0] = (h1r=array[0])+array[interval];
		array[interval] = h1r - array[interval];
	}
	else
	{
		array[0] = c1*((h1r=array[0])+array[interval]);
		array[interval] = c1*(h1r-array[interval]);
		ComplexFFT(array,&array[interval],num/2,2*interval,-1.0);
	}
}

void CosineTransform(float* array,int64 num,int64 interval,int64 inversion)
{
	// Ordering is such that derivative is zero midway between array[0] and array[-1]
	// and midway between array[num-1] and array[num]
	
	int64 i;
	float* y;
	float sum,sum1,y1,y2,ytemp;
	float theta,wi,wi1,wpi,wpr,wr,wr1,wtemp;
	
	// Copy data into an array of stride one if necessary
	// This can actually be faster than operating on original array
	
	if (interval>1)
	{
		y = new float[num];
		for (i=0;i<num;i++)
			y[i] = array[i*interval];
	}
	else
		y = array;
	
	wr = 1.0;
	wi = 0.0;
	theta = 0.5*pi/num;
	wr1 = cos(theta);
	wi1 = sin(theta);
	wpr = -2.0*wi1*wi1;
	wpi = sin(2.0*theta);
	
	if (inversion==1)
	{
		for (i=0;i<num/2;i++)
		{
			y1 = 0.5*(y[i] + y[num-i-1]);
			y2 = wi1*(y[i] - y[num-i-1]);
			y[i] = y1 + y2;
			y[num-i-1] = y1 - y2;
			wr1 = (wtemp=wr1)*wpr - wi1*wpi + wr1;
			wi1 = wi1*wpr + wtemp*wpi + wi1;
		}
		RealFFT(y,num,1,inversion);
		for (i=2;i<num;i+=2)
		{
			wr = (wtemp=wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
			y1 = y[i]*wr - y[i+1]*wi;
			y2 = y[i+1]*wr + y[i]*wi;
			y[i] = y1;
			y[i+1] = y2;
		}
		sum = 0.5*y[1];
		for (i=num-1;i>=1;i-=2)
		{
			sum1 = sum;
			sum += y[i];
			y[i] = sum1;
		}
	}
	else
	{
		ytemp = y[num-1];
		for (i=num-1;i>=3;i-=2)
			y[i] = y[i-2] - y[i];
		y[1] = 2.0*ytemp;
		for (i=2;i<num;i+=2)
		{
			wr = (wtemp=wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
			y1 = y[i]*wr + y[i+1]*wi;
			y2 = y[i+1]*wr - y[i]*wi;
			y[i] = y1;
			y[i+1] = y2;
		}
		RealFFT(y,num,1,inversion);
		for (i=0;i<num/2;i++)
		{
			y1 = y[i] + y[num-i-1];
			y2 = (0.5/wi1)*(y[i] - y[num-i-1]);
			y[i] = 0.5*(y1 + y2);
			y[num-i-1] = 0.5*(y1 - y2);
			wr1 = (wtemp=wr1)*wpr - wi1*wpi + wr1;
			wi1 = wi1*wpr + wtemp*wpi + wi1;
		}
	}
	
	if (interval>1)
	{
		for (i=0;i<num;i++)
			array[i*interval] = y[i];
		delete [] y;
	}
}

void SineTransform(float* array,int64 num,int64 interval,int64 inversion)
{
	// Ordering is such that array[0] = 0
	// array[num] is implicitly zero but is not computed
	
	int64 j;
	float sum,y1,y2;
	float theta,wi,wpi,wpr,wr,wtemp;
	float* y;

	// Copy data into an array of stride one if necessary
	// This can actually be faster than operating on original array
	
	if (interval>1)
	{
		y = new float[num];
		for (j=0;j<num;j++)
			y[j] = array[j*interval];
	}
	else
		y = array;
	
	theta = 3.141592653589793/float(num);
	wr = 1.0;
	wi = 0.0;
	wpr = -2.0*pow(sin(0.5*theta),2);
	wpi = sin(theta);
	y[0] = 0.0;
	for (j=1;j<=num/2;j++)
	{
		wtemp = wr;
		wr = wr*wpr - wi*wpi + wr;
		wi = wi*wpr + wtemp*wpi + wi;
		y1 = wi*(y[j] + y[num - j]);
		y2 = 0.5*(y[j] - y[num - j]);
		y[j] = y1 + y2;
		y[num - j] = y1 - y2;
	}
	RealFFT(y,num,1,1);
	sum = 0.0;
	y[0] = 0.5*y[0];
	y[1] = 0.0;
	for (j=1;j<=num-1;j+=2)
	{
		sum += y[j - 1];
		y[j - 1] = y[j];
		y[j] = sum;
	}
	
	if (inversion==-1)
		for (j=0;j<num;j++)
			y[j] *= 2.0/float(num);
			
	if (interval>1)
	{
		for (j=0;j<num;j++)
			array[j*interval] = y[j];
		delete [] y;
	}
}

void ShiftForSineTransform(float* array,int64 num,int64 interval)
{
	int64 j;
	for (j=num-1;j>=1;j--)
	{
		array[j*interval] = 0.5*(array[j*interval] + array[(j-1)*interval]);
	}
	array[0] = 0.0;
}

void ShiftForInverseSineTransform(float* array,int64 num,int64 interval)
{
	int64 j;
	for (j=0;j<num-1;j++)
	{
		array[j*interval] = 0.5*(array[j*interval] + array[(j+1)*interval]);
	}
	array[num-1] *= 0.5;
}

void ToggleOrdering(complx *array,int64 num,int64 interval)
{
	int64 i;
	
	for (i=0;i<num/2;i++)
			Swap(array[i*interval],array[(i+num/2)*interval]);
}

