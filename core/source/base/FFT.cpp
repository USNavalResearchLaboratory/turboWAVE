module;

#include "tw_includes.h"

export module fft;
import base;

export namespace fft {

	tw::Int CountTrueBits(tw::Int theBits) {
		tw::Int i,numTrue = 0;
		tw::Int numBits = sizeof(tw::Int)*8;
		tw::Int test = 1;

		for (i=0;i<numBits;i++)
		{
			if (test & theBits)
				numTrue++;
			theBits >>= 1;
		}
		return numTrue;
	}

	void ReverseBits(tw::Int& theBits,tw::Int numBits) {
		tw::Int i;
		tw::Int ans = 0;
		tw::Int test = 1;
		tw::Int currVal = 1;
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

	/// `inversion == 1.0` forward, `inversion == -1.0` reverse
	void ComplexFFT(tw::Float *realPart,tw::Float *imagPart,tw::Int numPoints,tw::Int interval,tw::Float inversion) {
		/*** ARRANGE DATA IN BIT REVERSED ORDER ***/

		tw::Int n = numPoints;
		tw::Int numBits = 0;
		while (n>1)
		{
			n >>= 1;
			numBits++;
		}
		for (tw::Int n=0; n<numPoints; n++)
		{
			tw::Int m = n;
			ReverseBits(m,numBits);
			if (m>n)
			{
				std::swap(realPart[m*interval],realPart[n*interval]);
				std::swap(imagPart[m*interval],imagPart[n*interval]);
			}
		}

		/*** CONSTRUCT TRANSFORM ***/

		tw::Int transforms,points,a,b;
		tw::Complex W,FW,recurrenceFactor;

		transforms = numPoints;
		points = 1;

		while (transforms > 1)
		{		
			recurrenceFactor = std::exp(ii*pi*inversion/tw::Float(points));
			W = tw::Complex(1,0);
			for (auto i=0;i<points;i++)
			{
				for (auto j=0;j<transforms;j+=2)
				{
					a = interval*(j*points + i);
					b = interval*((j+1)*points + i);
					FW = tw::Complex(realPart[b],imagPart[b])*W;
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
			for (auto i=0;i<numPoints;i++)
			{
				realPart[i*interval] /= tw::Float(numPoints);
				imagPart[i*interval] /= tw::Float(numPoints);
			}
		}
	}

	/// `inversion == 1` forward, `inversion == -1` reverse
	void RealFFT(tw::Float *array,tw::Int num,tw::Int interval,tw::Int inversion) {
		// Computes N/2-1 complex numbers and 2 real numbers
		// The two reals are the lowest and highest frequencies
		// the lowest freq is stored in array[0] and the highest in array[1]
		// Then complex(array[2],array[3]) is the second frequency, etc.
		// The negative frequencies are not calculated but can be obtained from symmetry

		tw::Int i,i1,i2,i3,i4;
		tw::Float c1=0.5;
		tw::Float c2,wtemp,theta;
		tw::Float h1r,h1i,h2r,h2i,wr,wi,wpr,wpi;

		theta = 3.141592653589793/tw::Float(num/2);
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
		wtemp = std::sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = std::sin(theta);
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

	void CosineTransform(tw::Float *array,tw::Int num,tw::Int interval,tw::Int inversion) {
		// Ordering is such that derivative is zero midway between array[0] and array[-1]
		// and midway between array[num-1] and array[num]

		tw::Int i;
		tw::Float sum,sum1,y1,y2,ytemp;
		tw::Float theta,wi,wi1,wpi,wpr,wr,wr1,wtemp;
		std::valarray<tw::Float> y(num);

		// Copy data into a valarray of stride one

		for (i=0;i<num;i++)
			y[i] = array[i*interval];

		wr = 1.0;
		wi = 0.0;
		theta = 0.5*pi/num;
		wr1 = std::cos(theta);
		wi1 = std::sin(theta);
		wpr = -2.0*wi1*wi1;
		wpi = std::sin(2.0*theta);

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
			RealFFT(&y[0],num,1,inversion);
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
			RealFFT(&y[0],num,1,inversion);
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

		for (i=0;i<num;i++)
			array[i*interval] = y[i];
	}

	void SineTransform(tw::Float *array,tw::Int num,tw::Int interval,tw::Int inversion) {
		// Ordering is such that array[0] = 0
		// array[num] is implicitly zero but is not computed

		tw::Int j;
		tw::Float sum,y1,y2;
		tw::Float theta,wi,wpi,wpr,wr,wtemp;
		std::valarray<tw::Float> y(num);

		// Copy data into valarray of stride one

		for (j=0;j<num;j++)
			y[j] = array[j*interval];

		theta = 3.141592653589793/tw::Float(num);
		wr = 1.0;
		wi = 0.0;
		wpr = -2.0*std::pow(std::sin(0.5*theta),2.0);
		wpi = std::sin(theta);
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
		RealFFT(&y[0],num,1,1);
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
				y[j] *= 2.0/tw::Float(num);

		for (j=0;j<num;j++)
			array[j*interval] = y[j];
	}

	void ShiftForSineTransform(tw::Float *array,tw::Int num,tw::Int interval) {
		tw::Int j;
		for (j=num-1;j>=1;j--)
		{
			array[j*interval] = 0.5*(array[j*interval] + array[(j-1)*interval]);
		}
		array[0] = 0.0;
	}

	void ShiftForInverseSineTransform(tw::Float *array,tw::Int num,tw::Int interval) {
		tw::Int j;
		for (j=0;j<num-1;j++)
		{
			array[j*interval] = 0.5*(array[j*interval] + array[(j+1)*interval]);
		}
		array[num-1] *= 0.5;
	}

	// eigenvalues of the laplacian

	inline tw::Float eigenvalue_CFFT(tw::Int i,tw::Int num,tw::Float dxi)
	{
		if (num>1)
			return 2.0*dxi*dxi*(std::cos(tw::Float(2.0*pi*i)/tw::Float(num)) - 1.0);
		else
			return 0.0;
	}

	inline tw::Float eigenvalue_RFFT(tw::Int i,tw::Int num,tw::Float dxi)
	{
		auto trunc = i/2;
		if (num>1)
		{
			if (i!=1)
				return 2.0*dxi*dxi*(std::cos(tw::Float(2.0*pi*trunc)/tw::Float(num)) - 1.0);
			else
				return -4.0*dxi*dxi;
		}
		else
			return 0.0;
	}

	inline tw::Float eigenvalue_FST(tw::Int i,tw::Int num,tw::Float dxi)
	{
		if (num>1)
			return 2.0*dxi*dxi*(std::cos(tw::Float(pi*i)/tw::Float(num)) - 1.0);
		else
			return 0.0;
	}

	inline tw::Float eigenvalue_FCT(tw::Int i,tw::Int num,tw::Float dxi)
	{
		if (num>1)
			return 2.0*dxi*dxi*(std::cos(tw::Float(pi*i)/tw::Float(num)) - 1.0);
		else
			return 0.0;
	}

	// utilities

	void ToggleOrdering(tw::Complex *array,tw::Int num,tw::Int interval) {
		for (auto i=0; i<num/2; i++)
			std::swap(array[i*interval],array[(i+num/2)*interval]);
	}
}
