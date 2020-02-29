// Special functions and special function dispatch.
// To deal with current compiler differences we wrap some special functions in the tw namespace.
// Definitions passed in from make:
// SF_STD : if defined C++17 special functions are in the std namespace
// SF_ROOT : if defined C++17 special functions are in the root namespace
// If neither are defined use our own special functions.

#ifndef SF_STD
#ifndef SF_ROOT

namespace tw
{

	inline tw::Float cyl_bessel_j(tw::Int n,tw::Float x)
	{
		if (n==0) return j0(x);
		if (n==1) return j1(x);
		return 0;
	}

	inline tw::Float hermite(tw::Int n,tw::Float x)
	{
		tw::Float y0 = 1;
		tw::Float y1 = 2*x;
		tw::Float y2 = n==0 ? 1 : 2*x;
		tw::Int i = 2;
		while (i<=n)
		{
			y2 = 2*x*y1 - 2*tw::Float(i-1)*y0;
			y0 = y1;
			y1 = y2;
			i++;
		}
		return y2;
	}

	inline tw::Float assoc_laguerre(tw::Int n,tw::Int m,tw::Float x)
	{
		if (n==0)
			return 1.0;
		if (n==1)
			return 1.0 + m - x;
		if (n==2)
			return 0.5*(x*x - 2.0*(m+2.0)*x + (m+1.0)*(m+2.0));
		if (n==3)
			return (6.0 + 11.0*m + 6.0*m*m + m*m*m - 18.0*x - 15.0*m*x - 3.0*m*m*x + 9.0*x*x + 3*m*x*x - x*x*x)/6.0;
		return 0.0;
	}

}

#endif
#endif

#ifdef SF_ROOT

namespace tw
{

	inline tw::Float cyl_bessel_j(tw::Int n,tw::Float x) { return cyl_bessel_j(n,x); }

	inline tw::Float hermite(tw::Int n,tw::Float x) { return hermite(n,x); }

	inline tw::Float assoc_laguerre(tw::Int n,tw::Int m,tw::Float x) { return assoc_laguerre(n,m,x); }

}

#endif

#ifdef SF_STD

namespace tw
{

	inline tw::Float cyl_bessel_j(tw::Int n,tw::Float x) { return std::cyl_bessel_j(n,x); }

	inline tw::Float hermite(tw::Int n,tw::Float x) { return std::hermite(n,x); }

	inline tw::Float assoc_laguerre(tw::Int n,tw::Int m,tw::Float x) { return std::assoc_laguerre(n,m,x); }

}

#endif

namespace tw
{

	inline tw::Complex erfcc(const tw::Complex x)
	{
		// modified from NR to support complex argument
		tw::Complex z,t,ans;
		const tw::Complex z0(1.26551223,0.0);
		const tw::Complex z1(1.00002368,0.0);
		const tw::Complex z2(0.37409196,0.0);
		const tw::Complex z3(0.09678418,0.0);
		const tw::Complex z4(0.18628806,0.0);
		const tw::Complex z5(0.27886807,0.0);
		const tw::Complex z6(1.13520398,0.0);
		const tw::Complex z7(1.48851587,0.0);
		const tw::Complex z8(0.82215223,0.0);
		const tw::Complex z9(0.17087277,0.0);
		z = tw::Complex(fabs(real(x)),imag(x));
		t=two/(two+z);
		ans=t*std::exp(-z*z-z0+t*(z1+t*(z2+t*(z3+t*(-z4+t*(z5+t*(-z6+t*(z7+t*(-z8+t*z9)))))))));
		return real(x) >= 0.0 ? ans : tw::Complex(2-real(ans),imag(ans));
	}

	inline tw::Float erfi(const tw::Float x)
	{
		return real((one - erfcc(ii*x))/ii);
	}

}

inline tw::Float Factorial(tw::Int n)
{
	tw::Float ans[20] = { 1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880.,3628800.,39916800.,479001600.,
		6227020800.,87178291200.,1307674368000.,20922789888000.,355687428096000.,
		6402373705728000.,121645100408832000. };
	return ans[n];
}

inline tw::Complex SphericalHarmonic(tw::Int l,tw::Int m,tw::Float phi,tw::Float theta)
{
	tw::Int ma = std::abs(m);
	tw::Complex ans;
	if (l==0)
		ans = 1.0/sqrt(4.0*pi);
	if (l==1 && ma==0)
		ans = ii*tw::Float(cos(theta)*sqrt(3.0/(4.0*pi)));
	if (l==1 && ma==1)
		ans = -ii*tw::Float(sin(theta)*sqrt(3.0/(8.0*pi)));
	if (l==2 && ma==0)
		ans = sqrt(5.0/(16.0*pi))*(1.0-3.0*sqr(cos(theta)));
	if (l==2 && ma==1)
		ans = sqrt(15.0/(8.0*pi))*cos(theta)*sin(theta);
	if (l==2 && ma==2)
		ans = -sqrt(15.0/(32.0*pi))*sqr(sin(theta));
	if (m<0)
		ans = conj(ans);
	ans *= std::exp(ii*tw::Float(m)*phi);
	return tw::Float(sqrt(4.0*pi))*ans;
}

inline tw::spinor SphericalHarmonicSpinor(tw::Float j,tw::Int l,tw::Float m,tw::Float phi,tw::Float theta)
{
	// L&L Omega_jlm spinor components
	// j must be l+1/2 or l-1/2
	// m must be half integral
	const tw::Complex Q_up = SphericalHarmonic(l,std::floor(m-0.49),phi,theta);
	const tw::Complex Q_down = SphericalHarmonic(l,std::floor(m+0.51),phi,theta);
	if (j>tw::Float(l))
		return tw::spinor(tw::Float(sqrt((j+m)/(2*j)))*Q_up,tw::Float(sqrt((j-m)/(2*j)))*Q_down);
	else
		return tw::spinor(-tw::Float(sqrt((j-m+1)/(2*j+2)))*Q_up,tw::Float(sqrt((j+m+1)/(2*j+2)))*Q_down);
}

inline tw::cvec3 SphericalHarmonicVectorMagnetic(tw::Int j,tw::Int m,tw::Float phi,tw::Float theta)
{
	// Magnetic multipole function
	// components returned as radial, azimuthal, polar
	if (j==0)
		return tw::cvec3(0.0,0.0,0.0);
	if (j==1 && m==0)
		return -sqrt(0.75/pi)*sin(theta)*tw::cvec3(0.0,1.0,0.0);
	if (j==1 && m==1)
		return sqrt(0.375/pi)*std::exp(ii*phi)*tw::cvec3(0.0,-cos(theta),ii);
	if (j==2 && m==0)
		return -1.5*sqrt(5.0/pi)*sin(theta)*cos(theta)*tw::cvec3(0.0,1.0,0.0);
	if (j==2 && m==1)
		return sqrt(1.875/pi)*std::exp(ii*phi)*tw::cvec3(0.0,-cos(2*theta),ii*cos(theta));
	return tw::cvec3(0.0);
}

inline tw::Float ConfluentHypergeometric(tw::Int n,tw::Float gam,tw::Float z)
{
	// assume n is 0 or negative integer
	// for an atom n is typically the negative of the radial quantum number
	// Mathematica calls this function Hypergeometric1F1
	if (n==0)
		return 1.0;
	if (n==-1)
		return 1.0 - z/gam;
	if (n==-2)
		return 1.0 - 2.0*z/gam + z*z/(gam+gam*gam);
	if (n==-3)
		return 1.0 + 3.0*z*z/(gam+gam*gam) - (3.0*z + z*z*z/(2.0+3.0*gam+gam*gam))/gam;
	return 0.0;
}

// rise between x=0 and x=1
inline tw::Float QuinticRise(tw::Float x)
{
	return 10.0*pow(x,tw::Float(3.0)) - 15.0*pow(x,tw::Float(4.0)) + 6.0*pow(x,tw::Float(5.0));
}

// base to base from x=0 to x=1
inline tw::Float QuinticPulse(tw::Float x)
{
	if (x>0.5)
		return QuinticRise((1.0-x)*2.0);
	else
		return QuinticRise(x*2.0);
}

inline tw::Float D1QuinticRise(tw::Float x)
{
	return 30.0*pow(x,tw::Float(2.0)) - 60.0*pow(x,tw::Float(3.0)) + 30.0*pow(x,tw::Float(4.0));
}

inline tw::Float D1QuinticPulse(tw::Float x)
{
	if (x>0.5)
		return 2.0*D1QuinticRise((1.0-x)*2.0);
	else
		return 2.0*D1QuinticRise(x*2.0);
}

inline tw::Float QuarticRise(tw::Float x)
{
	return x*x*(x-2)*(x-2); // same as 1/2 of D1QuinticRise, normalized
}

inline tw::Float QuarticPulse(tw::Float x)
{
	if (x>0.5)
		return QuarticRise((1.0-x)*2.0);
	else
		return QuarticRise(x*2.0);
}

inline tw::Float D1QuarticRise(tw::Float x)
{
	return 4.0*x*(2.0-3.0*x+x*x);
}

inline tw::Float D1QuarticPulse(tw::Float x)
{
	if (x>0.5)
		return 2.0*D1QuarticRise((1.0-x)*2.0);
	else
		return 2.0*D1QuarticRise(x*2.0);
}

inline tw::Float Sin2Rise(tw::Float x)
{
	return sqr(sin(0.5*pi*x));
}

inline tw::Float Sin2Pulse(tw::Float x)
{
	if (x>0.5)
		return Sin2Rise((1.0-x)*2.0);
	else
		return Sin2Rise(x*2.0);
}

inline tw::Float D1Sin2Rise(tw::Float x)
{
	return 0.5*pi*sin(pi*x);
}

inline tw::Float D1Sin2Pulse(tw::Float x)
{
	if (x>0.5)
		return 2.0*D1Sin2Rise((1.0-x)*2.0);
	else
		return 2.0*D1Sin2Rise(x*2.0);
}

inline tw::Float SechRise(tw::Float x)
{
		if (x>0.33)
			return 1.0/cosh((1.0-x)*6.0);
		return sqr(cos(((1.0-x)*6.0-4.0)*pi/4.0))/cosh((1.0-x)*6.0);
}

inline tw::Float SechPulse(tw::Float x)
{
	if (x>0.5)
		return SechRise((1.0-x)*2.0);
	else
		return SechRise(x*2.0);
}

inline tw::Float D1SechRise(tw::Float x)
{
		if (x>0.33)
			return 6.0*tanh(6.0-6.0*x)/cosh(6.0-6.0*x);
		return (3.0/cosh(6.0-6.0*x))*sin(1.5*pi*x)*(pi*cos(1.5*pi*x)+2.0*sin(1.5*pi*x)*tanh(6.0-6.0*x));
}

inline tw::Float D1SechPulse(tw::Float x)
{
	if (x>0.5)
		return 2.0*D1SechRise((1.0-x)*2.0);
	else
		return 2.0*D1SechRise(x*2.0);
}
