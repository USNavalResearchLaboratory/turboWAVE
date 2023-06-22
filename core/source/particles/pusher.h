struct BundlePusherBoris : virtual ParticleBundle
{
	alignas(AB) float F[6][N]; // q*F*dth/m
	alignas(AB) tw::Float t[3][N];
	alignas(AB) tw::Float s[3][N];

	BundlePusherBoris(Mover *owner) : ParticleBundle(owner) { ; }
	void impulse(tw::Float u[4][N],float F[6][N]);
	void rotation1(tw::Float t[3][N],tw::Float u[4][N],float F[6][N]);
	void rotation2(tw::Float s[3][N],tw::Float t[3][N]);
	void rotation3(tw::Float s[3][N],tw::Float t[3][N],tw::Float vel[3][N],tw::Float u[4][N]);
	void velocity(tw::Float vel[3][N],tw::Float u[4][N]);
	void Push();
};

struct BundlePusherPGC : BundlePusherBoris
{
	alignas(AB) float las[8][N]; // q*q*a*a/m*m
	alignas(AB) float chi[N];
	alignas(AB) tw::Float avgGam[N];

	BundlePusherPGC(Mover *owner) : ParticleBundle(owner), BundlePusherBoris(owner)  { ; }
	void avg_gam_1(tw::Float avgGam[N],tw::Float vel[3][N],tw::Float u[4][N],float F[6][N],float las[8][N]);
	void avg_gam_2(tw::Float avgGam[N],tw::Float u[4][N],float las[8][N]);
	void impulse(tw::Float u[4][N],float F[6][N],float las[8][N],tw::Float avgGam[N]);
	void rotation1(tw::Float t[3][N],float F[6][N],tw::Float avgGam[N]);
	void velocity(tw::Float vel[3][N],tw::Float u[4][N],tw::Float avgGam[N]);
	void load_chi(float chi[N],float number[N],tw::Float avgGam[N]);
	void Push();
};

struct BundlePusherUnitary : virtual ParticleBundle
{
	alignas(AB) float F[6][N]; // q*F*dth/m
	alignas(AB) float ds[N];
	alignas(AB) tw::Float a[3][N];
	alignas(AB) tw::Float z[2][2][2][N];
	alignas(AB) tw::Float zi[2][2][2][N];
	alignas(AB) tw::Float zf[2][2][2][N];
	BundlePusherUnitary(Mover *owner) : ParticleBundle(owner) { ; }
	void estimate_ds(float ds[N],float F[6][N],tw::Float u[4][N]);
	void copy_spinor(tw::Float dst[2][2][2][N],tw::Float src[2][2][2][N]);
	void add_spinor(tw::Float dst[2][2][2][N],tw::Float src[2][2][2][N]);
	void to_spinor(tw::Float u[4][N],tw::Float z[2][2][2][N]);
	void to_vector(tw::Float u[4][N],tw::Float z[2][2][2][N],tw::Float a[3][N]);
	void dagger(tw::Float z[2][2][2][N]);
	void scalar_mul(tw::Float a[3][N],tw::Float z[2][2][2][N]);
	void left_mul_sig1(tw::Float z[2][2][2][N]);
	void left_mul_sig2(tw::Float z[2][2][2][N]);
	void left_mul_sig3(tw::Float z[2][2][2][N]);
	void set_psi_1(tw::Float a[3][N],float F[6][N]);
	void set_psi_2(tw::Float a[3][N],float F[6][N]);
	void set_psi_3(tw::Float a[3][N],float F[6][N]);
	void set_psi24(tw::Float a[3][N],float F[6][N]);
	void velocity(tw::Float vel[3][N],tw::Float u[4][N]);
	void Lambda();
	void Push();
};

struct BundlePusherBohmian : virtual ParticleBundle
{
	BundlePusherBohmian(Mover *owner) : ParticleBundle(owner) { ; }
	void bohm_velocity(tw::Float vel[3][N],tw::Float u[4][N],float J[4][N]);
	void Push();
};


///////////////////////////////////////////
// UNITARY PUSHER
//////////////////////////////////////////


inline void BundlePusherUnitary::estimate_ds(float ds[N],float F[6][N],tw::Float u[4][N])
{
	// ds is the ratio of proper time step to lab step
	// the dilation factor is folded into F (which already contains dt)
	#pragma omp simd aligned(ds,F,u:AB)
	for (int i=0;i<N;i++)
	{
		ds[i] = 1.0/u[0][i];
		ds[i] *= 1.0 - (F[0][i]*u[1][i] + F[1][i]*u[2][i] + F[2][i]*u[3][i])*dt/sqr(u[0][i]);
	}
	#pragma omp simd aligned(ds,F:AB)
	for (int i=0;i<N;i++)
	{
		F[0][i] *= ds[i];
		F[1][i] *= ds[i];
		F[2][i] *= ds[i];
		F[3][i] *= ds[i];
		F[4][i] *= ds[i];
		F[5][i] *= ds[i];
	}
}
inline void BundlePusherUnitary::copy_spinor(tw::Float dst[2][2][2][N],tw::Float src[2][2][2][N])
{
	#pragma omp simd aligned(dst,src:AB)
	for (int i=0;i<N;i++)
	{
		dst[0][0][0][i] = src[0][0][0][i];
		dst[0][0][1][i] = src[0][0][1][i];
		dst[0][1][0][i] = src[0][1][0][i];
		dst[0][1][1][i] = src[0][1][1][i];
		dst[1][0][0][i] = src[1][0][0][i];
		dst[1][0][1][i] = src[1][0][1][i];
		dst[1][1][0][i] = src[1][1][0][i];
		dst[1][1][1][i] = src[1][1][1][i];
	}
}
inline void BundlePusherUnitary::add_spinor(tw::Float dst[2][2][2][N],tw::Float src[2][2][2][N])
{
	#pragma omp simd aligned(dst,src:AB)
	for (int i=0;i<N;i++)
	{
		dst[0][0][0][i] += src[0][0][0][i];
		dst[0][0][1][i] += src[0][0][1][i];
		dst[0][1][0][i] += src[0][1][0][i];
		dst[0][1][1][i] += src[0][1][1][i];
		dst[1][0][0][i] += src[1][0][0][i];
		dst[1][0][1][i] += src[1][0][1][i];
		dst[1][1][0][i] += src[1][1][0][i];
		dst[1][1][1][i] += src[1][1][1][i];
	}
}
inline void BundlePusherUnitary::to_spinor(tw::Float u[4][N],tw::Float z[2][2][2][N])
{
	// Form the spinor, z, from the four velocity, u
	#pragma omp simd aligned(u,z:AB)
	for (int i=0;i<N;i++)
	{
		// need the following line because u[0] is not part of the global particle state
		u[0][i] = sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i]);
		// real part
		z[0][0][0][i] = u[0][i] + u[3][i];
		z[0][1][0][i] = u[1][i];
		z[1][0][0][i] = u[1][i];
		z[1][1][0][i] = u[0][i] - u[3][i];
		// imag part
		z[0][0][1][i] = 0.0f;
		z[0][1][1][i] = -u[2][i];
		z[1][0][1][i] = u[2][i];
		z[1][1][1][i] = 0.0f;
	}
}
inline void BundlePusherUnitary::to_vector(tw::Float u[4][N],tw::Float z[2][2][2][N],tw::Float a[3][N])
{
	// Restore the four velocity, u, from the spinor, z
	set_psi24(a,F);
	#pragma omp simd aligned(u,z:AB)
	for (int i=0;i<N;i++)
	{
		a[0][i] = 1.0/((1.0 - a[0][i])*(1.0 - a[0][i]) + a[1][i]*a[1][i]);
		u[0][i] = a[0][i]*0.5f*(z[0][0][0][i] + z[1][1][0][i]);
		u[1][i] = a[0][i]*z[0][1][0][i];
		u[2][i] = a[0][i]*z[1][0][1][i];
		u[3][i] = a[0][i]*0.5f*(z[0][0][0][i] - z[1][1][0][i]);
	}
}
inline void BundlePusherUnitary::dagger(tw::Float z[2][2][2][N])
{
	// Form the Hermitian conjugate of z in place
	#pragma omp simd aligned(z:AB)
	for (int i=0;i<N;i++)
	{
		std::swap(z[0][1][0][i],z[1][0][0][i]);
		std::swap(z[0][1][1][i],z[1][0][1][i]);
		z[0][0][1][i] *= -1.0f;
		z[0][1][1][i] *= -1.0f;
		z[1][0][1][i] *= -1.0f;
		z[1][1][1][i] *= -1.0f;
	}
}
inline void BundlePusherUnitary::scalar_mul(tw::Float a[3][N],tw::Float z[2][2][2][N])
{
	// multiply a spinor, z, by a complex scalar, a.
	// the extra element of the scalar is used as a temporary.
	#pragma omp simd aligned(a,z:AB)
	for (int i=0;i<N;i++)
	{
		a[2][i] = z[0][0][0][i];
		z[0][0][0][i] = z[0][0][0][i]*a[0][i] - z[0][0][1][i]*a[1][i];
		z[0][0][1][i] = a[2][i]*a[1][i] + z[0][0][1][i]*a[0][i];

		a[2][i] = z[0][1][0][i];
		z[0][1][0][i] = z[0][1][0][i]*a[0][i] - z[0][1][1][i]*a[1][i];
		z[0][1][1][i] = a[2][i]*a[1][i] + z[0][1][1][i]*a[0][i];

		a[2][i] = z[1][0][0][i];
		z[1][0][0][i] = z[1][0][0][i]*a[0][i] - z[1][0][1][i]*a[1][i];
		z[1][0][1][i] = a[2][i]*a[1][i] + z[1][0][1][i]*a[0][i];

		a[2][i] = z[1][1][0][i];
		z[1][1][0][i] = z[1][1][0][i]*a[0][i] - z[1][1][1][i]*a[1][i];
		z[1][1][1][i] = a[2][i]*a[1][i] + z[1][1][1][i]*a[0][i];
	}
}
inline void BundlePusherUnitary::left_mul_sig1(tw::Float z[2][2][2][N])
{
	// Left multiply by Pauli matrix sigma_1 (swap rows)
	#pragma omp simd aligned(z:AB)
	for (int i=0;i<N;i++)
	{
		std::swap(z[0][0][0][i],z[1][0][0][i]);
		std::swap(z[0][0][1][i],z[1][0][1][i]);
		std::swap(z[0][1][0][i],z[1][1][0][i]);
		std::swap(z[0][1][1][i],z[1][1][1][i]);
	}
}
inline void BundlePusherUnitary::left_mul_sig2(tw::Float z[2][2][2][N])
{
	// Left multiply by Pauli matrix sigma_2
	#pragma omp simd aligned(z:AB)
	for (int i=0;i<N;i++)
	{
		// First multiply top row by i and bottom row by -i
		std::swap(z[0][0][0][i],z[0][0][1][i]);
		z[0][0][0][i] *= -1.0f;
		std::swap(z[0][1][0][i],z[0][1][1][i]);
		z[0][1][0][i] *= -1.0f;
		std::swap(z[1][0][0][i],z[1][0][1][i]);
		z[1][0][1][i] *= -1.0f;
		std::swap(z[1][1][0][i],z[1][1][1][i]);
		z[1][1][1][i] *= -1.0f;
		// Now swap rows
		std::swap(z[0][0][0][i],z[1][0][0][i]);
		std::swap(z[0][0][1][i],z[1][0][1][i]);
		std::swap(z[0][1][0][i],z[1][1][0][i]);
		std::swap(z[0][1][1][i],z[1][1][1][i]);
	}
}
inline void BundlePusherUnitary::left_mul_sig3(tw::Float z[2][2][2][N])
{
	// Left multiply by Pauli matrix sigma_3 (change sign of bottom row)
	#pragma omp simd aligned(z:AB)
	for (int i=0;i<N;i++)
	{
		z[1][0][0][i] *= -1.0f;
		z[1][0][1][i] *= -1.0f;
		z[1][1][0][i] *= -1.0f;
		z[1][1][1][i] *= -1.0f;
	}
}
inline void BundlePusherUnitary::set_psi_1(tw::Float a[3][N],float F[6][N])
{
	#pragma omp simd aligned(a,F:AB)
	for (int i=0;i<N;i++)
	{
		a[0][i] = F[0][i];
		a[1][i] = F[3][i];
	}
}
inline void BundlePusherUnitary::set_psi_2(tw::Float a[3][N],float F[6][N])
{
	#pragma omp simd aligned(a,F:AB)
	for (int i=0;i<N;i++)
	{
		a[0][i] = F[1][i];
		a[1][i] = F[4][i];
	}
}
inline void BundlePusherUnitary::set_psi_3(tw::Float a[3][N],float F[6][N])
{
	#pragma omp simd aligned(a,F:AB)
	for (int i=0;i<N;i++)
	{
		a[0][i] = F[2][i];
		a[1][i] = F[5][i];
	}
}
inline void BundlePusherUnitary::set_psi24(tw::Float a[3][N],float F[6][N])
{
	#pragma omp simd aligned(a,F:AB)
	for (int i=0;i<N;i++)
	{
		a[0][i] = 0.25*(F[0][i]*F[0][i] - F[3][i]*F[3][i] + F[1][i]*F[1][i] - F[4][i]*F[4][i] + F[2][i]*F[2][i] - F[5][i]*F[5][i]);
		a[1][i] = 0.5*(F[0][i]*F[3][i] + F[1][i]*F[4][i] + F[2][i]*F[5][i]);
	}
}
inline void BundlePusherUnitary::velocity(tw::Float vel[3][N],tw::Float u[4][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,u:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = u[c+1][i]/u[0][i];
}

///////////////////////////////////////////
// BORIS PUSHER
//////////////////////////////////////////


inline void BundlePusherBoris::impulse(tw::Float u[4][N],float F[6][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(u,F:AB)
		for (int i=0;i<N;i++)
			u[c+1][i] += F[c][i];
}
inline void BundlePusherBoris::rotation1(tw::Float t[3][N],tw::Float u[4][N],float F[6][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(t,u,F:AB)
		for (int i=0;i<N;i++)
			t[c][i] = F[c+3][i]/sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i]);
}
inline void BundlePusherBoris::rotation2(tw::Float s[3][N],tw::Float t[3][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(s,t:AB)
		for (int i=0;i<N;i++)
			s[c][i] = 2.0*t[c][i]/(1.0 + t[0][i]*t[0][i] + t[1][i]*t[1][i] + t[2][i]*t[2][i]);
}
inline void BundlePusherBoris::rotation3(tw::Float s[3][N],tw::Float t[3][N],tw::Float vel[3][N],tw::Float u[4][N])
{
	// Use vel as temporary while adding rotational impulse, cross(u + cross(u,t),s)
	#pragma omp simd aligned(s,t,vel,u:AB)
	for (int i=0;i<N;i++)
	{
		vel[0][i] = u[1][i] + u[2][i]*t[2][i] - u[3][i]*t[1][i];
		vel[1][i] = u[2][i] + u[3][i]*t[0][i] - u[1][i]*t[2][i];
		vel[2][i] = u[3][i] + u[1][i]*t[1][i] - u[2][i]*t[0][i];
		u[1][i] += vel[1][i]*s[2][i] - vel[2][i]*s[1][i];
		u[2][i] += vel[2][i]*s[0][i] - vel[0][i]*s[2][i];
		u[3][i] += vel[0][i]*s[1][i] - vel[1][i]*s[0][i];
	}
}
inline void BundlePusherBoris::velocity(tw::Float vel[3][N],tw::Float u[4][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,u:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = u[c+1][i]/sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i]);
}


///////////////////////////////////////////
// PGC PUSHER
//////////////////////////////////////////


inline void BundlePusherPGC::avg_gam_1(tw::Float avgGam[N],tw::Float vel[3][N],tw::Float u[4][N],float F[6][N],float las[8][N])
{
	const tw::Float dth = 0.5*dt;
	alignas(AB) tw::Float g1[N];
	alignas(AB) tw::Float dA[N];
	// estimate avgGam at level n using solution of <g>^2 = g1<g> - dA/4
	#pragma omp simd aligned(g1,u,las:AB)
	for (int i=0;i<N;i++)
		g1[i] = sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i] + 0.5*las[7][i]);
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,u,g1:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = u[c+1][i]/g1[i];
	#pragma omp simd aligned(g1,F,vel,dA,las,avgGam:AB)
	for (int i=0;i<N;i++)
	{
		g1[i] += F[0][i]*vel[0][i] + F[1][i]*vel[1][i] + F[2][i]*vel[2][i];
		dA[i] = dth*(las[0][i]*vel[0][i] + las[1][i]*vel[1][i] + las[2][i]*vel[2][i]);
		avgGam[i] = 0.5*g1[i] + 0.5*sqrt(g1[i]*g1[i] - dA[i]);
	}
}
inline void BundlePusherPGC::avg_gam_2(tw::Float avgGam[N],tw::Float u[4][N],float las[8][N])
{
	// estimate avgGam at level n+1/2 using <g>^2 = g1^2 + (g1/<g>)dA/2
	// (solve cubic and expand in powers of dA/g1^2)
	const tw::Float dth = 0.5*dt;
	alignas(AB) tw::Float g1[N];
	alignas(AB) tw::Float dA[N];
	#pragma omp simd aligned(g1,dA,u,las,avgGam:AB)
	for (int i=0;i<N;i++)
	{
		g1[i] = sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i] + 0.5*las[7][i]);
		dA[i] = dth*(las[3][i]*u[1][i] + las[4][i]*u[2][i] + las[5][i]*u[3][i])/g1[i];
		avgGam[i] = g1[i] + 0.25*dA[i]/g1[i] - 0.09375*dA[i]*dA[i]/(g1[i]*g1[i]*g1[i]);
	}
}
inline void BundlePusherPGC::impulse(tw::Float u[4][N],float F[6][N],float las[8][N],tw::Float avgGam[N])
{
	const tw::Float dth = 0.5*dt;
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(u,F,las,avgGam:AB)
		for (int i=0;i<N;i++)
			u[c+1][i] += F[c][i] - 0.25*las[c][i]*dth/avgGam[i];
}
inline void BundlePusherPGC::rotation1(tw::Float t[3][N],float F[6][N],tw::Float avgGam[N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(t,F,avgGam:AB)
		for (int i=0;i<N;i++)
			t[c][i] = F[c+3][i]/avgGam[i];
}
inline void BundlePusherPGC::velocity(tw::Float vel[3][N],tw::Float u[4][N],tw::Float avgGam[N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,u,avgGam:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = u[c+1][i]/avgGam[i];
}
inline void BundlePusherPGC::load_chi(float chi[N],float number[N],tw::Float avgGam[N])
{
	#pragma omp simd aligned(chi,number,avgGam:AB)
	for (int i=0;i<N;i++)
		chi[i] = -q0*q0*number[i]/(m0*avgGam[i]);
}


///////////////////////////////////////////
// BOHMIAN PUSHER
//////////////////////////////////////////


inline void BundlePusherBohmian::bohm_velocity(tw::Float vel[3][N],tw::Float u[4][N],float J[4][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,u,J:AB)
		for (int i=0;i<N;i++)
		{
			const tw::Float vn = J[c+1][i]/(tw::small_pos + J[0][i]);
			vel[c][i] = 1.5*vn - 0.5*u[c+1][i]; // extrapolate to n+1/2
			u[c+1][i] = vn;
		}
}
