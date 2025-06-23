module;

#include "tw_includes.h"

export module pusher;
import bundle;

/// @brief Select subcycling for this step based on the field
/// @param F0 the angle tensor **will be rescaled** for the substep duration
/// @return number of subcycles
export tw::Int set_subcycles(float F0[6][tw::max_bundle_size])
{
	const tw::Int N = tw::max_bundle_size;
	const tw::Int AB = tw::vec_align_bytes;
	float max2 = 0.0;
	for (int c = 0; c < 6; c++)
		for (int i = 0; i < N; i++)
		{
			if (sqr(2 * F0[c][i]) > max2)
				max2 = sqr(2 * F0[c][i]);
		}
	tw::Int subcycles = 1 + tw::Int(std::sqrt(max2));
	const float subi = 1.0f / subcycles;
	for (int c = 0; c < 6; c++)
#pragma omp simd aligned(F0:AB)
		for (int i = 0; i < N; i++)
			F0[c][i] *= subi;
	return subcycles;
}

export struct BundlePusherBoris : virtual ParticleBundle
{
	alignas(AB) float F[6][N]; // q*F*dth/m
	alignas(AB) tw::Float t[3][N];
	alignas(AB) tw::Float s[3][N];

	BundlePusherBoris(const MoverParams& mov) : ParticleBundle(mov) { ; }
	void impulse(tw::Float u[4][N], float F[6][N]);
	void rotation1(tw::Float t[3][N], tw::Float u[4][N], float F[6][N]);
	void rotation2(tw::Float s[3][N], tw::Float t[3][N]);
	void rotation3(tw::Float s[3][N], tw::Float t[3][N], tw::Float vel[4][N], tw::Float u[4][N]);
	void velocity(tw::Float vel[4][N], tw::Float u[4][N]);
	void Push(tw::Float dts);
};

export struct BundlePusherHC : BundlePusherBoris
{
	alignas(AB) float gammaNew[N];
	alignas(AB) float b[N];

	BundlePusherHC(const MoverParams& mov) : ParticleBundle(mov), BundlePusherBoris(mov) { ; }
	void impulse(tw::Float u[4][N], float F[6][N]);
	void gamma_new(float gammaNew[N], float b[N], tw::Float u[4][N], float F[6][N]);
	void rotation1(tw::Float t[3][N], float gammaNew[N], float F[6][N]);
	void rotation2(tw::Float s[3][N], tw::Float t[3][N]);
	void rotation3(tw::Float s[3][N], tw::Float t[3][N], tw::Float vel[4][N], tw::Float u[4][N]);
	void velocity(tw::Float vel[4][N], tw::Float u[4][N]);
	void Push(tw::Float dts);
};

export struct BundlePusherPGC : BundlePusherBoris
{
	alignas(AB) float las[8][N]; // q*q*a*a*dth/m*m
	alignas(AB) float chi[N];
	alignas(AB) tw::Float avgGam[N];

	BundlePusherPGC(const MoverParams& mov) : ParticleBundle(mov), BundlePusherBoris(mov) { ; }
	void avg_gam_1(tw::Float avgGam[N], tw::Float vel[4][N], tw::Float u[4][N], float F[6][N], float las[8][N]);
	void avg_gam_2(tw::Float avgGam[N], tw::Float u[4][N], float las[8][N]);
	void impulse(tw::Float u[4][N], float F[6][N], float las[8][N], tw::Float avgGam[N]);
	void rotation1(tw::Float t[3][N], float F[6][N], tw::Float avgGam[N]);
	void velocity(tw::Float vel[4][N], tw::Float u[4][N], tw::Float avgGam[N]);
	void load_chi(float chi[N], float number[N], tw::Float avgGam[N]);
	void Push(tw::Float dts);
};

export struct BundlePusherUnitary : virtual ParticleBundle
{
	alignas(AB) float F0[6][N]; // q*F*dth/m
	alignas(AB) float F[6][N]; // includes individual dilation factors
	alignas(AB) tw::Float a[3][N];
	alignas(AB) tw::Float z[2][2][2][N];
	alignas(AB) tw::Float zi[2][2][2][N];
	alignas(AB) tw::Float zf[2][2][2][N];
	BundlePusherUnitary(const MoverParams& mov) : ParticleBundle(mov) { ; }
	void dilate_timestep(float F0[6][N], float F[6][N], tw::Float u[4][N]);
	void copy_spinor(tw::Float dst[2][2][2][N], tw::Float src[2][2][2][N]);
	void add_spinor(tw::Float dst[2][2][2][N], tw::Float src[2][2][2][N]);
	void to_spinor(tw::Float u[4][N], tw::Float z[2][2][2][N]);
	void to_vector(tw::Float u[4][N], tw::Float z[2][2][2][N], tw::Float a[3][N]);
	void dagger(tw::Float z[2][2][2][N]);
	void scalar_mul(tw::Float a[3][N], tw::Float z[2][2][2][N]);
	void left_mul_sig1(tw::Float z[2][2][2][N]);
	void left_mul_sig2(tw::Float z[2][2][2][N]);
	void left_mul_sig3(tw::Float z[2][2][2][N]);
	void set_psi_1(tw::Float a[3][N], float F[6][N]);
	void set_psi_2(tw::Float a[3][N], float F[6][N]);
	void set_psi_3(tw::Float a[3][N], float F[6][N]);
	void set_psi24(tw::Float a[3][N], float F[6][N]);
	void velocity(tw::Float vel[4][N], tw::Float u[4][N]);
	void Lambda();
	void Push(tw::Float dts);
};

export struct BundlePusherBohmian : virtual ParticleBundle
{
	BundlePusherBohmian(const MoverParams& mov) : ParticleBundle(mov) { ; }
	void bohm_velocity(tw::Float vel[4][N], tw::Float u[4][N], float J[4][N]);
	void Push(tw::Float dts);
};

export struct BundlePusherPhoton : virtual ParticleBundle
{
	BundlePusherPhoton(const MoverParams& mov) : ParticleBundle(mov) { ; }
	void velocity(tw::Float vel[4][N], tw::Float u[4][N]);
	void Push(tw::Float dts);
};


///////////////////////////////////////////
// UNITARY PUSHER
//////////////////////////////////////////


inline void BundlePusherUnitary::copy_spinor(tw::Float dst[2][2][2][N], tw::Float src[2][2][2][N])
{
#pragma omp simd aligned(dst,src:AB)
	for (int i = 0; i < N; i++)
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
inline void BundlePusherUnitary::add_spinor(tw::Float dst[2][2][2][N], tw::Float src[2][2][2][N])
{
#pragma omp simd aligned(dst,src:AB)
	for (int i = 0; i < N; i++)
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
inline void BundlePusherUnitary::to_spinor(tw::Float u[4][N], tw::Float z[2][2][2][N])
{
	// Form the spinor, z, from the four velocity, u
#pragma omp simd aligned(u,z:AB)
	for (int i = 0; i < N; i++)
	{
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
inline void BundlePusherUnitary::to_vector(tw::Float u[4][N], tw::Float z[2][2][2][N], tw::Float a[3][N])
{
	// Restore the four velocity, u, from the spinor, z
	set_psi24(a, F);
#pragma omp simd aligned(u,z:AB)
	for (int i = 0; i < N; i++)
	{
		a[0][i] = 1.0 / ((1.0 - a[0][i]) * (1.0 - a[0][i]) + a[1][i] * a[1][i]);
		u[0][i] = a[0][i] * 0.5f * (z[0][0][0][i] + z[1][1][0][i]);
		u[1][i] = a[0][i] * z[0][1][0][i];
		u[2][i] = a[0][i] * z[1][0][1][i];
		u[3][i] = a[0][i] * 0.5f * (z[0][0][0][i] - z[1][1][0][i]);
	}
}
inline void BundlePusherUnitary::dagger(tw::Float z[2][2][2][N])
{
	// Form the Hermitian conjugate of z in place
#pragma omp simd aligned(z:AB)
	for (int i = 0; i < N; i++)
	{
		std::swap(z[0][1][0][i], z[1][0][0][i]);
		std::swap(z[0][1][1][i], z[1][0][1][i]);
		z[0][0][1][i] *= -1.0f;
		z[0][1][1][i] *= -1.0f;
		z[1][0][1][i] *= -1.0f;
		z[1][1][1][i] *= -1.0f;
	}
}
inline void BundlePusherUnitary::scalar_mul(tw::Float a[3][N], tw::Float z[2][2][2][N])
{
	// multiply a spinor, z, by a complex scalar, a.
	// the extra element of the scalar is used as a temporary.
#pragma omp simd aligned(a,z:AB)
	for (int i = 0; i < N; i++)
	{
		a[2][i] = z[0][0][0][i];
		z[0][0][0][i] = z[0][0][0][i] * a[0][i] - z[0][0][1][i] * a[1][i];
		z[0][0][1][i] = a[2][i] * a[1][i] + z[0][0][1][i] * a[0][i];

		a[2][i] = z[0][1][0][i];
		z[0][1][0][i] = z[0][1][0][i] * a[0][i] - z[0][1][1][i] * a[1][i];
		z[0][1][1][i] = a[2][i] * a[1][i] + z[0][1][1][i] * a[0][i];

		a[2][i] = z[1][0][0][i];
		z[1][0][0][i] = z[1][0][0][i] * a[0][i] - z[1][0][1][i] * a[1][i];
		z[1][0][1][i] = a[2][i] * a[1][i] + z[1][0][1][i] * a[0][i];

		a[2][i] = z[1][1][0][i];
		z[1][1][0][i] = z[1][1][0][i] * a[0][i] - z[1][1][1][i] * a[1][i];
		z[1][1][1][i] = a[2][i] * a[1][i] + z[1][1][1][i] * a[0][i];
	}
}
inline void BundlePusherUnitary::left_mul_sig1(tw::Float z[2][2][2][N])
{
	// Left multiply by Pauli matrix sigma_1 (swap rows)
#pragma omp simd aligned(z:AB)
	for (int i = 0; i < N; i++)
	{
		std::swap(z[0][0][0][i], z[1][0][0][i]);
		std::swap(z[0][0][1][i], z[1][0][1][i]);
		std::swap(z[0][1][0][i], z[1][1][0][i]);
		std::swap(z[0][1][1][i], z[1][1][1][i]);
	}
}
inline void BundlePusherUnitary::left_mul_sig2(tw::Float z[2][2][2][N])
{
	// Left multiply by Pauli matrix sigma_2
#pragma omp simd aligned(z:AB)
	for (int i = 0; i < N; i++)
	{
		// First multiply top row by i and bottom row by -i
		std::swap(z[0][0][0][i], z[0][0][1][i]);
		z[0][0][0][i] *= -1.0f;
		std::swap(z[0][1][0][i], z[0][1][1][i]);
		z[0][1][0][i] *= -1.0f;
		std::swap(z[1][0][0][i], z[1][0][1][i]);
		z[1][0][1][i] *= -1.0f;
		std::swap(z[1][1][0][i], z[1][1][1][i]);
		z[1][1][1][i] *= -1.0f;
		// Now swap rows
		std::swap(z[0][0][0][i], z[1][0][0][i]);
		std::swap(z[0][0][1][i], z[1][0][1][i]);
		std::swap(z[0][1][0][i], z[1][1][0][i]);
		std::swap(z[0][1][1][i], z[1][1][1][i]);
	}
}
inline void BundlePusherUnitary::left_mul_sig3(tw::Float z[2][2][2][N])
{
	// Left multiply by Pauli matrix sigma_3 (change sign of bottom row)
#pragma omp simd aligned(z:AB)
	for (int i = 0; i < N; i++)
	{
		z[1][0][0][i] *= -1.0f;
		z[1][0][1][i] *= -1.0f;
		z[1][1][0][i] *= -1.0f;
		z[1][1][1][i] *= -1.0f;
	}
}
inline void BundlePusherUnitary::set_psi_1(tw::Float a[3][N], float F[6][N])
{
#pragma omp simd aligned(a,F:AB)
	for (int i = 0; i < N; i++)
	{
		a[0][i] = F[0][i];
		a[1][i] = F[3][i];
	}
}
inline void BundlePusherUnitary::set_psi_2(tw::Float a[3][N], float F[6][N])
{
#pragma omp simd aligned(a,F:AB)
	for (int i = 0; i < N; i++)
	{
		a[0][i] = F[1][i];
		a[1][i] = F[4][i];
	}
}
inline void BundlePusherUnitary::set_psi_3(tw::Float a[3][N], float F[6][N])
{
#pragma omp simd aligned(a,F:AB)
	for (int i = 0; i < N; i++)
	{
		a[0][i] = F[2][i];
		a[1][i] = F[5][i];
	}
}
inline void BundlePusherUnitary::set_psi24(tw::Float a[3][N], float F[6][N])
{
#pragma omp simd aligned(a,F:AB)
	for (int i = 0; i < N; i++)
	{
		a[0][i] = 0.25 * (F[0][i] * F[0][i] - F[3][i] * F[3][i] + F[1][i] * F[1][i] - F[4][i] * F[4][i] + F[2][i] * F[2][i] - F[5][i] * F[5][i]);
		a[1][i] = 0.5 * (F[0][i] * F[3][i] + F[1][i] * F[4][i] + F[2][i] * F[5][i]);
	}
}
inline void BundlePusherUnitary::velocity(tw::Float vel[4][N], tw::Float u[4][N])
{
	for (int c = 0; c < 4; c++)
#pragma omp simd aligned(vel,u:AB)
		for (int i = 0; i < N; i++)
			vel[c][i] = u[c][i] / u[0][i];
}

void BundlePusherUnitary::Lambda()
{
	// Left multiply the spinor by Lambda (time translation operator)
	copy_spinor(zf, zi);

	copy_spinor(z, zi);
	set_psi24(a, F);
	scalar_mul(a, z);
	add_spinor(zf, z);

	copy_spinor(z, zi);
	left_mul_sig1(z);
	set_psi_1(a, F);
	scalar_mul(a, z);
	add_spinor(zf, z);

	copy_spinor(z, zi);
	left_mul_sig2(z);
	set_psi_2(a, F);
	scalar_mul(a, z);
	add_spinor(zf, z);

	copy_spinor(z, zi);
	left_mul_sig3(z);
	set_psi_3(a, F);
	scalar_mul(a, z);
	add_spinor(zf, z);
}

void BundlePusherUnitary::dilate_timestep(float F0[6][N], float F[6][N], tw::Float u[4][N])
{
	// ds[o][i] is the ratio of proper time step to lab step of order o+1 for particle i
	alignas(AB) tw::Float ds[4][N], E2[N], B2[N], udotE[N], udotB[N], EdotB[N], udotExB[N], tst[N];

	// TODO: is this needed/desirable
#pragma omp simd aligned(u:AB)
	for (int i = 0; i < N; i++)
		u[0][i] = std::sqrt(1.0 + u[1][i] * u[1][i] + u[2][i] * u[2][i] + u[3][i] * u[3][i]);

	for (int c = 0; c < 6; c++)
#pragma omp simd aligned(F0,F:AB)
		for (int i = 0; i < N; i++)
			F[c][i] = F0[c][i];

#pragma omp simd aligned(ds,F,u,E2,B2,udotE,udotB,EdotB,udotExB,tst:AB)
	for (int i = 0; i < N; i++)
	{
		udotE[i] = F[0][i] * u[1][i] + F[1][i] * u[2][i] + F[2][i] * u[3][i];
		udotB[i] = F[3][i] * u[1][i] + F[4][i] * u[2][i] + F[5][i] * u[3][i];
		EdotB[i] = F[0][i] * F[3][i] + F[1][i] * F[4][i] + F[2][i] * F[5][i];
		E2[i] = F[0][i] * F[0][i] + F[1][i] * F[1][i] + F[2][i] * F[2][i];
		B2[i] = F[3][i] * F[3][i] + F[4][i] * F[4][i] + F[5][i] * F[5][i];
		udotExB[i] =
			u[1][i] * (F[1][i] * F[5][i] - F[2][i] * F[4][i]) +
			u[2][i] * (F[2][i] * F[3][i] - F[0][i] * F[5][i]) +
			u[3][i] * (F[0][i] * F[4][i] - F[1][i] * F[3][i]);

		ds[0][i] = 1 / u[0][i];

		tst[i] = udotE[i] * sqr(ds[0][i]);
		ds[1][i] = ds[0][i] * (1 - tw::Float(tst[i] < 0) * tst[i]) / (1 + tw::Float(tst[i] >= 0) * tst[i]);
		tst[i] = tw::Float(ds[1][i] < 1.0);
		ds[1][i] = tst[i] * ds[1][i] + (1 - tst[i]);

		tst[i] = udotE[i] * sqr(ds[1][i]) + 0.66667 * (u[0][i] * E2[i] - udotExB[i]) * cub(ds[1][i]);
		ds[2][i] = ds[0][i] * (1 - tw::Float(tst[i] < 0) * tst[i]) / (1 + tw::Float(tst[i] >= 0) * tst[i]);
		tst[i] = tw::Float(ds[2][i] < 1.0);
		ds[2][i] = tst[i] * ds[2][i] + (1 - tst[i]);

		tst[i] = udotE[i] * sqr(ds[2][i]) + 0.66667 * (u[0][i] * E2[i] - udotExB[i]) * cub(ds[2][i])
			+ 0.125 * (3 * udotE[i] * (E2[i] - B2[i]) + 2 * udotB[i] * EdotB[i]) * quad(ds[2][i]);
		ds[3][i] = ds[0][i] * (1 - tw::Float(tst[i] < 0) * tst[i]) / (1 + tw::Float(tst[i] >= 0) * tst[i]);
		tst[i] = tw::Float(ds[3][i] < 1.0);
		ds[3][i] = tst[i] * ds[3][i] + (1 - tst[i]);

		F[0][i] *= ds[3][i];
		F[1][i] *= ds[3][i];
		F[2][i] *= ds[3][i];
		F[3][i] *= ds[3][i];
		F[4][i] *= ds[3][i];
		F[5][i] *= ds[3][i];
	}
}

void BundlePusherUnitary::Push(tw::Float dts)
{
	int subcycles = set_subcycles(F0);
	for (int sub = 0; sub < subcycles; sub++)
	{
		// We rewrite L*z*L^dag as (L*(L*z)^dag)^dag
		dilate_timestep(F0, F, u); // dilation gets buried in F
		to_spinor(u, zi);
		Lambda();
		copy_spinor(zi, zf);
		dagger(zi);
		Lambda();
		dagger(zf);
		to_vector(u, zf, a);
	}
	velocity(vel, u);
	translate(x, vel, dts);
	load_j4(J, number, vel);
}

///////////////////////////////////////////
// BORIS PUSHER
//////////////////////////////////////////


inline void BundlePusherBoris::impulse(tw::Float u[4][N], float F[6][N])
{
	for (int c = 0; c < 3; c++)
#pragma omp simd aligned(u,F:AB)
		for (int i = 0; i < N; i++)
			u[c + 1][i] += F[c][i];
}
inline void BundlePusherBoris::rotation1(tw::Float t[3][N], tw::Float u[4][N], float F[6][N])
{
	for (int c = 0; c < 3; c++)
#pragma omp simd aligned(t,u,F:AB)
		for (int i = 0; i < N; i++)
			t[c][i] = F[c + 3][i] / std::sqrt(1.0 + u[1][i] * u[1][i] + u[2][i] * u[2][i] + u[3][i] * u[3][i]);
}
inline void BundlePusherBoris::rotation2(tw::Float s[3][N], tw::Float t[3][N])
{
	for (int c = 0; c < 3; c++)
#pragma omp simd aligned(s,t:AB)
		for (int i = 0; i < N; i++)
			s[c][i] = 2.0 * t[c][i] / (1.0 + t[0][i] * t[0][i] + t[1][i] * t[1][i] + t[2][i] * t[2][i]);
}
inline void BundlePusherBoris::rotation3(tw::Float s[3][N], tw::Float t[3][N], tw::Float vel[4][N], tw::Float u[4][N])
{
	// Use vel as temporary while adding rotational impulse, cross(u + cross(u,t),s)
#pragma omp simd aligned(s,t,vel,u:AB)
	for (int i = 0; i < N; i++)
	{
		vel[1][i] = u[1][i] + u[2][i] * t[2][i] - u[3][i] * t[1][i];
		vel[2][i] = u[2][i] + u[3][i] * t[0][i] - u[1][i] * t[2][i];
		vel[3][i] = u[3][i] + u[1][i] * t[1][i] - u[2][i] * t[0][i];
		u[1][i] += vel[2][i] * s[2][i] - vel[3][i] * s[1][i];
		u[2][i] += vel[3][i] * s[0][i] - vel[1][i] * s[2][i];
		u[3][i] += vel[1][i] * s[1][i] - vel[2][i] * s[0][i];
	}
}
inline void BundlePusherBoris::velocity(tw::Float vel[4][N], tw::Float u[4][N])
{
#pragma omp simd aligned(u:AB)
	for (int i = 0; i < N; i++)
		u[0][i] = std::sqrt(1.0 + u[1][i] * u[1][i] + u[2][i] * u[2][i] + u[3][i] * u[3][i]);

	for (int c = 0; c < 4; c++)
#pragma omp simd aligned(vel,u:AB)
		for (int i = 0; i < N; i++)
			vel[c][i] = u[c][i] / u[0][i];
}

void BundlePusherBoris::Push(tw::Float dts)
{
	int subcycles = set_subcycles(F);
	for (int sub = 0; sub < subcycles; sub++)
	{
		impulse(u, F);
		rotation1(t, u, F);
		rotation2(s, t);
		rotation3(s, t, vel, u);
		impulse(u, F);
	}
	velocity(vel, u);
	translate(x, vel, dts);
	load_j4(J, number, vel);
}

///////////////////////////////////////////
// HC PUSHER
///////////////////////////////////////////


inline void BundlePusherHC::impulse(tw::Float u[4][N], float F[6][N])
{
	for (int c = 0; c < 3; c++)
#pragma omp simd aligned(u,F:AB)
		for (int i = 0; i < N; i++)
			u[c + 1][i] += F[c][i];
}
inline void BundlePusherHC::gamma_new(float gammaNew[N], float b[N], tw::Float u[4][N], float F[6][N])
{
#pragma omp simd aligned(b,u,F:AB)
	for (int i = 0; i < N; i++)
		b[i] = (1.0 + u[1][i] * u[1][i] + u[2][i] * u[2][i] + u[3][i] * u[3][i]) - (F[3][i] * F[3][i] + F[4][i] * F[4][i] + F[5][i] * F[5][i]);

#pragma omp simd aligned(gammaNew,b,u,F:AB)
	for (int i = 0; i < N; i++)
	{
		gammaNew[i] = F[3][i] * u[1][i] + F[4][i] * u[2][i] + F[5][i] * u[3][i];
		gammaNew[i] = 4.0 * (gammaNew[i] + F[3][i] * F[3][i] + F[4][i] * F[4][i] + F[5][i] * F[5][i]);
		gammaNew[i] = std::sqrt(0.5 * (b[i] + std::sqrt(b[i] * b[i] + gammaNew[i])));
	}
}
inline void BundlePusherHC::rotation1(tw::Float t[3][N], float gammaNew[N], float F[6][N])
{
	for (int c = 0; c < 3; c++)
#pragma omp simd aligned(t,gammaNew,F:AB)
		for (int i = 0; i < N; i++)
			t[c][i] = F[c + 3][i] / gammaNew[i];
}
inline void BundlePusherHC::rotation2(tw::Float s[3][N], tw::Float t[3][N])
{
	for (int c = 0; c < 3; c++)
#pragma omp simd aligned(s,t:AB)
		for (int i = 0; i < N; i++)
			s[c][i] = 2.0 * t[c][i] / (1.0 + t[0][i] * t[0][i] + t[1][i] * t[1][i] + t[2][i] * t[2][i]);
}
inline void BundlePusherHC::rotation3(tw::Float s[3][N], tw::Float t[3][N], tw::Float vel[4][N], tw::Float u[4][N])
{
	// Use vel as temporary while adding rotational impulse, cross(u + cross(u,t),s)
#pragma omp simd aligned(s,t,vel,u:AB)
	for (int i = 0; i < N; i++)
	{
		vel[1][i] = u[1][i] + u[2][i] * t[2][i] - u[3][i] * t[1][i];
		vel[2][i] = u[2][i] + u[3][i] * t[0][i] - u[1][i] * t[2][i];
		vel[3][i] = u[3][i] + u[1][i] * t[1][i] - u[2][i] * t[0][i];
		u[1][i] += vel[2][i] * s[2][i] - vel[3][i] * s[1][i];
		u[2][i] += vel[3][i] * s[0][i] - vel[1][i] * s[2][i];
		u[3][i] += vel[1][i] * s[1][i] - vel[2][i] * s[0][i];
	}
}

inline void BundlePusherHC::velocity(tw::Float vel[4][N], tw::Float u[4][N])
{
#pragma omp simd aligned(u:AB)
	for (int i = 0; i < N; i++)
		u[0][i] = std::sqrt(1.0 + u[1][i] * u[1][i] + u[2][i] * u[2][i] + u[3][i] * u[3][i]);

	for (int c = 0; c < 4; c++)
#pragma omp simd aligned(vel,u:AB)
		for (int i = 0; i < N; i++)
			vel[c][i] = u[c][i] / u[0][i];
}

void BundlePusherHC::Push(tw::Float dts)
{
	int subcycles = set_subcycles(F);
	for (int sub = 0; sub < subcycles; sub++)
	{
		impulse(u, F);
		gamma_new(gammaNew, b, u, F);
		rotation1(t, gammaNew, F);
		rotation2(s, t);
		rotation3(s, t, vel, u);
		impulse(u, F);
	}
	velocity(vel, u);
	translate(x, vel, dts);
	load_j4(J, number, vel);
}


///////////////////////////////////////////
// PGC PUSHER
//////////////////////////////////////////


inline void BundlePusherPGC::avg_gam_1(tw::Float avgGam[N], tw::Float vel[4][N], tw::Float u[4][N], float F[6][N], float las[8][N])
{
	// las = [q2m2dth*grad(a^2(n)), q2m2dth*grad(a^2(n+1/2)), q2m2h*a^2(n), q2m2h*a^2(n+1/2)]
	// we are estimating std::sqrt(1 + u^2 + 0.5*q^2*a^2/m^2) at level n, using u at level n-1/2
#pragma omp simd aligned(u,las:AB)
	for (int i = 0; i < N; i++)
		// this uses las[7]=a^2(n+1/2) - arguably we should use las[6]=a^2(n)
		avgGam[i] = std::sqrt(1.0 + u[1][i] * u[1][i] + u[2][i] * u[2][i] + u[3][i] * u[3][i] + las[7][i]);
	for (int c = 0; c < 4; c++)
#pragma omp simd aligned(avgGam,vel,u:AB)
		for (int i = 0; i < N; i++)
			vel[c][i] = u[c][i] / avgGam[i];
#pragma omp simd aligned(F,vel,avgGam:AB)
	for (int i = 0; i < N; i++)
		avgGam[i] += F[0][i] * vel[1][i] + F[1][i] * vel[2][i] + F[2][i] * vel[3][i];
#pragma omp simd aligned(vel,las,avgGam:AB)
	for (int i = 0; i < N; i++)
		avgGam[i] = 0.5 * (avgGam[i] + std::sqrt(sqr(avgGam[i]) - las[0][i] * vel[1][i] - las[1][i] * vel[2][i] - las[2][i] * vel[3][i]));
}
inline void BundlePusherPGC::avg_gam_2(tw::Float avgGam[N], tw::Float u[4][N], float las[8][N])
{
	// las = [q2m2dth*grad(a^2(n)), q2m2dth*grad(a^2(n+1/2)), q2m2h*a^2(n), q2m2h*a^2(n+1/2)]
	// We are estimating std::sqrt(1 + u^2 + 0.5*q^2*a^2/m^2) at level n+1/2, using u at level n+1/2.
	// This is nontrivial because the *spatial* position where a^2 is evaluated has changed.
	alignas(AB) tw::Float g1[N];
	alignas(AB) tw::Float dA[N];
#pragma omp simd aligned(g1,dA,u,las,avgGam:AB)
	for (int i = 0; i < N; i++)
		g1[i] = std::sqrt(1.0 + u[1][i] * u[1][i] + u[2][i] * u[2][i] + u[3][i] * u[3][i] + las[7][i]);
#pragma omp simd aligned(g1,dA,u,las,avgGam:AB)
	for (int i = 0; i < N; i++)
		dA[i] = (las[3][i] * u[1][i] + las[4][i] * u[2][i] + las[5][i] * u[3][i]) / g1[i];
#pragma omp simd aligned(g1,dA,u,las,avgGam:AB)
	for (int i = 0; i < N; i++)
		avgGam[i] = g1[i] + 0.25 * dA[i] / g1[i] - 0.09375 * dA[i] * dA[i] / (g1[i] * g1[i] * g1[i]);
}
inline void BundlePusherPGC::impulse(tw::Float u[4][N], float F[6][N], float las[8][N], tw::Float avgGam[N])
{
	for (int c = 0; c < 3; c++)
#pragma omp simd aligned(u,F,las,avgGam:AB)
		for (int i = 0; i < N; i++)
			u[c + 1][i] += F[c][i] - 0.25 * las[c][i] / avgGam[i];
}
inline void BundlePusherPGC::rotation1(tw::Float t[3][N], float F[6][N], tw::Float avgGam[N])
{
	for (int c = 0; c < 3; c++)
#pragma omp simd aligned(t,F,avgGam:AB)
		for (int i = 0; i < N; i++)
			t[c][i] = F[c + 3][i] / avgGam[i];
}
inline void BundlePusherPGC::velocity(tw::Float vel[4][N], tw::Float u[4][N], tw::Float avgGam[N])
{
#pragma omp simd aligned(vel,u,avgGam:AB)
	for (int i = 0; i < N; i++)
		u[0][i] = avgGam[i];
	for (int c = 0; c < 4; c++)
#pragma omp simd aligned(vel,u,avgGam:AB)
		for (int i = 0; i < N; i++)
			vel[c][i] = u[c][i] / u[0][i];
}
inline void BundlePusherPGC::load_chi(float chi[N], float number[N], tw::Float avgGam[N])
{
	const float q0 = mov.q0;
	const float m0 = mov.m0;
#pragma omp simd aligned(chi,number,avgGam:AB)
	for (int i = 0; i < N; i++)
		chi[i] = -q0 * q0 * number[i] / (m0 * avgGam[i]);
}

void BundlePusherPGC::Push(tw::Float dts)
{
	avg_gam_1(avgGam, vel, u, F, las);
	impulse(u, F, las, avgGam);
	rotation1(t, F, avgGam);
	rotation2(s, t);
	rotation3(s, t, vel, u);
	impulse(u, F, las, avgGam);
	avg_gam_2(avgGam, u, las);
	velocity(vel, u, avgGam);
	translate(x, vel, dts);
	load_j4(J, number, vel);
	load_chi(chi, number, avgGam);
}


///////////////////////////////////////////
// BOHMIAN PUSHER
//////////////////////////////////////////


inline void BundlePusherBohmian::bohm_velocity(tw::Float vel[4][N], tw::Float u[4][N], float J[4][N])
{
	for (int c = 0; c < 4; c++)
#pragma omp simd aligned(vel,u,J:AB)
		for (int i = 0; i < N; i++)
		{
			const tw::Float vn = J[c][i] / (tw::small_pos + J[0][i]);
			vel[c][i] = 1.5 * vn - 0.5 * u[c][i]; // extrapolate to n+1/2
			u[c][i] = vn;
		}
}

void BundlePusherBohmian::Push(tw::Float dts)
{
	// estimate vel(n+1/2) using u(n-1) and J(n), and update u(n-1) to u(n)
	bohm_velocity(vel, u, J);
	// update x(n) to x(n+1) using vel(n+1/2)
	translate(x, vel, dts);
}


///////////////////////////////////////////
// PHOTON PUSHER
//////////////////////////////////////////


inline void BundlePusherPhoton::velocity(tw::Float vel[4][N], tw::Float u[4][N])
{
	for (int c = 0; c < 4; c++)
#pragma omp simd aligned(vel,u:AB)
		for (int i = 0; i < N; i++)
			vel[c][i] = u[c][i] / u[0][i];
}

void BundlePusherPhoton::Push(tw::Float dts)
{
	velocity(vel,u);
	translate(x,vel,dts);
}
