#include "meta_base.h"
#include "computeTool.h"
#include "bundle.h"
#include "pusher.h"
#include "slicer.h"
#include "tiler.h"
#include "mover.h"

/// @brief Select subcycling for this step based on the field
/// @param F0 the angle tensor *will be rescaled** for the substep duration
/// @return number of subcycles
tw::Int set_subcycles(float F0[6][tw::max_bundle_size])
{
	const tw::Int N = tw::max_bundle_size;
	const tw::Int AB = tw::vec_align_bytes;
	float max2 = 0.0;
	for (int c=0;c<6;c++)
		for (int i=0;i<N;i++)
		{
			if (sqr(F0[c][i])>max2)
				max2 = sqr(F0[c][i]);
		}
	tw::Int subcycles = 1 + tw::Int(sqrt(max2));
	const float subi = 1.0f/subcycles;
	for (int c=0;c<6;c++)
		#pragma omp simd aligned(F0:AB)
		for (int i=0;i<N;i++)
			F0[c][i] *= subi;
	return subcycles;
}

void BundlePusherBoris::Push(tw::Float dts)
{
	int subcycles = set_subcycles(F);
	for (int sub=0;sub<subcycles;sub++)
	{
		impulse(u,F);
		rotation1(t,u,F);
		rotation2(s,t);
		rotation3(s,t,vel,u);
		impulse(u,F);
	}
	velocity(vel,u);
	translate(x,vel,dts);
	load_j4(J,number,vel);
}

void BundlePusherUnitary::Lambda()
{
	// Left multiply the spinor by Lambda (time translation operator)
	copy_spinor(zf,zi);

	copy_spinor(z,zi);
	set_psi24(a,F);
	scalar_mul(a,z);
	add_spinor(zf,z);

	copy_spinor(z,zi);
	left_mul_sig1(z);
	set_psi_1(a,F);
	scalar_mul(a,z);
	add_spinor(zf,z);

	copy_spinor(z,zi);
	left_mul_sig2(z);
	set_psi_2(a,F);
	scalar_mul(a,z);
	add_spinor(zf,z);

	copy_spinor(z,zi);
	left_mul_sig3(z);
	set_psi_3(a,F);
	scalar_mul(a,z);
	add_spinor(zf,z);
}

void BundlePusherUnitary::dilate_timestep(float F0[6][N],float F[6][N],tw::Float u[4][N])
{
	// ds[o][i] is the ratio of proper time step to lab step of order o+1 for particle i
	alignas(AB) tw::Float ds[4][N],E2[N],B2[N],udotE[N],udotB[N],EdotB[N],udotExB[N],tst[N];

	// TODO: is this needed/desirable
	#pragma omp simd aligned(u:AB)
	for (int i=0;i<N;i++)
		u[0][i] = sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i]);

	for (int c=0;c<6;c++)
		#pragma omp simd aligned(F0,F:AB)
		for (int i=0;i<N;i++)
			F[c][i] = F0[c][i];
	
	#pragma omp simd aligned(ds,F,u,E2,B2,udotE,udotB,EdotB,udotExB,tst:AB)
	for (int i=0;i<N;i++)
	{
		udotE[i] = F[0][i]*u[1][i] + F[1][i]*u[2][i] + F[2][i]*u[3][i];
		udotB[i] = F[3][i]*u[1][i] + F[4][i]*u[2][i] + F[5][i]*u[3][i];
		EdotB[i] = F[0][i]*F[3][i] + F[1][i]*F[4][i] + F[2][i]*F[5][i];
		E2[i] = F[0][i]*F[0][i] + F[1][i]*F[1][i] + F[2][i]*F[2][i];
		B2[i] = F[3][i]*F[3][i] + F[4][i]*F[4][i] + F[5][i]*F[5][i];
		udotExB[i] =
			u[1][i]*(F[1][i]*F[5][i] - F[2][i]*F[4][i]) +
			u[2][i]*(F[2][i]*F[3][i] - F[0][i]*F[5][i]) +
			u[3][i]*(F[0][i]*F[4][i] - F[1][i]*F[3][i]);

		ds[0][i] = 1/u[0][i];

		tst[i] = udotE[i]*sqr(ds[0][i]);
		ds[1][i] = ds[0][i] * (1 - tw::Float(tst[i]<0)*tst[i]) / (1 + tw::Float(tst[i]>=0)*tst[i]);
		tst[i] = tw::Float(ds[1][i] < 1.0);
		ds[1][i] = tst[i]*ds[1][i] + (1-tst[i]);

		tst[i] = udotE[i]*sqr(ds[1][i]) + 0.66667*(u[0][i]*E2[i] - udotExB[i])*cub(ds[1][i]);
		ds[2][i] = ds[0][i] * (1 - tw::Float(tst[i]<0)*tst[i]) / (1 + tw::Float(tst[i]>=0)*tst[i]);
		tst[i] = tw::Float(ds[2][i] < 1.0);
		ds[2][i] = tst[i]*ds[2][i] + (1-tst[i]);

		tst[i] = udotE[i]*sqr(ds[2][i]) + 0.66667*(u[0][i]*E2[i] - udotExB[i])*cub(ds[2][i])
			+ 0.125*(3*udotE[i]*(E2[i]-B2[i]) + 2*udotB[i]*EdotB[i])*quad(ds[2][i]);
		ds[3][i] = ds[0][i] * (1 - tw::Float(tst[i]<0)*tst[i]) / (1 + tw::Float(tst[i]>=0)*tst[i]);
		tst[i] = tw::Float(ds[3][i] < 1.0);
		ds[3][i] = tst[i]*ds[3][i] + (1-tst[i]);

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
	for (int sub=0;sub<subcycles;sub++)
	{
		// We rewrite L*z*L^dag as (L*(L*z)^dag)^dag
		dilate_timestep(F0,F,u); // dilation gets buried in F
		to_spinor(u,zi);
		Lambda();
		copy_spinor(zi,zf);
		dagger(zi);
		Lambda();
		dagger(zf);
		to_vector(u,zf,a);
	}
	velocity(vel,u);
	translate(x,vel,dts);
	load_j4(J,number,vel);
}

void BundlePusherPGC::Push(tw::Float dts)
{
	avg_gam_1(avgGam,vel,u,F,las);
	impulse(u,F,las,avgGam);
	rotation1(t,F,avgGam);
	rotation2(s,t);
	rotation3(s,t,vel,u);
	impulse(u,F,las,avgGam);
	avg_gam_2(avgGam,u,las);
	velocity(vel,u,avgGam);
	translate(x,vel,dts);
	load_j4(J,number,vel);
	load_chi(chi,number,avgGam);
}

void BundlePusherBohmian::Push(tw::Float dts)
{
	// estimate vel(n+1/2) using u(n-1) and J(n), and update u(n-1) to u(n)
	bohm_velocity(vel,u,J);
	// update x(n) to x(n+1) using vel(n+1/2)
	translate(x,vel,dts);
}

void BundlePusherPhoton::Push(tw::Float dts)
{
	velocity(vel,u);
	translate(x,vel,dts);
}
