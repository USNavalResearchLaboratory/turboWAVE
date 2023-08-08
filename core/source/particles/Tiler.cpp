#include "meta_base.h"
#include "computeTool.h"
#include "bundle.h"
#include "pusher.h"
#include "slicer.h"
#include "tiler.h"
#include "mover.h"

void BundleTilerEM2D::GatherF(float F[6][N],const float w[3][3][N],const float l[3][3][N],const float qmdth)
{
	// Assumes every particle in bundle is in the same cell and data is packed with Yee fields
	ZeroArray(F,0,5);
	for (tw::Int i=0;i<3;i++)
		for (tw::Int k=0;k<3;k++)
			#pragma omp simd aligned(F,w,l:AB)
			for (tw::Int n=0;n<N;n++)
			{
				F[0][n] += l[i][0][n]*w[k][2][n]*F_tile[i][k][0]*qmdth;
				F[1][n] += w[i][0][n]*w[k][2][n]*F_tile[i][k][1]*qmdth;
				F[2][n] += w[i][0][n]*l[k][2][n]*F_tile[i][k][2]*qmdth;
				F[3][n] += w[i][0][n]*l[k][2][n]*F_tile[i][k][3]*qmdth;
				F[4][n] += l[i][0][n]*l[k][2][n]*F_tile[i][k][4]*qmdth;
				F[5][n] += l[i][0][n]*w[k][2][n]*F_tile[i][k][5]*qmdth;
			}
}

void BundleTilerEM2D::ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti)
{
	// J_tile will be loaded with charge and current
	// if J_tile is nonzero on entry, it must contain the density decomposition, not the current
	// Current deposition from T.Zh. Esirkepov, Comp. Phys. Comm. 135, 144 (2001)
	tw::Int i1,i2,k1,k2,dc[3];
	float xw0[5][3];
	float xw1[5][3];
	float xdw[5][3];
	alignas(AB) float dw[3][3][N];
	float sum;

	// Scalar code to scatter current from inter-cell particles.
	for (tw::Int n=0;n<num;n++)
	{
		if (cellMask[n]!=0.0)
			continue;

		get_cell_displ(dc,n);

		// Build Extended Weights
		for (tw::Int i=0;i<5;i++)
			for (tw::Int j=0;j<3;j++)
			{
				xw0[i][j] = 0.0f;
				xw1[i][j] = 0.0f;
			}
		for (tw::Int i=0;i<3;i++)
			for (tw::Int j=0;j<3;j++)
			{
				xw0[i+1][j] = w0[i][j][n];
				xw1[i+1+dc[j]][j] = w1[i][j][n];
			}
		for (tw::Int i=0;i<5;i++)
			for (tw::Int j=0;j<3;j++)
				xdw[i][j] = xw1[i][j] - xw0[i][j];

		// Charge Deposition (as in coulombs in the cell)
		for (tw::Int i=0;i<5;i++)
			for (tw::Int k=0;k<5;k++)
				J_tile[i][k][0] += J[0][n]*xw1[i][0]*xw1[k][2];

		// Optimize Current Deposition Loops
		i1 = dc[0]<0 ? 0 : 1;
		i2 = dc[0]>0 ? 4 : 3;
		k1 = dc[2]<0 ? 0 : 1;
		k2 = dc[2]>0 ? 4 : 3;

		// Current Deposition (as in amps through the wall)
		for (tw::Int i=1;i<=i2;i++)
			for (tw::Int k=k1;k<=k2;k++)
				J_tile[i][k][1] += J[0][n]*dti*xdw[i][0]*(xw0[k][2]+0.5f*xdw[k][2]);
		for (tw::Int i=i1;i<=i2;i++)
			for (tw::Int k=k1;k<=k2;k++)
				J_tile[i][k][2] += J[2][n]*(xw0[i][0]*xw0[k][2]+0.5f*xdw[i][0]*xw0[k][2]+0.5f*xw0[i][0]*xdw[k][2]+0.3333333f*xdw[i][0]*xdw[k][2]);
		for (tw::Int i=i1;i<=i2;i++)
			for (tw::Int k=1;k<=k2;k++)
				J_tile[i][k][3] += J[0][n]*dti*xdw[k][2]*(xw0[i][0]+0.5f*xdw[i][0]);
	}

	// Vector code to scatter current from intra-cell particles
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			#pragma omp simd aligned(dw,w0,w1:AB)
			for (tw::Int n=0;n<N;n++)
				dw[i][j][n] = w1[i][j][n] - w0[i][j][n];

	for (tw::Int i=0;i<3;i++)
		for (tw::Int k=0;k<3;k++)
		{
			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,w1,cellMask:AB)
			for (tw::Int n=0;n<N;n++)
				sum += J[0][n]*cellMask[n]*w1[i][0][n]*w1[k][2][n];
			J_tile[i+1][k+1][0] += sum;

			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
			for (tw::Int n=0;n<N;n++)
				sum += J[0][n]*dti*cellMask[n]*dw[i][0][n]*(w0[k][2][n]+0.5f*dw[k][2][n]);
			J_tile[i+1][k+1][1] += sum;

			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
			for (tw::Int n=0;n<N;n++)
				sum += J[2][n]*cellMask[n]*(w0[i][0][n]*w0[k][2][n]+0.5f*dw[i][0][n]*w0[k][2][n]+0.5f*w0[i][0][n]*dw[k][2][n]+0.3333333f*dw[i][0][n]*dw[k][2][n]);
			J_tile[i+1][k+1][2] += sum;

			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
			for (tw::Int n=0;n<N;n++)
				sum += J[0][n]*dti*cellMask[n]*dw[k][2][n]*(w0[i][0][n]+0.5f*dw[i][0][n]);
			J_tile[i+1][k+1][3] += sum;
		}

	// Integrate to get current from density decomposition
	// The relationship is J_i = j_i+1 + W_i
	for (tw::Int k=0;k<5;k++)
	{
		sum = J_tile[4][k][1];
		sum += J_tile[3][k][1]; J_tile[3][k][1] = sum;
		sum += J_tile[2][k][1]; J_tile[2][k][1] = sum;
		sum += J_tile[1][k][1]; J_tile[1][k][1] = sum;
		J_tile[0][k][1] = 0.0;

		sum = J_tile[k][4][3];
		sum += J_tile[k][3][3]; J_tile[k][3][3] = sum;
		sum += J_tile[k][2][3]; J_tile[k][2][3] = sum;
		sum += J_tile[k][1][3]; J_tile[k][1][3] = sum;
		J_tile[k][0][3] = 0.0;
	}
}


void BundleTilerEM3D::GatherF(float F[6][N],const float w[3][3][N],const float l[3][3][N],const float qmdth)
{
	// Assumes every particle in bundle is in the same cell and data is packed with Yee fields
	ZeroArray(F,0,5);
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
				#pragma omp simd aligned(F,w,l:AB)
				for (tw::Int n=0;n<N;n++)
				{
					F[0][n] += l[i][0][n]*w[j][1][n]*w[k][2][n]*F_tile[i][j][k][0]*qmdth;
					F[1][n] += w[i][0][n]*l[j][1][n]*w[k][2][n]*F_tile[i][j][k][1]*qmdth;
					F[2][n] += w[i][0][n]*w[j][1][n]*l[k][2][n]*F_tile[i][j][k][2]*qmdth;
					F[3][n] += w[i][0][n]*l[j][1][n]*l[k][2][n]*F_tile[i][j][k][3]*qmdth;
					F[4][n] += l[i][0][n]*w[j][1][n]*l[k][2][n]*F_tile[i][j][k][4]*qmdth;
					F[5][n] += l[i][0][n]*l[j][1][n]*w[k][2][n]*F_tile[i][j][k][5]*qmdth;
				}
}

void BundleTilerEM3D::ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti)
{
	// J_tile will be loaded with charge and current
	// if J_tile is nonzero on entry, it must contain the density decomposition, not the current
	// Current deposition from T.Zh. Esirkepov, Comp. Phys. Comm. 135, 144 (2001)
	tw::Int i1,i2,j1,j2,k1,k2,dc[3];
	float xw0[5][3];
	float xw1[5][3];
	float xdw[5][3];
	alignas(AB) float dw[3][3][N];
	float sum;

	// Scalar code to scatter current from inter-cell particles.
	for (tw::Int n=0;n<num;n++)
	{
		if (cellMask[n]!=0.0)
			continue;

		get_cell_displ(dc,n);

		// Build Extended Weights
		for (tw::Int i=0;i<5;i++)
			for (tw::Int j=0;j<3;j++)
			{
				xw0[i][j] = 0.0f;
				xw1[i][j] = 0.0f;
			}
		for (tw::Int i=0;i<3;i++)
			for (tw::Int j=0;j<3;j++)
			{
				xw0[i+1][j] = w0[i][j][n];
				xw1[i+1+dc[j]][j] = w1[i][j][n];
			}
		for (tw::Int i=0;i<5;i++)
			for (tw::Int j=0;j<3;j++)
				xdw[i][j] = xw1[i][j] - xw0[i][j];

		// Charge Deposition (as in coulombs in the cell)
		for (tw::Int i=0;i<5;i++)
			for (tw::Int j=0;j<5;j++)
				for (tw::Int k=0;k<5;k++)
					J_tile[i][j][k][0] += J[0][n]*xw1[i][0]*xw1[j][1]*xw1[k][2];

		// Optimize Current Deposition Loops
		i1 = dc[0]<0 ? 0 : 1;
		i2 = dc[0]>0 ? 4 : 3;
		j1 = dc[1]<0 ? 0 : 1;
		j2 = dc[1]>0 ? 4 : 3;
		k1 = dc[2]<0 ? 0 : 1;
		k2 = dc[2]>0 ? 4 : 3;

		// Current Deposition (as in amps through the wall)

		for (tw::Int i=1;i<=i2;i++)
			for (tw::Int j=j1;j<=j2;j++)
				for (tw::Int k=k1;k<=k2;k++)
					J_tile[i][j][k][1] += J[0][n]*dti*xdw[i][0]*(xw0[j][1]*xw0[k][2]+
						0.5f*xdw[j][1]*xw0[k][2]+0.5f*xw0[j][1]*xdw[k][2]+0.3333333f*xdw[j][1]*xdw[k][2]);
		for (tw::Int i=i1;i<=i2;i++)
			for (tw::Int j=1;j<=j2;j++)
				for (tw::Int k=k1;k<=k2;k++)
					J_tile[i][j][k][2] += J[0][n]*dti*xdw[j][1]*(xw0[i][0]*xw0[k][2]+
						0.5f*xdw[i][0]*xw0[k][2]+0.5f*xw0[i][0]*xdw[k][2]+0.3333333f*xdw[i][0]*xdw[k][2]);
		for (tw::Int i=i1;i<=i2;i++)
			for (tw::Int j=j1;j<=j2;j++)
				for (tw::Int k=1;k<=k2;k++)
					J_tile[i][j][k][3] += J[0][n]*dti*xdw[k][2]*(xw0[i][0]*xw0[j][1]+
						0.5f*xdw[i][0]*xw0[j][1]+0.5f*xw0[i][0]*xdw[j][1]+0.3333333f*xdw[i][0]*xdw[j][1]);
	}

	// Vector code to scatter current from intra-cell particles
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			#pragma omp simd aligned(dw,w0,w1:AB)
			for (tw::Int n=0;n<N;n++)
				dw[i][j][n] = w1[i][j][n] - w0[i][j][n];

	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
			{
				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,w1,cellMask:AB)
				for (tw::Int n=0;n<N;n++)
					sum += J[0][n]*cellMask[n]*w1[i][0][n]*w1[j][1][n]*w1[k][2][n];
				J_tile[i+1][j+1][k+1][0] += sum;

				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
				for (tw::Int n=0;n<N;n++)
					sum += J[0][n]*dti*cellMask[n]*dw[i][0][n]*(w0[j][1][n]*w0[k][2][n]+
						0.5f*dw[j][1][n]*w0[k][2][n]+0.5f*w0[j][1][n]*dw[k][2][n]+0.3333333f*dw[j][1][n]*dw[k][2][n]);
				J_tile[i+1][j+1][k+1][1] += sum;

				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
				for (tw::Int n=0;n<N;n++)
					sum += J[0][n]*dti*cellMask[n]*dw[j][1][n]*(w0[i][0][n]*w0[k][2][n]+
						0.5f*dw[i][0][n]*w0[k][2][n]+0.5f*w0[i][0][n]*dw[k][2][n]+0.3333333f*dw[i][0][n]*dw[k][2][n]);
				J_tile[i+1][j+1][k+1][2] += sum;

				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
				for (tw::Int n=0;n<N;n++)
					sum += J[0][n]*dti*cellMask[n]*dw[k][2][n]*(w0[i][0][n]*w0[j][1][n]+
						0.5f*dw[i][0][n]*w0[j][1][n]+0.5f*w0[i][0][n]*dw[j][1][n]+0.3333333f*dw[i][0][n]*dw[j][1][n]);
				J_tile[i+1][j+1][k+1][3] += sum;
			}

	// Integrate to get current from density decomposition
	// The relationship is J_i = j_i+1 + W_i
	for (tw::Int j=0;j<5;j++)
		for (tw::Int k=0;k<5;k++)
		{
			sum = J_tile[4][j][k][1];
			sum += J_tile[3][j][k][1]; J_tile[3][j][k][1] = sum;
			sum += J_tile[2][j][k][1]; J_tile[2][j][k][1] = sum;
			sum += J_tile[1][j][k][1]; J_tile[1][j][k][1] = sum;
			J_tile[0][j][k][1] = 0.0;

			sum = J_tile[j][4][k][2];
			sum += J_tile[j][3][k][2]; J_tile[j][3][k][2] = sum;
			sum += J_tile[j][2][k][2]; J_tile[j][2][k][2] = sum;
			sum += J_tile[j][1][k][2]; J_tile[j][1][k][2] = sum;
			J_tile[j][0][k][2] = 0.0;

			sum = J_tile[j][k][4][3];
			sum += J_tile[j][k][3][3]; J_tile[j][k][3][3] = sum;
			sum += J_tile[j][k][2][3]; J_tile[j][k][2][3] = sum;
			sum += J_tile[j][k][1][3]; J_tile[j][k][1][3] = sum;
			J_tile[j][k][0][3] = 0.0;
		}
}

void BundleTilerPGC2D::GatherLaser(float las[8][N],const float w[3][3][N],const float q2m2dth)
{
	alignas(AB) float factorNow[N];
	ZeroArray(las,0,7);
	for (tw::Int i=0;i<3;i++)
		for (tw::Int k=0;k<3;k++)
		{
			#pragma omp simd aligned(w,factorNow:AB)
			for (tw::Int n=0;n<N;n++)
				factorNow[n] = w[i][0][n]*w[k][2][n];
			#pragma omp simd aligned(las,factorNow:AB)
			for (tw::Int n=0;n<N;n++)
			{
				las[0][n] += factorNow[n]*las_tile[i][k][0]*q2m2dth;
				las[1][n] += factorNow[n]*las_tile[i][k][1]*q2m2dth;
				las[2][n] += factorNow[n]*las_tile[i][k][2]*q2m2dth;
				las[3][n] += factorNow[n]*las_tile[i][k][3]*q2m2dth;
				las[4][n] += factorNow[n]*las_tile[i][k][4]*q2m2dth;
				las[5][n] += factorNow[n]*las_tile[i][k][5]*q2m2dth;
				las[6][n] += factorNow[n]*las_tile[i][k][6]*q2m2dth;
				las[7][n] += factorNow[n]*las_tile[i][k][7]*q2m2dth;
			}
		}
}

void BundleTilerPGC2D::ScatterChi(const float chi[N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N])
{
	tw::Int dc[3];
	float xw0[5][3];
	float xw1[5][3];
	float sum;

	// Scalar code to scatter current from inter-cell particles.
	for (tw::Int n=0;n<num;n++)
	{
		if (cellMask[n]!=0.0)
			continue;

		get_cell_displ(dc,n);

		for (tw::Int i=0;i<5;i++)
			for (tw::Int j=0;j<3;j++)
			{
				xw0[i][j] = 0.0f;
				xw1[i][j] = 0.0f;
			}
		for (tw::Int i=0;i<3;i++)
			for (tw::Int j=0;j<3;j++)
			{
				xw0[i+1][j] = w0[i][j][n];
				xw1[i+1+dc[j]][j] = w1[i][j][n];
			}

		for (tw::Int i=0;i<5;i++)
			for (tw::Int k=0;k<5;k++)
				chi_tile[i][k] += 0.5f*chi[n]*(xw0[i][0]*xw0[k][2] + xw1[i][0]*xw1[k][2]);
	}

	// Vector code to scatter current from intra-cell particles
	for (tw::Int i=0;i<3;i++)
		for (tw::Int k=0;k<3;k++)
		{
			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(chi,w0,w1,cellMask:AB)
			for (tw::Int n=0;n<N;n++)
				sum += 0.5f*chi[n]*cellMask[n]*(w0[i][0][n]*w0[k][2][n] + w1[i][0][n]*w1[k][2][n]);
			chi_tile[i+1][k+1] += sum;
		}
}

void BundleTilerPGC3D::GatherLaser(float las[8][N],const float w[3][3][N],const float q2m2dth)
{
	alignas(AB) float factorNow[N];
	ZeroArray(las,0,7);
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
			{
				#pragma omp simd aligned(w,factorNow:AB)
				for (tw::Int n=0;n<N;n++)
					factorNow[n] = w[i][0][n]*w[j][1][n]*w[k][2][n];
				#pragma omp simd aligned(las,factorNow:AB)
				for (tw::Int n=0;n<N;n++)
				{
					las[0][n] += factorNow[n]*las_tile[i][j][k][0]*q2m2dth;
					las[1][n] += factorNow[n]*las_tile[i][j][k][1]*q2m2dth;
					las[2][n] += factorNow[n]*las_tile[i][j][k][2]*q2m2dth;
					las[3][n] += factorNow[n]*las_tile[i][j][k][3]*q2m2dth;
					las[4][n] += factorNow[n]*las_tile[i][j][k][4]*q2m2dth;
					las[5][n] += factorNow[n]*las_tile[i][j][k][5]*q2m2dth;
					las[6][n] += factorNow[n]*las_tile[i][j][k][6]*q2m2dth;
					las[7][n] += factorNow[n]*las_tile[i][j][k][7]*q2m2dth;
				}
			}
}

void BundleTilerPGC3D::ScatterChi(const float chi[N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N])
{
	tw::Int dc[3];
	float xw0[5][3];
	float xw1[5][3];
	float sum;

	// Scalar code to scatter current from inter-cell particles.
	for (tw::Int n=0;n<num;n++)
	{
		if (cellMask[n]!=0.0)
			continue;

		get_cell_displ(dc,n);

		for (tw::Int i=0;i<5;i++)
			for (tw::Int j=0;j<3;j++)
			{
				xw0[i][j] = 0.0f;
				xw1[i][j] = 0.0f;
			}
		for (tw::Int i=0;i<3;i++)
			for (tw::Int j=0;j<3;j++)
			{
				xw0[i+1][j] = w0[i][j][n];
				xw1[i+1+dc[j]][j] = w1[i][j][n];
			}

		for (tw::Int i=0;i<5;i++)
			for (tw::Int j=0;j<5;j++)
				for (tw::Int k=0;k<5;k++)
					chi_tile[i][j][k] += 0.5f*chi[n]*(xw0[i][0]*xw0[j][1]*xw0[k][2] + xw1[i][0]*xw1[j][1]*xw1[k][2]);
	}

	// Vector code to scatter current from intra-cell particles
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
			{
				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(chi,w0,w1,cellMask:AB)
				for (tw::Int n=0;n<N;n++)
					sum += 0.5f*chi[n]*cellMask[n]*(w0[i][0][n]*w0[j][1][n]*w0[k][2][n] + w1[i][0][n]*w1[j][1][n]*w1[k][2][n]);
				chi_tile[i+1][j+1][k+1] += sum;
			}
}

void BundleTilerBohmian2D::GatherJ4(float J[4][N],const float w[3][3][N])
{
	ZeroArray(J,0,3);
	for (tw::Int i=0;i<3;i++)
		for (tw::Int k=0;k<3;k++)
			#pragma omp simd aligned(J,w:AB)
			for (tw::Int n=0;n<N;n++)
			{
				J[0][n] += w[i][0][n]*w[k][2][n]*tile[i][k][0];
				J[1][n] += w[i][0][n]*w[k][2][n]*tile[i][k][1];
				J[2][n] += w[i][0][n]*w[k][2][n]*tile[i][k][2];
				J[3][n] += w[i][0][n]*w[k][2][n]*tile[i][k][3];
			}
}

void BundleTilerBohmian3D::GatherJ4(float J[4][N],const float w[3][3][N])
{
	ZeroArray(J,0,3);
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
				#pragma omp simd aligned(J,w:AB)
				for (tw::Int n=0;n<N;n++)
				{
					J[0][n] += w[i][0][n]*w[j][1][n]*w[k][2][n]*tile[i][j][k][0];
					J[1][n] += w[i][0][n]*w[j][1][n]*w[k][2][n]*tile[i][j][k][1];
					J[2][n] += w[i][0][n]*w[j][1][n]*w[k][2][n]*tile[i][j][k][2];
					J[3][n] += w[i][0][n]*w[j][1][n]*w[k][2][n]*tile[i][j][k][3];
				}
}
