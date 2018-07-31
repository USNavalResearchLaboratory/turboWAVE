#include "sim.h"
#include "particles.h"
#include "fieldSolve.h"
#include "laserSolve.h"

void Species::GetSubarrayBounds(std::vector<ParticleRef>& sorted,tw::Int low[4],tw::Int high[4],tw::Int layers)
{
	// Assumes particles are sorted in increasing memory order
	// This depends in turn on the cell encoding respecting memory order

	DecodeCell(sorted.front().cell,low);
	DecodeCell(sorted.back().cell,high);
	for (int i=1;i<=3;i++)
		if (low[i]>high[i])
			std::swap(low[i],high[i]);
	// Assume z-packing and maximal loading
	if (low[1]!=high[1])
	{
		low[2] = 1;
		high[2] = Dim(2);
		low[3] = 1;
		high[3] = Dim(3);
	}
	if (low[2]!=high[2])
	{
		low[3] = 1;
		high[3] = Dim(3);
	}
	// Expand subarray to allow for motion and particle shape
	// For gather this induces 1 ghost cell layers
	// For scatter this induces 2 ghost cell layers
	for (int i=1;i<=3;i++)
	{
		if (Dim(i)>1)
		{
			low[i] -= layers;
			high[i] += layers;
		}
		if (low[i]<owner->N0(i))
			low[i] = owner->N0(i);
		if (high[i]>owner->N1(i))
			high[i] = owner->N1(i);
	}
}

template <class BundleType>
void Species::Push()
{
	if (particle.size()==0)
		return;

	const tw::Int min_particles_per_task = 256;
	const tw::Int num = particle.size();
	const tw::Int max_tasks = 1 + num / min_particles_per_task;
	const tw::Int preferred_tasks = 4*tw::GetOMPMaxThreads();
	const tw::Int num_tasks = preferred_tasks > max_tasks ? max_tasks : preferred_tasks;

	#pragma omp parallel for
	for (tw::Int tsk=0;tsk<num_tasks;tsk++)
	{
		#pragma omp task firstprivate(tsk)
		{
			tw::Int i,first,last,next,low[4],high[4];
			std::vector<ParticleRef> map;
			BundleType b;
			tw::GetOMPTaskLoopRange(tsk,num,num_tasks,&first,&last);
			map.resize(last-first+1);
			for (i=first;i<=last;i++)
				map[i-first] = ParticleRef(i,particle[i]);
			std::sort(map.begin(),map.end());
			GetSubarrayBounds(map,low,high,1);
			b.LoadFieldSlice(this,low,high,ignorable);
			GetSubarrayBounds(map,low,high,2);
			b.InitSourceSlice(this,low,high,ignorable);
			for (i=first;i<=last;i++)
			{
				b.Append(particle[map[i-first].idx]);
				next = i==last ? i : i+1;
				if (i==last || b.Complete(particle[map[next-first].idx]))
				{
					b.Push(this);
					b.CopyBack(this);
					b.Reset();
				}
			}
			b.DepositSourceSlice(this);
		}
	}
}

void Species::DispatchPush()
{
	if (qo_j4)
	{
		Push<ParticleBundleBohmian>();
		return;
	}
	if (laser)
	{
		Push<ParticleBundlePGC>();
		return;
	}

	if (ignorable[2])
		Push<ParticleBundle2D>();
	else
		Push<ParticleBundle3D>();
}


///////////////////////////////
//                           //
//     BUNDLE OPERATIONS     //
//                           //
///////////////////////////////


void ParticleBundleBohmian::Push(Species *owner)
{
	// The empty particles in the bundle can be operated on harmlessly
	// except in gather scatter.
	// N = particles in a full bundle
	// num = particles in this bundle
	Grid *theGrid = owner->owner;
	const tw::Float q0 = owner->charge;
	const tw::Float m0 = owner->restMass;
	const tw::Float dth = theGrid->dth;
	const tw::Float dt = theGrid->dt;
	const tw::Float k[3] = { dxi(*theGrid) , dyi(*theGrid) , dzi(*theGrid) };
	const float dti = 1.0/dt;

	PadBundle();

	cell0 = cell[0];
	theGrid->DecodeCell(cell0,&ijk0[0],&ijk0[1],&ijk0[2]);
	theGrid->GetWeights(w0,x);
	LoadTile();
	GatherJ4(J,w0);

	// estimate vel(n+1/2) using p(n-1) and J(n), and update p(n-1) to p(n)
	bohm_velocity(vel,p,J);
	// update x(n) to x(n+1) using vel(n+1/2)
	translate(x,vel,k,dt);
	theGrid->MinimizePrimitive(cell,ijk,x,domainMask);
}

void ParticleBundleBohmian::GatherJ4(float J[4][N],const float w[3][3][N])
{
	tw::Int i,j,k,n;
	ZeroArray(J,0,3);
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			for (k=0;k<3;k++)
				#pragma omp simd aligned(J,w:AB)
				for (n=0;n<N;n++)
				{
					J[0][n] += w[i][0][n]*w[j][1][n]*w[k][2][n]*tile[i][j][k][0];
					J[1][n] += w[i][0][n]*w[j][1][n]*w[k][2][n]*tile[i][j][k][1];
					J[2][n] += w[i][0][n]*w[j][1][n]*w[k][2][n]*tile[i][j][k][2];
					J[3][n] += w[i][0][n]*w[j][1][n]*w[k][2][n]*tile[i][j][k][3];
				}
}


void ParticleBundle2D::Push(Species *owner)
{
	// The empty particles in the bundle can be operated on harmlessly
	// except in gather scatter.
	// N = particles in a full bundle
	// num = particles in this bundle
	Grid *theGrid = owner->owner;
	const tw::Float q0 = owner->charge;
	const tw::Float m0 = owner->restMass;
	const tw::Float dth = theGrid->dth;
	const tw::Float dt = theGrid->dt;
	const tw::Float k[3] = { dxi(*theGrid) , dyi(*theGrid) , dzi(*theGrid) };
	const float dti = 1.0/dt;

	PadBundle();

	cell0 = cell[0];
	theGrid->DecodeCell(cell0,&ijk0[0],&ijk0[1],&ijk0[2]);
	theGrid->GetWeights(w0,x);
	theGrid->GetWallWeights(l0,x);
	LoadTile();
	GatherF(F,w0,l0);

	impulse(p,F,q0,dth);
	rotation1(t,p,F,q0,m0,dth);
	rotation2(s,t);
	rotation3(s,t,vel,p);
	impulse(p,F,q0,dth);
	velocity(vel,p,m0);
	translate(x,vel,k,dt);
	load_j4(J,number,vel,k,q0);

	theGrid->MinimizePrimitive(cell,ijk,x,domainMask);
	theGrid->GetWeights(w1,x);
	set_cell_mask(cellMask,cell0,cell);
	ResetXTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreXTile();
}

void ParticleBundle2D::GatherF(float F[6][N],const float w[3][3][N],const float l[3][3][N])
{
	// Assumes every particle in bundle is in the same cell and data is packed with Yee fields
	tw::Int i,k,n;
	ZeroArray(F,0,5);
	for (i=0;i<3;i++)
		for (k=0;k<3;k++)
			#pragma omp simd aligned(F,w,l:AB)
			for (n=0;n<N;n++)
			{
				F[0][n] += l[i][0][n]*w[k][2][n]*tile[i][k][0];
				F[1][n] += w[i][0][n]*w[k][2][n]*tile[i][k][1];
				F[2][n] += w[i][0][n]*l[k][2][n]*tile[i][k][2];
				F[3][n] += w[i][0][n]*l[k][2][n]*tile[i][k][3];
				F[4][n] += l[i][0][n]*l[k][2][n]*tile[i][k][4];
				F[5][n] += l[i][0][n]*w[k][2][n]*tile[i][k][5];
			}
}

void ParticleBundle2D::ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti)
{
	// xtile will be loaded with charge and current
	// if xtile is nonzero on entry, it must contain the density decomposition, not the current
	// Current deposition from T.Zh. Esirkepov, Comp. Phys. Comm. 135, 144 (2001)
	tw::Int i,j,k,i1,i2,k1,k2,dc[3],n;
	float xw0[5][3];
	float xw1[5][3];
	float xdw[5][3];
	alignas(AB) float dw[3][3][N];
	float sum;

	// Scalar code to scatter current from inter-cell particles.
	for (n=0;n<num;n++)
	{
		if (cellMask[n]==0.0)
		{
			get_cell_displ(dc,n);

			// Build Extended Weights
			for (i=0;i<5;i++)
				for (j=0;j<3;j++)
				{
					xw0[i][j] = 0.0f;
					xw1[i][j] = 0.0f;
				}
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
				{
					xw0[i+1][j] = w0[i][j][n];
					xw1[i+1+dc[j]][j] = w1[i][j][n];
				}
			for (i=0;i<5;i++)
				for (j=0;j<3;j++)
					xdw[i][j] = xw1[i][j] - xw0[i][j];

			// Charge Deposition (as in coulombs in the cell)
			for (i=0;i<5;i++)
				for (k=0;k<5;k++)
					xtile[i][k][0] += J[0][n]*xw1[i][0]*xw1[k][2];

			// Optimize Current Deposition Loops
			i1 = dc[0]<0 ? 0 : 1;
			i2 = dc[0]>0 ? 4 : 3;
			k1 = dc[2]<0 ? 0 : 1;
			k2 = dc[2]>0 ? 4 : 3;

			// Current Deposition (as in amps through the wall)
			for (i=1;i<=i2;i++)
				for (k=k1;k<=k2;k++)
					xtile[i][k][1] += J[0][n]*dti*xdw[i][0]*(xw0[k][2]+0.5f*xdw[k][2]);
			for (i=i1;i<=i2;i++)
				for (k=k1;k<=k2;k++)
					xtile[i][k][2] += J[2][n]*(xw0[i][0]*xw0[k][2]+0.5f*xdw[i][0]*xw0[k][2]+0.5f*xw0[i][0]*xdw[k][2]+0.3333333f*xdw[i][0]*xdw[k][2]);
			for (i=i1;i<=i2;i++)
				for (k=1;k<=k2;k++)
					xtile[i][k][3] += J[0][n]*dti*xdw[k][2]*(xw0[i][0]+0.5f*xdw[i][0]);
		}
	}

	// Vector code to scatter current from intra-cell particles
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			#pragma omp simd aligned(dw,w0,w1:AB)
			for (n=0;n<N;n++)
				dw[i][j][n] = w1[i][j][n] - w0[i][j][n];

	for (i=0;i<3;i++)
		for (k=0;k<3;k++)
		{
			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,w1,cellMask:AB)
			for (n=0;n<N;n++)
				sum += J[0][n]*cellMask[n]*w1[i][0][n]*w1[k][2][n];
			xtile[i+1][k+1][0] += sum;

			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
			for (n=0;n<N;n++)
				sum += J[0][n]*dti*cellMask[n]*dw[i][0][n]*(w0[k][2][n]+0.5f*dw[k][2][n]);
			xtile[i+1][k+1][1] += sum;

			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
			for (n=0;n<N;n++)
				sum += J[2][n]*cellMask[n]*(w0[i][0][n]*w0[k][2][n]+0.5f*dw[i][0][n]*w0[k][2][n]+0.5f*w0[i][0][n]*dw[k][2][n]+0.3333333f*dw[i][0][n]*dw[k][2][n]);
			xtile[i+1][k+1][2] += sum;

			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
			for (n=0;n<N;n++)
				sum += J[0][n]*dti*cellMask[n]*dw[k][2][n]*(w0[i][0][n]+0.5f*dw[i][0][n]);
			xtile[i+1][k+1][3] += sum;
		}

	// Integrate to get current from density decomposition
	// The relationship is J_i = j_i+1 + W_i
	for (k=0;k<5;k++)
	{
		sum = xtile[4][k][1];
		sum += xtile[3][k][1]; xtile[3][k][1] = sum;
		sum += xtile[2][k][1]; xtile[2][k][1] = sum;
		sum += xtile[1][k][1]; xtile[1][k][1] = sum;
		xtile[0][k][1] = 0.0;

		sum = xtile[k][4][3];
		sum += xtile[k][3][3]; xtile[k][3][3] = sum;
		sum += xtile[k][2][3]; xtile[k][2][3] = sum;
		sum += xtile[k][1][3]; xtile[k][1][3] = sum;
		xtile[k][0][3] = 0.0;
	}
}

void ParticleBundle3D::Push(Species *owner)
{
	// The empty particles in the bundle can be operated on harmlessly
	// except in gather scatter.
	// N = particles in a full bundle
	// num = particles in this bundle
	Grid *theGrid = owner->owner;
	const tw::Float q0 = owner->charge;
	const tw::Float m0 = owner->restMass;
	const tw::Float dth = theGrid->dth;
	const tw::Float dt = theGrid->dt;
	const tw::Float k[3] = { dxi(*theGrid) , dyi(*theGrid) , dzi(*theGrid) };
	const float dti = 1.0/dt;

	PadBundle();

	cell0 = cell[0];
	theGrid->DecodeCell(cell0,&ijk0[0],&ijk0[1],&ijk0[2]);
	theGrid->GetWeights(w0,x);
	theGrid->GetWallWeights(l0,x);
	LoadTile();
	GatherF(F,w0,l0);

	impulse(p,F,q0,dth);
	rotation1(t,p,F,q0,m0,dth);
	rotation2(s,t);
	rotation3(s,t,vel,p);
	impulse(p,F,q0,dth);
	velocity(vel,p,m0);
	translate(x,vel,k,dt);
	load_j4(J,number,vel,k,q0);

	theGrid->MinimizePrimitive(cell,ijk,x,domainMask);
	theGrid->GetWeights(w1,x);
	set_cell_mask(cellMask,cell0,cell);
	ResetXTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreXTile();
}

void ParticleBundle3D::GatherF(float F[6][N],const float w[3][3][N],const float l[3][3][N])
{
	// Assumes every particle in bundle is in the same cell and data is packed with Yee fields
	tw::Int i,j,k,n;
	ZeroArray(F,0,5);
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			for (k=0;k<3;k++)
				#pragma omp simd aligned(F,w,l:AB)
				for (n=0;n<N;n++)
				{
					F[0][n] += l[i][0][n]*w[j][1][n]*w[k][2][n]*tile[i][j][k][0];
					F[1][n] += w[i][0][n]*l[j][1][n]*w[k][2][n]*tile[i][j][k][1];
					F[2][n] += w[i][0][n]*w[j][1][n]*l[k][2][n]*tile[i][j][k][2];
					F[3][n] += w[i][0][n]*l[j][1][n]*l[k][2][n]*tile[i][j][k][3];
					F[4][n] += l[i][0][n]*w[j][1][n]*l[k][2][n]*tile[i][j][k][4];
					F[5][n] += l[i][0][n]*l[j][1][n]*w[k][2][n]*tile[i][j][k][5];
				}
}

void ParticleBundle3D::ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti)
{
	// xtile will be loaded with charge and current
	// if xtile is nonzero on entry, it must contain the density decomposition, not the current
	// Current deposition from T.Zh. Esirkepov, Comp. Phys. Comm. 135, 144 (2001)
	tw::Int i,j,k,i1,i2,j1,j2,k1,k2,dc[3],n;
	float xw0[5][3];
	float xw1[5][3];
	float xdw[5][3];
	alignas(AB) float dw[3][3][N];
	float sum;

	// Scalar code to scatter current from inter-cell particles.
	for (n=0;n<num;n++)
	{
		if (cellMask[n]==0.0)
		{
			get_cell_displ(dc,n);

			// Build Extended Weights
			for (i=0;i<5;i++)
				for (j=0;j<3;j++)
				{
					xw0[i][j] = 0.0f;
					xw1[i][j] = 0.0f;
				}
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
				{
					xw0[i+1][j] = w0[i][j][n];
					xw1[i+1+dc[j]][j] = w1[i][j][n];
				}
			for (i=0;i<5;i++)
				for (j=0;j<3;j++)
					xdw[i][j] = xw1[i][j] - xw0[i][j];

			// Charge Deposition (as in coulombs in the cell)
			for (i=0;i<5;i++)
				for (j=0;j<5;j++)
					for (k=0;k<5;k++)
						xtile[i][j][k][0] += J[0][n]*xw1[i][0]*xw1[j][1]*xw1[k][2];

			// Optimize Current Deposition Loops
			i1 = dc[0]<0 ? 0 : 1;
			i2 = dc[0]>0 ? 4 : 3;
			j1 = dc[1]<0 ? 0 : 1;
			j2 = dc[1]>0 ? 4 : 3;
			k1 = dc[2]<0 ? 0 : 1;
			k2 = dc[2]>0 ? 4 : 3;

			// Current Deposition (as in amps through the wall)

			for (i=1;i<=i2;i++)
				for (j=j1;j<=j2;j++)
					for (k=k1;k<=k2;k++)
						xtile[i][j][k][1] += J[0][n]*dti*xdw[i][0]*(xw0[j][1]*xw0[k][2]+
							0.5f*xdw[j][1]*xw0[k][2]+0.5f*xw0[j][1]*xdw[k][2]+0.3333333f*xdw[j][1]*xdw[k][2]);
			for (i=i1;i<=i2;i++)
				for (j=1;j<=j2;j++)
					for (k=k1;k<=k2;k++)
						xtile[i][j][k][2] += J[0][n]*dti*xdw[j][1]*(xw0[i][0]*xw0[k][2]+
							0.5f*xdw[i][0]*xw0[k][2]+0.5f*xw0[i][0]*xdw[k][2]+0.3333333f*xdw[i][0]*xdw[k][2]);
			for (i=i1;i<=i2;i++)
				for (j=j1;j<=j2;j++)
					for (k=1;k<=k2;k++)
						xtile[i][j][k][3] += J[0][n]*dti*xdw[k][2]*(xw0[i][0]*xw0[j][1]+
							0.5f*xdw[i][0]*xw0[j][1]+0.5f*xw0[i][0]*xdw[j][1]+0.3333333f*xdw[i][0]*xdw[j][1]);
		}
	}

	// Vector code to scatter current from intra-cell particles
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			#pragma omp simd aligned(dw,w0,w1:AB)
			for (n=0;n<N;n++)
				dw[i][j][n] = w1[i][j][n] - w0[i][j][n];

	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			for (k=0;k<3;k++)
			{
				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,w1,cellMask:AB)
				for (n=0;n<N;n++)
					sum += J[0][n]*cellMask[n]*w1[i][0][n]*w1[j][1][n]*w1[k][2][n];
				xtile[i+1][j+1][k+1][0] += sum;

				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
				for (n=0;n<N;n++)
					sum += J[0][n]*dti*cellMask[n]*dw[i][0][n]*(w0[j][1][n]*w0[k][2][n]+
						0.5f*dw[j][1][n]*w0[k][2][n]+0.5f*w0[j][1][n]*dw[k][2][n]+0.3333333f*dw[j][1][n]*dw[k][2][n]);
				xtile[i+1][j+1][k+1][1] += sum;

				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
				for (n=0;n<N;n++)
					sum += J[0][n]*dti*cellMask[n]*dw[j][1][n]*(w0[i][0][n]*w0[k][2][n]+
						0.5f*dw[i][0][n]*w0[k][2][n]+0.5f*w0[i][0][n]*dw[k][2][n]+0.3333333f*dw[i][0][n]*dw[k][2][n]);
				xtile[i+1][j+1][k+1][2] += sum;

				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
				for (n=0;n<N;n++)
					sum += J[0][n]*dti*cellMask[n]*dw[k][2][n]*(w0[i][0][n]*w0[j][1][n]+
						0.5f*dw[i][0][n]*w0[j][1][n]+0.5f*w0[i][0][n]*dw[j][1][n]+0.3333333f*dw[i][0][n]*dw[j][1][n]);
				xtile[i+1][j+1][k+1][3] += sum;
			}

	// Integrate to get current from density decomposition
	// The relationship is J_i = j_i+1 + W_i
	for (j=0;j<5;j++)
		for (k=0;k<5;k++)
		{
			sum = xtile[4][j][k][1];
			sum += xtile[3][j][k][1]; xtile[3][j][k][1] = sum;
			sum += xtile[2][j][k][1]; xtile[2][j][k][1] = sum;
			sum += xtile[1][j][k][1]; xtile[1][j][k][1] = sum;
			xtile[0][j][k][1] = 0.0;

			sum = xtile[j][4][k][2];
			sum += xtile[j][3][k][2]; xtile[j][3][k][2] = sum;
			sum += xtile[j][2][k][2]; xtile[j][2][k][2] = sum;
			sum += xtile[j][1][k][2]; xtile[j][1][k][2] = sum;
			xtile[j][0][k][2] = 0.0;

			sum = xtile[j][k][4][3];
			sum += xtile[j][k][3][3]; xtile[j][k][3][3] = sum;
			sum += xtile[j][k][2][3]; xtile[j][k][2][3] = sum;
			sum += xtile[j][k][1][3]; xtile[j][k][1][3] = sum;
			xtile[j][k][0][3] = 0.0;
		}
}

void ParticleBundlePGC::Push(Species *owner)
{
	// The empty particles in the bundle can be operated on harmlessly
	// except in gather scatter.
	// N = particles in a full bundle
	// num = particles in this bundle
	Grid *theGrid = owner->owner;
	const tw::Float q0 = owner->charge;
	const tw::Float m0 = owner->restMass;
	const tw::Float dth = theGrid->dth;
	const tw::Float dt = theGrid->dt;
	const tw::Float k[3] = { dxi(*theGrid) , dyi(*theGrid) , dzi(*theGrid) };
	const float dti = 1.0/dt;

	PadBundle();

	cell0 = cell[0];
	theGrid->DecodeCell(cell0,&ijk0[0],&ijk0[1],&ijk0[2]);
	theGrid->GetWeights(w0,x);
	theGrid->GetWallWeights(l0,x);
	LoadTile();
	GatherF(F,w0,l0);
	GatherLaser(las,w0);

	avg_mass_1(avgMass,vel,p,F,las,q0,m0,dth);
	impulse(p,F,las,avgMass,q0,dth);
	rotation1(t,F,avgMass,q0,dth);
	rotation2(s,t);
	rotation3(s,t,vel,p);
	impulse(p,F,las,avgMass,q0,dth);
	avg_mass_2(avgMass,p,las,q0,m0,dth);
	velocity(vel,p,avgMass);
	translate(x,vel,k,dt);
	load_j4(J,number,vel,k,q0);
	load_chi(chi,number,avgMass,q0);

	theGrid->MinimizePrimitive(cell,ijk,x,domainMask);
	theGrid->GetWeights(w1,x);
	set_cell_mask(cellMask,cell0,cell);
	ResetXTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	ScatterChi(chi,w0,w1,cellMask);
	StoreXTile();
}

void ParticleBundlePGC::GatherLaser(float las[8][N],const float w[3][3][N])
{
	tw::Int i,j,k,n;
	alignas(AB) float factorNow[N];
	ZeroArray(las,0,7);
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			for (k=0;k<3;k++)
			{
				#pragma omp simd aligned(w,factorNow:AB)
				for (n=0;n<N;n++)
					factorNow[n] = w[i][0][n]*w[j][1][n]*w[k][2][n];
				#pragma omp simd aligned(las,factorNow:AB)
				for (n=0;n<N;n++)
				{
					las[0][n] += factorNow[n]*las_tile[i][j][k][0];
					las[1][n] += factorNow[n]*las_tile[i][j][k][1];
					las[2][n] += factorNow[n]*las_tile[i][j][k][2];
					las[3][n] += factorNow[n]*las_tile[i][j][k][3];
					las[4][n] += factorNow[n]*las_tile[i][j][k][4];
					las[5][n] += factorNow[n]*las_tile[i][j][k][5];
					las[6][n] += factorNow[n]*las_tile[i][j][k][6];
					las[7][n] += factorNow[n]*las_tile[i][j][k][7];
				}
			}
}

void ParticleBundlePGC::ScatterChi(const float chi[N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N])
{
	tw::Int i,j,k,dc[3],n;
	float xw0[5][3];
	float xw1[5][3];
	float sum;

	// Scalar code to scatter current from inter-cell particles.
	for (n=0;n<num;n++)
	{
		if (cellMask[n]==0.0)
		{
			get_cell_displ(dc,n);

			for (i=0;i<5;i++)
				for (j=0;j<3;j++)
				{
					xw0[i][j] = 0.0f;
					xw1[i][j] = 0.0f;
				}
			for (i=0;i<3;i++)
				for (j=0;j<3;j++)
				{
					xw0[i+1][j] = w0[i][j][n];
					xw1[i+1+dc[j]][j] = w1[i][j][n];
				}

			for (i=0;i<5;i++)
				for (j=0;j<5;j++)
					for (k=0;k<5;k++)
						chi_tile[i][j][k] += 0.5f*chi[n]*(xw0[i][0]*xw0[j][1]*xw0[k][2] + xw1[i][0]*xw1[j][1]*xw1[k][2]);
		}
	}

	// Vector code to scatter current from intra-cell particles
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			for (k=0;k<3;k++)
			{
				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(chi,w0,w1,cellMask:AB)
				for (n=0;n<N;n++)
					sum += 0.5f*chi[n]*cellMask[n]*(w0[i][0][n]*w0[j][1][n]*w0[k][2][n] + w1[i][0][n]*w1[j][1][n]*w1[k][2][n]);
				chi_tile[i+1][j+1][k+1] += sum;
			}
}
