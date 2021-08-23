#include "simulation.h"
#include "particles.h"
#include "pusher.h"
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
		// if (low[i]<owner->N0(i))
		// 	throw tw::FatalError("Particle unexpectedly low");
		// if (high[i]>owner->N1(i))
		// 	throw tw::FatalError("Particle unexpectedly high");
	}
}

void Species::BunchTasks(std::vector<tw::Int>& task_map)
{
	for (tw::Int i=0;i<task_map.size();i++)
		task_map[i] = i;
}

void Species::SpreadTasks(std::vector<tw::Int>& task_map)
{
	// Try to order tasks so they are spread out during concurrent execution.
	// Gives us a chance of non-overlapping memory in the main source field.
	const tw::Int num_tasks = task_map.size();
	std::vector<tw::Int> assigned_tasks;
	std::vector<tw::Int> unassigned_tasks(num_tasks);
	for (tw::Int i=0;i<num_tasks;i++)
		unassigned_tasks[i] = i;
	// Assing the first task to the first slot
	task_map[0] = *unassigned_tasks.begin();
	assigned_tasks.push_back(*unassigned_tasks.begin());
	unassigned_tasks.erase(unassigned_tasks.begin());
	// Assign remaining tasks maximizing the distance to previously assigned tasks
	for (tw::Int slot=1;slot<num_tasks;slot++)
	{
		tw::Int ans;
		tw::Float distance = 0.0;
		for (auto utask : unassigned_tasks)
		{
			tw::Float test = 0.0;
			for (auto atask : assigned_tasks)
				test += 1.0/sqr(tw::Float(atask-utask));
			test = 1.0/test;
			if (test>distance)
			{
				ans = utask;
				distance = test;
			}
		}
		task_map[slot] = ans;
		assigned_tasks.push_back(ans);
		unassigned_tasks.erase(std::find(unassigned_tasks.begin(),unassigned_tasks.end(),ans));
	}
}

template <class BundleType>
void Species::PushSlice(tw::Int tasks,tw::Int tid,tw::Int bounds_data[][8])
{
	tw::Int first,last,next,low[4],high[4];
	std::vector<ParticleRef> map;
	BundleType b;
	first = bounds_data[tid][0];
	last = bounds_data[tid][1];
	map.resize(last-first+1);
	for (tw::Int i=first;i<=last;i++)
		map[i-first] = ParticleRef(i,particle[i]);
	std::sort(map.begin(),map.end());
	GetSubarrayBounds(map,low,high,1);
	b.LoadFieldSlice(this,low,high,ignorable);
	GetSubarrayBounds(map,low,high,2);
	b.InitSourceSlice(this,low,high,ignorable);
	// Save the bounds information
	for (tw::Int i=1;i<=3;i++)
	{
		bounds_data[tid][i*2] = low[i];
		bounds_data[tid][i*2+1] = high[i];
	}
	for (tw::Int i=first;i<=last;i++)
	{
		b.Append(particle[map[i-first].idx],restMass);
		next = i==last ? i : i+1;
		if (i==last || b.Complete(particle[map[next-first].idx]))
		{
			b.Advance(this);
			b.CopyBack(this);
			b.Reset();
		}
	}
	// Work out whether atomic operations are needed or not
	bool needs_atomic = false;
	#pragma omp barrier
	for (tw::Int i=0;i<tasks;i++)
	{
		if (i!=tid)
		{
			bool disjoint =
				bounds_data[tid][2]>bounds_data[i][3] ||
				bounds_data[tid][4]>bounds_data[i][5] ||
				bounds_data[tid][6]>bounds_data[i][7] ||
				bounds_data[tid][3]<bounds_data[i][2] ||
				bounds_data[tid][5]<bounds_data[i][4] ||
				bounds_data[tid][7]<bounds_data[i][6];
			needs_atomic = needs_atomic || !disjoint;
		}
	}
	needs_atomic = needs_atomic && tw::GetOMPMaxThreads()>1;
	b.DepositSourceSlice(this,needs_atomic);
}

template <class BundleType>
void Species::Push()
{
	if (particle.size()==0)
		return;

	const tw::Int min_particles_per_task = 256;
	const tw::Int num_par = particle.size();
	const tw::Int concurrent_tasks = tw::GetOMPMaxThreads();
	const tw::Int max_tasks = 1 + num_par / min_particles_per_task;
	const tw::Int preferred_tasks = 32*concurrent_tasks;
	const tw::Int num_tasks = preferred_tasks > max_tasks ? max_tasks : preferred_tasks;
	const tw::Int concurrent_sets = num_tasks / concurrent_tasks;
	const tw::Int remainder_tasks = num_tasks % concurrent_tasks;
	const tw::Int total_sets = concurrent_sets + (remainder_tasks?1:0);

	std::vector<tw::Int> task_map(num_tasks);
	SpreadTasks(task_map);
	//BunchTasks(task_map);

	// following has first,last,x0,x1,y0,y1,z0,z1
	tw::Int bounds_data[concurrent_tasks][8];
	for (tw::Int c=0;c<total_sets;c++)
	{
		tw::Int tasks_in_set = c<concurrent_sets ? concurrent_tasks : remainder_tasks;
		#pragma omp parallel num_threads(tasks_in_set)
		{
			tw::Int t = tw::GetOMPThreadNum();
			tw::Int task_idx = task_map[c*concurrent_tasks + t];
			tw::GetOMPTaskLoopRange(task_idx,num_par,num_tasks,&bounds_data[t][0],&bounds_data[t][1]);
			PushSlice<BundleType>(tasks_in_set,t,bounds_data);
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
	if (ESField)
	{
		//Push<ParticleBundleElectrostatic>();
		return;
	}
	if (laser)
	{
		Push<ParticleBundlePGC>();
		return;
	}

	if (ignorable[2])
		Push<ParticleBundle2D>();
		//Push<ParticleBundleUnitary>();
	else
		Push<ParticleBundle3D>();
}


///////////////////////////////
//                           //
//     BUNDLE OPERATIONS     //
//                           //
///////////////////////////////

void ParticleBundle2D::Advance(Species *owner)
{
	Gather(owner,F);
	Push(owner);
	Scatter(owner);
}

void ParticleBundle3D::Advance(Species *owner)
{
	Gather(owner,F);
	Push(owner);
	Scatter(owner);
}

void ParticleBundlePGC::Advance(Species *owner)
{
	Gather(owner,F,las);
	Push(owner);
	Scatter(owner,chi);
}

void ParticleBundleUnitary::Advance(Species *owner)
{
	Gather(owner,F);
	Push(owner);
	Scatter(owner);
}

void ParticleBundleBohmian::Advance(Species *owner)
{
	Gather(owner);
	Push(owner);
}

void BundlePusherBoris::Push(Species *owner)
{
	Simulation *sim = owner->owner;
	const tw::Float q0 = owner->charge;
	const tw::Float dt = timestep(*sim);
	const tw::Float k[3] = { dxi(*sim) , dyi(*sim) , dzi(*sim) };

	impulse(u,F);
	rotation1(t,u,F);
	rotation2(s,t);
	rotation3(s,t,vel,u);
	impulse(u,F);
	velocity(vel,u);
	translate(x,vel,k,dt);
	load_j4(J,number,vel,k,q0);
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

void BundlePusherUnitary::Push(Species *owner)
{
	// We rewrite L*z*L^dag as (L*(L*z)^dag)^dag
	Simulation *sim = owner->owner;
	const tw::Float q0 = owner->charge;
	const tw::Float dt = timestep(*sim);
	const tw::Float k[3] = { dxi(*sim) , dyi(*sim) , dzi(*sim) };

	to_spinor(u,zi);
	estimate_ds(ds,F,u,dt); // ds gets buried in F; at present must go after to_spinor
	Lambda();
	copy_spinor(zi,zf);
	dagger(zi);
	Lambda();
	dagger(zf);
	to_vector(u,zf,a);

	velocity(vel,u);
	translate(x,vel,k,dt);
	load_j4(J,number,vel,k,q0);
}

void BundlePusherPGC::Push(Species *owner)
{
	Simulation *sim = owner->owner;
	const tw::Float q0 = owner->charge;
	const tw::Float m0 = owner->restMass;
	const tw::Float dth = 0.5*timestep(*sim);
	const tw::Float dt = timestep(*sim);
	const tw::Float k[3] = { dxi(*sim) , dyi(*sim) , dzi(*sim) };

	avg_gam_1(avgGam,vel,u,F,las,dth);
	impulse(u,F,las,avgGam,dth);
	rotation1(t,F,avgGam);
	rotation2(s,t);
	rotation3(s,t,vel,u);
	impulse(u,F,las,avgGam,dth);
	avg_gam_2(avgGam,u,las,dth);
	velocity(vel,u,avgGam);
	translate(x,vel,k,dt);
	load_j4(J,number,vel,k,q0);
	load_chi(chi,number,avgGam,q0,m0);
}

void BundlePusherBohmian::Push(Species *owner)
{
	Simulation *sim = owner->owner;
	const tw::Float dt = timestep(*sim);
	const tw::Float k[3] = { dxi(*sim) , dyi(*sim) , dzi(*sim) };
	// estimate vel(n+1/2) using u(n-1) and J(n), and update u(n-1) to u(n)
	bohm_velocity(vel,u,J);
	// update x(n) to x(n+1) using vel(n+1/2)
	translate(x,vel,k,dt);
	sim->MinimizePrimitive(cell,ijk,x,domainMask);
}

void BundleTiler2D::Gather(Species *owner,float F[6][N])
{
	Simulation *sim = owner->owner;
	const float qmdth = 0.5 * owner->charge * timestep(*sim) / owner->restMass;
	PadBundle();
	cell0 = cell[0];
	sim->DecodeCell(cell0,&ijk0[0],&ijk0[1],&ijk0[2]);
	sim->GetWeights(w0,x);
	sim->GetWallWeights(l0,x);
	LoadFTile();
	GatherF(F,w0,l0,qmdth);
}

void BundleTiler2D::Scatter(Species *owner)
{
	Simulation *sim = owner->owner;
	const tw::Float dti = 1.0/timestep(*sim);
	sim->MinimizePrimitive(cell,ijk,x,domainMask);
	sim->GetWeights(w1,x);
	set_cell_mask(cellMask,cell0,cell);
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile();
}

void BundleTiler2D::GatherF(float F[6][N],const float w[3][3][N],const float l[3][3][N],const float qmdth)
{
	// Assumes every particle in bundle is in the same cell and data is packed with Yee fields
	tw::Int i,k,n;
	ZeroArray(F,0,5);
	for (i=0;i<3;i++)
		for (k=0;k<3;k++)
			#pragma omp simd aligned(F,w,l:AB)
			for (n=0;n<N;n++)
			{
				F[0][n] += l[i][0][n]*w[k][2][n]*F_tile[i][k][0]*qmdth;
				F[1][n] += w[i][0][n]*w[k][2][n]*F_tile[i][k][1]*qmdth;
				F[2][n] += w[i][0][n]*l[k][2][n]*F_tile[i][k][2]*qmdth;
				F[3][n] += w[i][0][n]*l[k][2][n]*F_tile[i][k][3]*qmdth;
				F[4][n] += l[i][0][n]*l[k][2][n]*F_tile[i][k][4]*qmdth;
				F[5][n] += l[i][0][n]*w[k][2][n]*F_tile[i][k][5]*qmdth;
			}
}

void BundleTiler2D::ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti)
{
	// J_tile will be loaded with charge and current
	// if J_tile is nonzero on entry, it must contain the density decomposition, not the current
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
		if (cellMask[n]!=0.0)
			continue;

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
				J_tile[i][k][0] += J[0][n]*xw1[i][0]*xw1[k][2];

		// Optimize Current Deposition Loops
		i1 = dc[0]<0 ? 0 : 1;
		i2 = dc[0]>0 ? 4 : 3;
		k1 = dc[2]<0 ? 0 : 1;
		k2 = dc[2]>0 ? 4 : 3;

		// Current Deposition (as in amps through the wall)
		for (i=1;i<=i2;i++)
			for (k=k1;k<=k2;k++)
				J_tile[i][k][1] += J[0][n]*dti*xdw[i][0]*(xw0[k][2]+0.5f*xdw[k][2]);
		for (i=i1;i<=i2;i++)
			for (k=k1;k<=k2;k++)
				J_tile[i][k][2] += J[2][n]*(xw0[i][0]*xw0[k][2]+0.5f*xdw[i][0]*xw0[k][2]+0.5f*xw0[i][0]*xdw[k][2]+0.3333333f*xdw[i][0]*xdw[k][2]);
		for (i=i1;i<=i2;i++)
			for (k=1;k<=k2;k++)
				J_tile[i][k][3] += J[0][n]*dti*xdw[k][2]*(xw0[i][0]+0.5f*xdw[i][0]);
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
			J_tile[i+1][k+1][0] += sum;

			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
			for (n=0;n<N;n++)
				sum += J[0][n]*dti*cellMask[n]*dw[i][0][n]*(w0[k][2][n]+0.5f*dw[k][2][n]);
			J_tile[i+1][k+1][1] += sum;

			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
			for (n=0;n<N;n++)
				sum += J[2][n]*cellMask[n]*(w0[i][0][n]*w0[k][2][n]+0.5f*dw[i][0][n]*w0[k][2][n]+0.5f*w0[i][0][n]*dw[k][2][n]+0.3333333f*dw[i][0][n]*dw[k][2][n]);
			J_tile[i+1][k+1][2] += sum;

			sum = 0.0;
			#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
			for (n=0;n<N;n++)
				sum += J[0][n]*dti*cellMask[n]*dw[k][2][n]*(w0[i][0][n]+0.5f*dw[i][0][n]);
			J_tile[i+1][k+1][3] += sum;
		}

	// Integrate to get current from density decomposition
	// The relationship is J_i = j_i+1 + W_i
	for (k=0;k<5;k++)
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

void BundleTiler3D::Gather(Species *owner,float F[6][N])
{
	Simulation *sim = owner->owner;
	const float qmdth = 0.5 * owner->charge * timestep(*sim) / owner->restMass;
	PadBundle();
	cell0 = cell[0];
	sim->DecodeCell(cell0,&ijk0[0],&ijk0[1],&ijk0[2]);
	sim->GetWeights(w0,x);
	sim->GetWallWeights(l0,x);
	LoadFTile();
	GatherF(F,w0,l0,qmdth);
}

void BundleTiler3D::Scatter(Species *owner)
{
	Simulation *sim = owner->owner;
	const tw::Float dti = 1.0/timestep(*sim);
	sim->MinimizePrimitive(cell,ijk,x,domainMask);
	sim->GetWeights(w1,x);
	set_cell_mask(cellMask,cell0,cell);
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile();
}

void BundleTiler3D::GatherF(float F[6][N],const float w[3][3][N],const float l[3][3][N],const float qmdth)
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
					F[0][n] += l[i][0][n]*w[j][1][n]*w[k][2][n]*F_tile[i][j][k][0]*qmdth;
					F[1][n] += w[i][0][n]*l[j][1][n]*w[k][2][n]*F_tile[i][j][k][1]*qmdth;
					F[2][n] += w[i][0][n]*w[j][1][n]*l[k][2][n]*F_tile[i][j][k][2]*qmdth;
					F[3][n] += w[i][0][n]*l[j][1][n]*l[k][2][n]*F_tile[i][j][k][3]*qmdth;
					F[4][n] += l[i][0][n]*w[j][1][n]*l[k][2][n]*F_tile[i][j][k][4]*qmdth;
					F[5][n] += l[i][0][n]*l[j][1][n]*w[k][2][n]*F_tile[i][j][k][5]*qmdth;
				}
}

void BundleTiler3D::ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti)
{
	// J_tile will be loaded with charge and current
	// if J_tile is nonzero on entry, it must contain the density decomposition, not the current
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
		if (cellMask[n]!=0.0)
			continue;

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
					J_tile[i][j][k][0] += J[0][n]*xw1[i][0]*xw1[j][1]*xw1[k][2];

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
					J_tile[i][j][k][1] += J[0][n]*dti*xdw[i][0]*(xw0[j][1]*xw0[k][2]+
						0.5f*xdw[j][1]*xw0[k][2]+0.5f*xw0[j][1]*xdw[k][2]+0.3333333f*xdw[j][1]*xdw[k][2]);
		for (i=i1;i<=i2;i++)
			for (j=1;j<=j2;j++)
				for (k=k1;k<=k2;k++)
					J_tile[i][j][k][2] += J[0][n]*dti*xdw[j][1]*(xw0[i][0]*xw0[k][2]+
						0.5f*xdw[i][0]*xw0[k][2]+0.5f*xw0[i][0]*xdw[k][2]+0.3333333f*xdw[i][0]*xdw[k][2]);
		for (i=i1;i<=i2;i++)
			for (j=j1;j<=j2;j++)
				for (k=1;k<=k2;k++)
					J_tile[i][j][k][3] += J[0][n]*dti*xdw[k][2]*(xw0[i][0]*xw0[j][1]+
						0.5f*xdw[i][0]*xw0[j][1]+0.5f*xw0[i][0]*xdw[j][1]+0.3333333f*xdw[i][0]*xdw[j][1]);
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
				J_tile[i+1][j+1][k+1][0] += sum;

				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
				for (n=0;n<N;n++)
					sum += J[0][n]*dti*cellMask[n]*dw[i][0][n]*(w0[j][1][n]*w0[k][2][n]+
						0.5f*dw[j][1][n]*w0[k][2][n]+0.5f*w0[j][1][n]*dw[k][2][n]+0.3333333f*dw[j][1][n]*dw[k][2][n]);
				J_tile[i+1][j+1][k+1][1] += sum;

				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
				for (n=0;n<N;n++)
					sum += J[0][n]*dti*cellMask[n]*dw[j][1][n]*(w0[i][0][n]*w0[k][2][n]+
						0.5f*dw[i][0][n]*w0[k][2][n]+0.5f*w0[i][0][n]*dw[k][2][n]+0.3333333f*dw[i][0][n]*dw[k][2][n]);
				J_tile[i+1][j+1][k+1][2] += sum;

				sum = 0.0;
				#pragma omp simd reduction(+:sum) aligned(J,dw,w0,cellMask:AB)
				for (n=0;n<N;n++)
					sum += J[0][n]*dti*cellMask[n]*dw[k][2][n]*(w0[i][0][n]*w0[j][1][n]+
						0.5f*dw[i][0][n]*w0[j][1][n]+0.5f*w0[i][0][n]*dw[j][1][n]+0.3333333f*dw[i][0][n]*dw[j][1][n]);
				J_tile[i+1][j+1][k+1][3] += sum;
			}

	// Integrate to get current from density decomposition
	// The relationship is J_i = j_i+1 + W_i
	for (j=0;j<5;j++)
		for (k=0;k<5;k++)
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

void BundleTilerPGC::Gather(Species *owner,float F[6][N],float las[8][N])
{
	const float q2m2 = sqr(owner->charge / owner->restMass);
	BundleTiler3D::Gather(owner,F);
	LoadLaserTile();
	GatherLaser(las,w0,q2m2);
}

void BundleTilerPGC::Scatter(Species *owner,float chi[N])
{
	BundleTiler3D::Scatter(owner);
	ResetChiTile();
	ScatterChi(chi,w0,w1,cellMask);
	StoreChiTile();
}

void BundleTilerPGC::GatherLaser(float las[8][N],const float w[3][3][N],const float q2m2)
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
					las[0][n] += factorNow[n]*las_tile[i][j][k][0]*q2m2;
					las[1][n] += factorNow[n]*las_tile[i][j][k][1]*q2m2;
					las[2][n] += factorNow[n]*las_tile[i][j][k][2]*q2m2;
					las[3][n] += factorNow[n]*las_tile[i][j][k][3]*q2m2;
					las[4][n] += factorNow[n]*las_tile[i][j][k][4]*q2m2;
					las[5][n] += factorNow[n]*las_tile[i][j][k][5]*q2m2;
					las[6][n] += factorNow[n]*las_tile[i][j][k][6]*q2m2;
					las[7][n] += factorNow[n]*las_tile[i][j][k][7]*q2m2;
				}
			}
}

void BundleTilerPGC::ScatterChi(const float chi[N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N])
{
	tw::Int i,j,k,dc[3],n;
	float xw0[5][3];
	float xw1[5][3];
	float sum;

	// Scalar code to scatter current from inter-cell particles.
	for (n=0;n<num;n++)
	{
		if (cellMask[n]!=0.0)
			continue;

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

void BundleTilerBohmian::Gather(Species *owner)
{
	Simulation *sim = owner->owner;
	PadBundle();
	cell0 = cell[0];
	sim->DecodeCell(cell0,&ijk0[0],&ijk0[1],&ijk0[2]);
	sim->GetWeights(w0,x);
	LoadTile();
	GatherJ4(J,w0);
}

void BundleTilerBohmian::GatherJ4(float J[4][N],const float w[3][3][N])
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
