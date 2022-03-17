#include "meta_base.h"
#include "computeTool.h"
#include "bundle.h"
#include "pusher.h"
#include "slicer.h"
#include "tiler.h"
#include "mover.h"

Mover::Mover(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	for (tw::Int i=0;i<4;i++)
		ignorable[i] = space->Ignorable(i);
}

void Mover::AddTransferParticle(const Particle& src)
{
	tw::Int ijk[4];
	TransferParticle dest;
	space->DecodeCell(src.q,ijk);
	dest.dst[0] = task->strip[0].Get_rank();
	for (tw::Int i=1;i<=3;i++)
		dest.dst[i] = tw::Int(ijk[i]>space->Dim(i)) - tw::Int(ijk[i]<1);
	dest.x = tw::vec4(0.0,space->PositionFromPrimitive(src.q));
	dest.p = tw::vec4(0.0,src.p);
	dest.number = src.number;
	dest.aux1 = src.aux1;
	dest.aux2 = src.aux2;
	transfer->push_back(dest);
}

void Mover::GetSubarrayBounds(std::vector<ParticleRef>& sorted,tw::Int low[4],tw::Int high[4],tw::Int layers)
{
	// Assumes particles are sorted in increasing memory order
	// This depends in turn on the cell encoding respecting memory order

	space->DecodeCell(sorted.front().cell,low);
	space->DecodeCell(sorted.back().cell,high);
	// Sorting by cell only sorts the outermost topological index in general.
	// Hence the following loop.
	for (int i=1;i<=3;i++)
		if (low[i]>high[i])
			std::swap(low[i],high[i]);
	// Assume z-packing and maximal loading
	if (low[1]!=high[1])
	{
		low[2] = 1;
		high[2] = space->Dim(2);
		low[3] = 1;
		high[3] = space->Dim(3);
	}
	if (low[2]!=high[2])
	{
		low[3] = 1;
		high[3] = space->Dim(3);
	}
	// Expand subarray to allow for motion and particle shape
	// For gather this induces 1 ghost cell layers
	// For scatter this induces 2 ghost cell layers
	for (int i=1;i<=3;i++)
	{
		if (space->Dim(i)>1)
		{
			low[i] -= layers;
			high[i] += layers;
		}
		// ASSERT_GTREQ(low[i],space->LFG(i));
		// ASSERT_LESSEQ(high[i],space->UFG(i));
	}
}

void Mover::BunchTasks(std::vector<tw::Int>& task_map)
{
	for (tw::Int i=0;i<task_map.size();i++)
		task_map[i] = i;
}

void Mover::SpreadTasks(std::vector<tw::Int>& task_map)
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
void Mover::MoveSlice(tw::Int tasks,tw::Int tid,tw::Int bounds_data[][8])
{
	tw::Int first,last,next,low[4],high[4];
	std::vector<ParticleRef> map;
	BundleType b(this);
	first = bounds_data[tid][0];
	last = bounds_data[tid][1];
	map.reserve(last-first+1);
	for (tw::Int i=first;i<=last;i++)
	{
		assert(space->IsRefCellWithin((*particle)[i].q,2));
		map.push_back(ParticleRef(i,(*particle)[i]));
	}
	std::sort(map.begin(),map.end());
	GetSubarrayBounds(map,low,high,1);
	b.LoadFieldSlice(low,high,ignorable);
	GetSubarrayBounds(map,low,high,2);
	b.InitSourceSlice(low,high,ignorable);
	// Save the bounds information
	for (tw::Int i=1;i<=3;i++)
	{
		bounds_data[tid][i*2] = low[i];
		bounds_data[tid][i*2+1] = high[i];
	}
	for (tw::Int i=first;i<=last;i++)
	{
		b.Append((*particle)[map[i-first].idx]);
		next = i==last ? i : i+1;
		if (i==last || b.Complete((*particle)[map[next-first].idx]))
		{
			b.Move();
			b.CopyBack();
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
	b.DepositSourceSlice(needs_atomic);
}

template <class BundleType>
void Mover::DoTasks()
{
	if (particle->size()==0)
		return;

	const tw::Int min_particles_per_task = 256;
	const tw::Int num_par = particle->size();
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
			MoveSlice<BundleType>(tasks_in_set,t,bounds_data);
		}
	}
}

void Mover::Advance() {}

void BorisMover::Advance()
{
	if (space->Dim(2)==1)
		DoTasks<BundleMoverBoris2D>();
	else
		DoTasks<BundleMoverBoris3D>();
}

void UnitaryMover::Advance()
{
	if (space->Dim(2)==1)
		DoTasks<BundleMoverUnitary2D>();
	else
		DoTasks<BundleMoverUnitary3D>();
}

void PGCMover::Advance()
{
	if (space->Dim(2)==1)
		DoTasks<BundleMoverPGC2D>();
	else
		DoTasks<BundleMoverPGC3D>();
}

void BohmianMover::Advance()
{
	if (space->Dim(2)==1)
		DoTasks<BundleMoverBohmian2D>();
	else
		DoTasks<BundleMoverBohmian3D>();
}

///////////////////////////
//                       //
//     BUNDLE MOVERS     //
//                       //
///////////////////////////

void BundleMoverBoris2D::Move()
{
	const float qmdth = 0.5*q0*dt/m0;
	const tw::Float dti = 1.0/dt;
	PrepareGather();
	LoadFTile();
	GatherF(F,w0,l0,qmdth);
	Push();
	PrepareScatter();
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile();
}

void BundleMoverBoris3D::Move()
{
	const float qmdth = 0.5*q0*dt/m0;
	const tw::Float dti = 1.0/dt;
	PrepareGather();
	LoadFTile();
	GatherF(F,w0,l0,qmdth);
	Push();
	PrepareScatter();
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile();
}

void BundleMoverUnitary2D::Move()
{
	const float qmdth = 0.5*q0*dt/m0;
	const tw::Float dti = 1.0/dt;
	PrepareGather();
	LoadFTile();
	GatherF(F,w0,l0,qmdth);
	Push();
	PrepareScatter();
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile();
}

void BundleMoverUnitary3D::Move()
{
	const float qmdth = 0.5*q0*dt/m0;
	const tw::Float dti = 1.0/dt;
	PrepareGather();
	LoadFTile();
	GatherF(F,w0,l0,qmdth);
	Push();
	PrepareScatter();
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile();
}

void BundleMoverPGC2D::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	BundleSlicerEM::LoadFieldSlice(low,high,ignorable);
	BundleSlicerPGC::LoadFieldSlice(low,high,ignorable);
}

void BundleMoverPGC2D::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	BundleSlicerEM::InitSourceSlice(low,high,ignorable);
	BundleSlicerPGC::InitSourceSlice(low,high,ignorable);
}

void BundleMoverPGC2D::DepositSourceSlice(bool needsAtomic)
{
	BundleSlicerEM::DepositSourceSlice(needsAtomic);
	BundleSlicerPGC::DepositSourceSlice(needsAtomic);
}

void BundleMoverPGC2D::Move()
{
	const float qmdth = 0.5*q0*dt/m0;
	const tw::Float dti = 1.0/dt;
	const float q2m2 = sqr(q0/m0);
	PrepareGather();
	LoadFTile();
	GatherF(F,w0,l0,qmdth);
	LoadLaserTile();
	GatherLaser(las,w0,q2m2);
	Push();
	PrepareScatter();
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile();
	ResetChiTile();
	ScatterChi(chi,w0,w1,cellMask);
	StoreChiTile();
}

void BundleMoverPGC3D::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	BundleSlicerEM::LoadFieldSlice(low,high,ignorable);
	BundleSlicerPGC::LoadFieldSlice(low,high,ignorable);
}

void BundleMoverPGC3D::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	BundleSlicerEM::InitSourceSlice(low,high,ignorable);
	BundleSlicerPGC::InitSourceSlice(low,high,ignorable);
}

void BundleMoverPGC3D::DepositSourceSlice(bool needsAtomic)
{
	BundleSlicerEM::DepositSourceSlice(needsAtomic);
	BundleSlicerPGC::DepositSourceSlice(needsAtomic);
}

void BundleMoverPGC3D::Move()
{
	const float qmdth = 0.5*q0*dt/m0;
	const tw::Float dti = 1.0/dt;
	const float q2m2 = sqr(q0/m0);
	PrepareGather();
	LoadFTile();
	GatherF(F,w0,l0,qmdth);
	LoadLaserTile();
	GatherLaser(las,w0,q2m2);
	Push();
	PrepareScatter();
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile();
	ResetChiTile();
	ScatterChi(chi,w0,w1,cellMask);
	StoreChiTile();
}

void BundleMoverBohmian2D::Move()
{
	PadBundle();
	cell0 = cell[0];
	owner->space->DecodeCell(cell0,ijk0);
	owner->space->GetWeights(w0,x);
	LoadTile();
	GatherJ4(J,w0);
	Push();
	owner->space->MinimizePrimitive(cell,ijk,x,domainMask);
}

void BundleMoverBohmian3D::Move()
{
	PadBundle();
	cell0 = cell[0];
	owner->space->DecodeCell(cell0,ijk0);
	owner->space->GetWeights(w0,x);
	LoadTile();
	GatherJ4(J,w0);
	Push();
	owner->space->MinimizePrimitive(cell,ijk,x,domainMask);
}
