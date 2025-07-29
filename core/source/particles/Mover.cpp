module;

#include "tw_includes.h"
#include <algorithm>

export module mover;
import compute_tool;
import fields;
import bundle;
import pusher;
import tiler;

/// @brief Orchestrates moving particle bundles.  The `Advance` function associated with
/// a specific subclass will invoke `DoTasks` with a template argument that selects
/// a particular bundle mover.
export struct Mover:ComputeTool
{
	tw::Float q0,m0;
	tw::Int ignorable[4];
	// Following particles and fields are owned by the parent Species module
	std::vector<Particle> *particle;
	std::vector<TransferParticle> *transfer;
	Field *ESField,*EM,*sources,*laser,*chi,*qo_j4;
	Mover(const std::string& name,MetricSpace *m,Task *tsk);
	MoverParams GetParams();
	void AddTransferParticle(const Particle& src);
	template <class MoverType>
	void CopyBack(MoverType *b);
	std::pair<tw::node5,tw::node5> GetSubarrayBounds(std::vector<ParticleRef>& sorted,tw::Int layers);
	void SpreadTasks(std::vector<tw::Int>& task_map);
	void BunchTasks(std::vector<tw::Int>& task_map);
	template <class MoverType>
	void MoveSlice(tw::Int tasks,tw::Int tid,tw::Int bounds_data[][8]);
	template <class MoverType>
	void DoTasks();
	virtual void Advance() {}

	virtual void InitTest();
	virtual void EncodingTest();
	virtual void MinimizePrimitiveScalarTest();
	virtual void MinimizePrimitiveVectorTest();
	virtual void TranslationTest();
	virtual void UniformETest();
	virtual void UniformBTest();
	virtual void PlaneWaveTest();
	virtual bool Test(tw::Int& id);
	virtual void CloseTest();
};

// Subclasses of the Mover Tool
// These are simple dispatchers. Their only purpose is to make the
// appropriate templated call to DoTasks().  The template argument is one
// of the BundleMover* classes.

export struct BorisMover:Mover
{
	BorisMover(const std::string& name,MetricSpace *m,Task *tsk): Mover(name,m,tsk) {}
	virtual void Advance();
	virtual void InitTest();
};

export struct HCMover:Mover
{
	HCMover(const std::string& name,MetricSpace *m,Task *tsk): Mover(name,m,tsk) {}
	virtual void Advance();
	virtual void InitTest();
};

export struct UnitaryMover:Mover
{
	UnitaryMover(const std::string& name,MetricSpace *m,Task *tsk): Mover(name,m,tsk) {}
	virtual void Advance();
	virtual void InitTest();
};

export struct PGCMover:Mover
{
	PGCMover(const std::string& name,MetricSpace *m,Task *tsk): Mover(name,m,tsk) {}
	virtual void Advance();
	virtual void InitTest();
};

export struct BohmianMover:Mover
{
	BohmianMover(const std::string& name,MetricSpace *m,Task *tsk): Mover(name,m,tsk) {}
	virtual void Advance();
	virtual bool Test(tw::Int& id) { return false; }
};

export struct PhotonMover:Mover
{
	PhotonMover(const std::string& name,MetricSpace *m,Task *tsk): Mover(name,m,tsk) {}
	virtual void Advance();
	virtual void InitTest();
};

/////////////////////////////////////
// Bundle Movers - pusher + tiler  //
// These are created by the mover  //
// tool, and invoked by template   //
/////////////////////////////////////

struct BundleMoverBoris2D : BundleTilerEM2D,BundlePusherBoris
{
	Slice<float> Fx,Jx;
	BundleMoverBoris2D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerEM2D(mov), BundlePusherBoris(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverBoris3D : BundleTilerEM3D,BundlePusherBoris
{
	Slice<float> Fx,Jx;
	BundleMoverBoris3D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerEM3D(mov), BundlePusherBoris(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverHC2D : BundleTilerEM2D,BundlePusherHC
{
	Slice<float> Fx,Jx;
	BundleMoverHC2D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerEM2D(mov), BundlePusherHC(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverHC3D : BundleTilerEM3D,BundlePusherHC
{
	Slice<float> Fx,Jx;
	BundleMoverHC3D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerEM3D(mov), BundlePusherHC(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverUnitary2D : BundleTilerEM2D,BundlePusherUnitary
{
	Slice<float> Fx,Jx;
	BundleMoverUnitary2D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerEM2D(mov), BundlePusherUnitary(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverUnitary3D : BundleTilerEM3D,BundlePusherUnitary
{
	Slice<float> Fx,Jx;
	BundleMoverUnitary3D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerEM3D(mov), BundlePusherUnitary(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverPGC2D : BundleTilerPGC2D,BundleTilerEM2D,BundlePusherPGC
{
	Slice<float> Fx,Jx,lasx,chix;
	BundleMoverPGC2D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerPGC2D(mov), BundleTilerEM2D(mov), BundlePusherPGC(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverPGC3D : BundleTilerPGC3D,BundleTilerEM3D,BundlePusherPGC
{
	Slice<float> Fx,Jx,lasx,chix;
	BundleMoverPGC3D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerPGC3D(mov), BundleTilerEM3D(mov), BundlePusherPGC(mov) {}
	// For PGC we have to call two slicers explicitly
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverBohmian2D : BundleTilerBohmian2D,BundlePusherBohmian
{
	Slice<float> Jx;
	BundleMoverBohmian2D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerBohmian2D(mov), BundlePusherBohmian(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverBohmian3D : BundleTilerBohmian3D,BundlePusherBohmian
{
	Slice<float> Jx;
	BundleMoverBohmian3D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerBohmian3D(mov), BundlePusherBohmian(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverPhoton2D : BundleTilerEM2D,BundlePusherPhoton
{
	Slice<float> Fx,Jx;
	BundleMoverPhoton2D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerEM2D(mov), BundlePusherPhoton(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

struct BundleMoverPhoton3D : BundleTilerEM3D,BundlePusherPhoton
{
	Slice<float> Fx,Jx;
	BundleMoverPhoton3D(const MoverParams& mov) : ParticleBundle(mov), BundleTilerEM3D(mov), BundlePusherPhoton(mov) {}
	void LoadFieldSlice(tw::node5& beg,tw::node5& end);
	void InitSourceSlice(tw::node5& beg,tw::node5& end);
	void DepositSourceSlice(bool needsAtomic);
	void Move(tw::Float dts);
};

Mover::Mover(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	for (auto i=0;i<4;i++)
		ignorable[i] = m->Ignorable(i);
}

MoverParams Mover::GetParams() {
	MoverParams ans;
	ans.chi = chi;
	ans.EM = EM;
	ans.ESField = ESField;
	for (auto i=0;i<4;i++)
		ans.ignorable[i] = ignorable[i];
	ans.laser = laser;
	ans.m0 = m0;
	ans.q0 = q0;
	ans.ms = space;
	ans.qo_j4 = qo_j4;
	ans.sources = sources;
	return ans;
}

/// @brief Add particle to transfer list, assuming it is in extended domain
/// @param src the particle to add
void Mover::AddTransferParticle(const Particle& src)
{
	// we could have x out of normal range if there is an algorithm that moves particles
	// far enough to escape the extended domain
	TransferParticle dest;
	space->DecodeCell(src.q,dest.ijk);
	dest.dst[0] = task->strip[0].Get_rank();
	for (tw::Int i=1;i<=3;i++)
		dest.dst[i] = tw::Int(dest.ijk[i]>space->Dim(i)) - tw::Int(dest.ijk[i]<1);
	for (tw::Int i=0;i<=3;i++)
		dest.x[i] = src.q.x[i];
	dest.p = src.p;
	dest.s = src.s;
	dest.number = src.number;
	dest.tag = src.tag;
	transfer->push_back(dest);
}

template <class MoverType>
void Mover::CopyBack(MoverType *b)
{
	for (int i=0;i<b->num;i++)
	{
		b->refs[i]->q.cell = b->cell[i];
		b->refs[i]->q.x[0] = b->x[0][i];
		b->refs[i]->q.x[1] = b->x[1][i];
		b->refs[i]->q.x[2] = b->x[2][i];
		b->refs[i]->q.x[3] = b->x[3][i];
		b->refs[i]->p[0] = b->u[0][i]*(sqr(m0)+tw::tiny)/(m0+tw::tiny);
		b->refs[i]->p[1] = b->u[1][i]*(sqr(m0)+tw::tiny)/(m0+tw::tiny);
		b->refs[i]->p[2] = b->u[2][i]*(sqr(m0)+tw::tiny)/(m0+tw::tiny);
		b->refs[i]->p[3] = b->u[3][i]*(sqr(m0)+tw::tiny)/(m0+tw::tiny);
		// If particle left MPI domain, add to transfer list, and mark for disposal
		if (b->domainMask[i]==0.0)
		{
			#pragma omp critical
			{
				AddTransferParticle(*(b->refs[i]));
			}
			// mark for disposal only after copying to transfer list
			b->refs[i]->number = 0.0;
		}
		// else
		// {
		// 	tw::Int ijk[4];
		// 	space->DecodeCell(b->refs[i]->q,ijk);
		// 	ASSERT_GTREQ(ijk[1],1);
		// 	ASSERT_LESSEQ(ijk[1],space->Dim(1));
		// 	ASSERT_GTREQ(ijk[2],1);
		// 	ASSERT_LESSEQ(ijk[2],space->Dim(2));
		// 	ASSERT_GTREQ(ijk[3],1);
		// 	ASSERT_LESSEQ(ijk[3],space->Dim(3));
		// }
	}
}

/// return a <beg,end> pair of node5, the internal axis will select component 0.
std::pair<tw::node5,tw::node5> Mover::GetSubarrayBounds(std::vector<ParticleRef>& sorted,tw::Int layers)
{
	// Assumes particles are sorted in increasing memory order
	// This depends in turn on the cell encoding respecting memory order

	tw::Int low[4],high[4];
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
	return std::pair<tw::node5,tw::node5>(
		tw::node5 {low[0],low[1],low[2],low[3],0},
		tw::node5 {high[0]+1,high[1]+1,high[2]+1,high[3]+1,1}
	);
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

template <class MoverType>
void Mover::MoveSlice(tw::Int tasks,tw::Int tid,tw::Int bounds_data[][8])
{
	tw::Int first,last,next;
	std::vector<ParticleRef> map;
	MoverType b(GetParams());
	first = bounds_data[tid][0];
	last = bounds_data[tid][1];
	map.reserve(last-first+1);
	for (tw::Int i=first;i<=last;i++)
	{
		assert(space->IsRefCellWithin((*particle)[i].q,2));
		map.push_back(ParticleRef(i,(*particle)[i]));
	}
	std::sort(map.begin(),map.end());
	auto bounds = GetSubarrayBounds(map,1);
	b.LoadFieldSlice(bounds.first,bounds.second);
	bounds = GetSubarrayBounds(map,2);
	b.InitSourceSlice(bounds.first,bounds.second);
	// Save the bounds information
	for (tw::Int i=1;i<=3;i++)
	{
		bounds_data[tid][i*2] = bounds.first[i];
		bounds_data[tid][i*2+1] = bounds.second[i];
	}
	for (tw::Int i=first;i<=last;i++)
	{
		b.Append((*particle)[map[i-first].idx]);
		next = i==last ? i : i+1;
		if (i==last || b.Complete((*particle)[map[next-first].idx]))
		{
			b.Move(space->dx(0));
			CopyBack(&b);
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

template <class MoverType>
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
	assert(concurrent_tasks<=128);
	tw::Int bounds_data[128][8];
	for (tw::Int c=0;c<total_sets;c++)
	{
		tw::Int tasks_in_set = c<concurrent_sets ? concurrent_tasks : remainder_tasks;
		#pragma omp parallel num_threads(tasks_in_set)
		{
			tw::Int t = tw::GetOMPThreadNum();
			tw::Int task_idx = task_map[c*concurrent_tasks + t];
			tw::GetOMPTaskLoopRange(task_idx,num_par,num_tasks,&bounds_data[t][0],&bounds_data[t][1]);
			MoveSlice<MoverType>(tasks_in_set,t,bounds_data);
		}
	}
}

