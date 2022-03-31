#include "meta_base.h"
#include "computeTool.h"
#include "bundle.h"
#include "pusher.h"
#include "slicer.h"
#include "tiler.h"
#include "mover.h"

ParticleBundle::ParticleBundle(Mover *owner)
{
	this->owner = owner;
	q0 = owner->q0;
	m0 = owner->m0;
	dt = timestep(*owner->space);
	k[0] = dxi(*owner->space);
	k[1] = dyi(*owner->space);
	k[2] = dzi(*owner->space);
	refs.resize(N);
	num = 0;
}

void ParticleBundle::CopyBack()
{
	for (int i=0;i<num;i++)
	{
		refs[i]->q.cell = cell[i];
		refs[i]->q.x[0] = x[0][i];
		refs[i]->q.x[1] = x[1][i];
		refs[i]->q.x[2] = x[2][i];
		refs[i]->q.x[3] = x[3][i];
		refs[i]->p[0] = u[0][i]*m0;
		refs[i]->p[1] = u[1][i]*m0;
		refs[i]->p[2] = u[2][i]*m0;
		refs[i]->p[3] = u[3][i]*m0;
		// If particle left MPI domain, add to transfer list, and mark for disposal
		if (domainMask[i]==0.0)
		{
			#pragma omp critical
			{
				owner->AddTransferParticle(*refs[i]);
			}
			// mark for disposal only after copying to transfer list
			refs[i]->number = 0.0;
		}
		// else
		// {
		// 	tw::Int ijk[4];
		// 	owner->space->DecodeCell(refs[i]->q,ijk);
		// 	ASSERT_GTREQ(ijk[1],1);
		// 	ASSERT_LESSEQ(ijk[1],owner->space->Dim(1));
		// 	ASSERT_GTREQ(ijk[2],1);
		// 	ASSERT_LESSEQ(ijk[2],owner->space->Dim(2));
		// 	ASSERT_GTREQ(ijk[3],1);
		// 	ASSERT_LESSEQ(ijk[3],owner->space->Dim(3));
		// }
	}
}

void ParticleBundle::PrepareGather()
{
	PadBundle();
	cell0 = cell[0];
	owner->space->DecodeCell(cell0,ijk0);
	owner->space->GetWeights(w0,x);
	owner->space->GetWallWeights(l0,x);
}

void ParticleBundle::PrepareScatter()
{
	owner->space->MinimizePrimitive(cell,ijk,x,domainMask);
	owner->space->GetWeights(w1,x);
	set_cell_mask(cellMask,cell0,cell);
}
