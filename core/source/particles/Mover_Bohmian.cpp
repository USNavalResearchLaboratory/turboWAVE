module;

#include "tw_includes.h"

module mover;

void BohmianMover::Advance()
{
	if (space->Dim(2)==1)
		DoTasks<BundleMoverBohmian2D>();
	else
		DoTasks<BundleMoverBohmian3D>();
}

void BundleMoverBohmian2D::LoadFieldSlice(tw::node5& beg,tw::node5& end)
{
	beg[4] = 0; end[4] = 4;
	Jx.Resize(beg,end);
	mov.qo_j4->LoadDataIntoSlice<float>(&Jx);
}

void BundleMoverBohmian2D::InitSourceSlice(tw::node5& beg,tw::node5& end)
{
}

void BundleMoverBohmian2D::DepositSourceSlice(bool needsAtomic)
{
}

void BundleMoverBohmian2D::Move(tw::Float dts)
{
	PadBundle();
	cell0 = cell[0];
	mov.ms->DecodeCell(cell0,ijk0);
	mov.ms->StaticSpace::GetWeights(w0,x);
	LoadTile(Jx);
	GatherJ4(J,w0);
	Push(dts);
	mov.ms->MinimizePrimitive(cell,ijk,x,domainMask);
}

void BundleMoverBohmian3D::LoadFieldSlice(tw::node5& beg,tw::node5& end)
{
	beg[4] = 0; end[4] = 4;
	Jx.Resize(beg,end);
	mov.qo_j4->LoadDataIntoSlice<float>(&Jx);
}

void BundleMoverBohmian3D::InitSourceSlice(tw::node5& beg,tw::node5& end)
{
}

void BundleMoverBohmian3D::DepositSourceSlice(bool needsAtomic)
{
}

void BundleMoverBohmian3D::Move(tw::Float dts)
{
	PadBundle();
	cell0 = cell[0];
	mov.ms->DecodeCell(cell0,ijk0);
	mov.ms->StaticSpace::GetWeights(w0,x);
	LoadTile(Jx);
	GatherJ4(J,w0);
	Push(dts);
	mov.ms->MinimizePrimitive(cell,ijk,x,domainMask);
}

