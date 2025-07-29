module;

#include "tw_includes.h"

module mover;

void PhotonMover::Advance()
{
	if (space->Dim(2)==1)
		DoTasks<BundleMoverPhoton2D>();
	else
		DoTasks<BundleMoverPhoton3D>();
}

void BundleMoverPhoton2D::LoadFieldSlice(tw::node5& beg,tw::node5& end)
{
	beg[4] = 0; end[4] = 6;
	Fx.Resize(beg,end);
	mov.EM->LoadDataIntoSlice<float>(&Fx);
}

void BundleMoverPhoton2D::InitSourceSlice(tw::node5& beg,tw::node5& end)
{
}

void BundleMoverPhoton2D::DepositSourceSlice(bool needsAtomic)
{
}

void BundleMoverPhoton2D::Move(tw::Float dts)
{
	PrepareGather();
	Push(dts);
	PrepareScatter();
}

void BundleMoverPhoton3D::LoadFieldSlice(tw::node5& beg,tw::node5& end)
{
	beg[4] = 0; end[4] = 6;
	Fx.Resize(beg,end);
	mov.EM->LoadDataIntoSlice<float>(&Fx);
}

void BundleMoverPhoton3D::InitSourceSlice(tw::node5& beg,tw::node5& end)
{
}

void BundleMoverPhoton3D::DepositSourceSlice(bool needsAtomic)
{
}

void BundleMoverPhoton3D::Move(tw::Float dts)
{
	PrepareGather();
	Push(dts);
	PrepareScatter();
}
