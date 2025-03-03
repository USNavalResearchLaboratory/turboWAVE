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

void BundleMoverPhoton2D::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	mov.EM->LoadDataIntoImage<float>(&Fx);
}

void BundleMoverPhoton2D::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
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

void BundleMoverPhoton3D::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	mov.EM->LoadDataIntoImage<float>(&Fx);
}

void BundleMoverPhoton3D::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
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
