module;

#include "tw_includes.h"

module mover;

void HCMover::Advance()
{
	if (space->Dim(2)==1)
		DoTasks<BundleMoverHC2D>();
	else
		DoTasks<BundleMoverHC3D>();
}

void BundleMoverHC2D::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	mov.EM->LoadDataIntoImage<float>(&Fx);
}

void BundleMoverHC2D::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	Jx = 0.0f;
}

void BundleMoverHC2D::DepositSourceSlice(bool needsAtomic)
{
	if (needsAtomic) {
		mov.sources->AddDataFromImageAtomic<float>(&Jx);
	} else {
		mov.sources->AddDataFromImage<float>(&Jx);
	}
}

void BundleMoverHC2D::Move(tw::Float dts)
{
	const float qmdth = 0.5*mov.q0*dts/mov.m0;
	const tw::Float dti = 1.0/dts;
	PrepareGather();
	LoadFTile(Fx);
	GatherF(F,w0,l0,qmdth);
	Push(dts);
	PrepareScatter();
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile(Jx);
}

void BundleMoverHC3D::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	mov.EM->LoadDataIntoImage<float>(&Fx);
}

void BundleMoverHC3D::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	Jx = 0.0f;
}

void BundleMoverHC3D::DepositSourceSlice(bool needsAtomic)
{
	if (needsAtomic) {
		mov.sources->AddDataFromImageAtomic<float>(&Jx);
	} else {
		mov.sources->AddDataFromImage<float>(&Jx);
	}
}

void BundleMoverHC3D::Move(tw::Float dts)
{
	const float qmdth = 0.5*mov.q0*dts/mov.m0;
	const tw::Float dti = 1.0/dts;
	PrepareGather();
	LoadFTile(Fx);
	GatherF(F,w0,l0,qmdth);
	Push(dts);
	PrepareScatter();
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile(Jx);
}

