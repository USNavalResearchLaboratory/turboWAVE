module;

#include "tw_includes.h"

module mover;

void PGCMover::Advance()
{
	if (space->Dim(2)==1)
		DoTasks<BundleMoverPGC2D>();
	else
		DoTasks<BundleMoverPGC3D>();
}

void BundleMoverPGC2D::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	mov.EM->LoadDataIntoSlice<float>(&Fx);
	lasx.Resize(Element(0,7),low,high,ignorable);
	mov.laser->LoadDataIntoSlice<float>(&lasx);
}

void BundleMoverPGC2D::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	Jx = 0.0f;
	chix.Resize(Element(0),low,high,ignorable);
	chix = 0.0f;
}

void BundleMoverPGC2D::DepositSourceSlice(bool needsAtomic)
{
	if (needsAtomic) {
		mov.sources->AddDataFromSliceAtomic<float>(&Jx);
		mov.chi->AddDataFromSliceAtomic<float>(&chix);
	} else {
		mov.sources->AddDataFromSlice<float>(&Jx);
		mov.chi->AddDataFromSlice<float>(&chix);
	}
}

void BundleMoverPGC2D::Move(tw::Float dts)
{
	const float qmdth = 0.5*mov.q0*dts/mov.m0;
	const tw::Float dti = 1.0/dts;
	const float q2m2dth = 0.5*dts*sqr(mov.q0/mov.m0);
	const float q2m2h = 0.5*sqr(mov.q0/mov.m0);
	PrepareGather();
	LoadFTile(Fx);
	GatherF(F,w0,l0,qmdth);
	LoadLaserTile(lasx);
	GatherLaser(las,w0,q2m2dth,q2m2h);
	Push(dts);
	PrepareScatter();
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile(Jx);
	ResetChiTile();
	ScatterChi(chi,w0,w1,cellMask);
	StoreChiTile(chix);
}

void BundleMoverPGC3D::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	mov.EM->LoadDataIntoSlice<float>(&Fx);
	lasx.Resize(Element(0,7),low,high,ignorable);
	mov.laser->LoadDataIntoSlice<float>(&lasx);
}

void BundleMoverPGC3D::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	Jx = 0.0f;
	chix.Resize(Element(0),low,high,ignorable);
	chix = 0.0f;
}

void BundleMoverPGC3D::DepositSourceSlice(bool needsAtomic)
{
	if (needsAtomic) {
		mov.sources->AddDataFromSliceAtomic<float>(&Jx);
		mov.chi->AddDataFromSliceAtomic<float>(&chix);
	} else {
		mov.sources->AddDataFromSlice<float>(&Jx);
		mov.chi->AddDataFromSlice<float>(&chix);
	}
}

void BundleMoverPGC3D::Move(tw::Float dts)
{
	const float qmdth = 0.5*mov.q0*dts/mov.m0;
	const tw::Float dti = 1.0/dts;
	const float q2m2dth = 0.5*dts*sqr(mov.q0/mov.m0);
	const float q2m2h = 0.5*sqr(mov.q0/mov.m0);
	PrepareGather();
	LoadFTile(Fx);
	GatherF(F,w0,l0,qmdth);
	LoadLaserTile(lasx);
	GatherLaser(las,w0,q2m2dth,q2m2h);
	Push(dts);
	PrepareScatter();
	ResetJTile();
	ScatterJ4(J,w0,w1,cellMask,dti);
	StoreJTile(Jx);
	ResetChiTile();
	ScatterChi(chi,w0,w1,cellMask);
	StoreChiTile(chix);
}

