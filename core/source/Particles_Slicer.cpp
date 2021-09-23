#include "meta_base.h"
#include "computeTool.h"
#include "particles_bundle.h"
#include "particles_pusher.h"
#include "particles_slicer.h"
#include "particles_tiler.h"
#include "particles_mover.h"

void BundleSlicerEM::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	owner->EM->LoadDataIntoImage<float>(&Fx);
}
void BundleSlicerEM::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	Jx = 0.0f;
}
void BundleSlicerEM::DepositSourceSlice(bool needsAtomic)
{
	if (needsAtomic)
		owner->sources->AddDataFromImageAtomic<float>(&Jx);
	else
		owner->sources->AddDataFromImage<float>(&Jx);
}

void BundleSlicerPGC::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	lasx.Resize(Element(0,7),low,high,ignorable);
	owner->laser->LoadDataIntoImage<float>(&lasx);
}
void BundleSlicerPGC::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	chix.Resize(Element(0),low,high,ignorable);
	chix = 0.0f;
}
void BundleSlicerPGC::DepositSourceSlice(bool needsAtomic)
{
	if (needsAtomic)
		owner->chi->AddDataFromImageAtomic<float>(&chix);
	else
		owner->chi->AddDataFromImage<float>(&chix);
}

void BundleSlicerBohmian::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	owner->qo_j4->LoadDataIntoImage<float>(&Jx);
}
void BundleSlicerBohmian::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
}
void BundleSlicerBohmian::DepositSourceSlice(bool needsAtomic)
{
}
