struct BundleSlicerEM : virtual ParticleBundle
{
	Slice<float> Fx,Jx;
	BundleSlicerEM(Mover *owner) : ParticleBundle(owner) {}
	void LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void DepositSourceSlice(bool needsAtomic);
};

struct BundleSlicerPGC : virtual ParticleBundle
{
	Slice<float> lasx,chix;
	BundleSlicerPGC(Mover *owner) : ParticleBundle(owner) {}
	void LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void DepositSourceSlice(bool needsAtomic);
};

struct BundleSlicerBohmian : virtual ParticleBundle
{
	Slice<float> Jx;
	BundleSlicerBohmian(Mover *owner) : ParticleBundle(owner) {}
	void LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void DepositSourceSlice(bool needsAtomic);
};
