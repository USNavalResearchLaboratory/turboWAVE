struct Pusher:ComputeTool
{
	tw::Float q0,m0;
	tw::Int ignorable[4];
	// Following particles and fields are owned by the parent Species module
	std::vector<Particle> *particle;
	std::vector<TransferParticle> *transfer;
	Field *EM,*sources,*laser,*chi,*qo_j4;
	Pusher(const std::string& name,MetricSpace *m,Task *tsk);
	void AddTransferParticle(const Particle& src);
	void GetSubarrayBounds(std::vector<ParticleRef>& sorted,tw::Int low[4],tw::Int high[4],tw::Int layers);
	void SpreadTasks(std::vector<tw::Int>& task_map);
	void BunchTasks(std::vector<tw::Int>& task_map);
	template <class BundleType>
	void PushSlice(tw::Int tasks,tw::Int tid,tw::Int bounds_data[][8]);
	template <class BundleType>
	void Push();
	virtual void Advance();
};

struct BorisPusher2D:Pusher
{
	BorisPusher2D(const std::string& name,MetricSpace *m,Task *tsk): Pusher(name,m,tsk) {}
	virtual void Advance();
};

struct BorisPusher3D:Pusher
{
	BorisPusher3D(const std::string& name,MetricSpace *m,Task *tsk): Pusher(name,m,tsk) {}
	virtual void Advance();
};

struct UnitaryPusher2D:Pusher
{
	UnitaryPusher2D(const std::string& name,MetricSpace *m,Task *tsk): Pusher(name,m,tsk) {}
	virtual void Advance();
};

struct UnitaryPusher3D:Pusher
{
	UnitaryPusher3D(const std::string& name,MetricSpace *m,Task *tsk): Pusher(name,m,tsk) {}
	virtual void Advance();
};

struct PGCPusher:Pusher
{
	PGCPusher(const std::string& name,MetricSpace *m,Task *tsk): Pusher(name,m,tsk) {}
	virtual void Advance();
};

struct BohmianPusher:Pusher
{
	BohmianPusher(const std::string& name,MetricSpace *m,Task *tsk): Pusher(name,m,tsk) {}
	virtual void Advance();
};

struct ParticleBundle
{
	// The purpose of the particle bundle is to allow for efficient vectorization of the pusher
	// Several inlined functions are needed due to OpenMP inability to work with member variables
	// Evidently this is related to uninstantiated variables being tagged as aligned
	// The workaround is to pass members back in as arguments
	Pusher *owner;
	tw::Float q0,m0,dt,k[3];
	tw::Int num,cell0,ijk0[3];
	static const tw::Int N = tw::max_bundle_size;
	static const tw::Int AB = tw::vec_align_bytes;
	// arrays storing particle state
	alignas(AB) tw::Int cell[N];
	alignas(AB) float x[4][N];
	alignas(AB) tw::Float u[4][N];
	alignas(AB) float number[N];
	// temporary arrays used during calculation
	alignas(AB) float J[4][N];
	alignas(AB) float domainMask[N];
	alignas(AB) float cellMask[N];
	alignas(AB) tw::Int ijk[3][N];
	alignas(AB) tw::Float vel[3][N];
	alignas(AB) float w0[3][3][N];
	alignas(AB) float w1[3][3][N];
	alignas(AB) float l0[3][3][N];

	std::valarray<Particle*> refs;

	ParticleBundle(Pusher *owner);
	void PadBundle();
	void Reset();
	bool Complete(const Particle& par);
	void Append(Particle& par);
	void CopyBack();
	void translate(float x[4][N],tw::Float vel[3][N]);
	void set_cell_mask(float cellMask[N],tw::Int cell0,tw::Int cell[N]);
	void get_cell_displ(tw::Int dc[3],tw::Int n);
	void load_j4(float J[4][N],const float number[N],const tw::Float vel[3][N]);
	void ZeroArray(float q[][N],tw::Int s1,tw::Int s2);
};

struct BundlePusherBoris : virtual ParticleBundle
{
	alignas(AB) float F[6][N]; // q*F*dth/m
	alignas(AB) tw::Float t[3][N];
	alignas(AB) tw::Float s[3][N];

	BundlePusherBoris(Pusher *owner) : ParticleBundle(owner) { ; }
	void impulse(tw::Float u[4][N],float F[6][N]);
	void rotation1(tw::Float t[3][N],tw::Float u[4][N],float F[6][N]);
	void rotation2(tw::Float s[3][N],tw::Float t[3][N]);
	void rotation3(tw::Float s[3][N],tw::Float t[3][N],tw::Float vel[3][N],tw::Float u[4][N]);
	void velocity(tw::Float vel[3][N],tw::Float u[4][N]);
	void Push();
};

struct BundlePusherPGC : BundlePusherBoris
{
	alignas(AB) float las[8][N]; // q*q*a*a/m*m
	alignas(AB) float chi[N];
	alignas(AB) tw::Float avgGam[N];

	BundlePusherPGC(Pusher *owner) : BundlePusherBoris(owner) , ParticleBundle(owner) { ; }
	void avg_gam_1(tw::Float avgGam[N],tw::Float vel[3][N],tw::Float u[4][N],float F[6][N],float las[8][N]);
	void avg_gam_2(tw::Float avgGam[N],tw::Float u[4][N],float las[8][N]);
	void impulse(tw::Float u[4][N],float F[6][N],float las[8][N],tw::Float avgGam[N]);
	void rotation1(tw::Float t[3][N],float F[6][N],tw::Float avgGam[N]);
	void velocity(tw::Float vel[3][N],tw::Float u[4][N],tw::Float avgGam[N]);
	void load_chi(float chi[N],float number[N],tw::Float avgGam[N]);
	void Push();
};

struct BundlePusherUnitary : virtual ParticleBundle
{
	alignas(AB) float ds[N];
	alignas(AB) float F[6][N];
	alignas(AB) tw::Float a[3][N];
	alignas(AB) tw::Float z[2][2][2][N];
	alignas(AB) tw::Float zi[2][2][2][N];
	alignas(AB) tw::Float zf[2][2][2][N];
	BundlePusherUnitary(Pusher *owner) : ParticleBundle(owner) { ; }
	void estimate_ds(float ds[N],float F[6][N],tw::Float u[4][N]);
	void copy_spinor(tw::Float dst[2][2][2][N],tw::Float src[2][2][2][N]);
	void add_spinor(tw::Float dst[2][2][2][N],tw::Float src[2][2][2][N]);
	void to_spinor(tw::Float u[4][N],tw::Float z[2][2][2][N]);
	void to_vector(tw::Float u[4][N],tw::Float z[2][2][2][N],tw::Float a[3][N]);
	void dagger(tw::Float z[2][2][2][N]);
	void scalar_mul(tw::Float a[3][N],tw::Float z[2][2][2][N]);
	void left_mul_sig1(tw::Float z[2][2][2][N]);
	void left_mul_sig2(tw::Float z[2][2][2][N]);
	void left_mul_sig3(tw::Float z[2][2][2][N]);
	void set_psi_1(tw::Float a[3][N],float F[6][N]);
	void set_psi_2(tw::Float a[3][N],float F[6][N]);
	void set_psi_3(tw::Float a[3][N],float F[6][N]);
	void set_psi24(tw::Float a[3][N],float F[6][N]);
	void velocity(tw::Float vel[3][N],tw::Float u[4][N]);
	void Lambda();
	void Push();
};

struct BundlePusherBohmian : virtual ParticleBundle
{
	BundlePusherBohmian(Pusher *owner) : ParticleBundle(owner) { ; }
	void bohm_velocity(tw::Float vel[3][N],tw::Float u[4][N],float J[4][N]);
	void Push();
};

struct BundleTilerEM : virtual ParticleBundle
{
	Slice<float> Fx,Jx;
	BundleTilerEM(Pusher *owner) : ParticleBundle(owner) { ; }
	void LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void DepositSourceSlice(bool needsAtomic);
};

struct BundleTiler2D : BundleTilerEM
{
	float F_tile[3][3][6];
	float J_tile[5][5][4];

	BundleTiler2D(Pusher *owner) : BundleTilerEM(owner) , ParticleBundle(owner) { ; }
	void LoadFTile();
	void ResetJTile();
	void StoreJTile();
	void GatherF(float F[6][N],const float w0[3][3][N],const float l0[3][3][N],const float qmdth);
	void ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti);
	void Gather(float F[6][N]);
	void Scatter();
};

struct BundleTiler3D : BundleTilerEM
{
	float F_tile[3][3][3][6];
	float J_tile[5][5][5][4];

	BundleTiler3D(Pusher *owner) : BundleTilerEM(owner) , ParticleBundle(owner) { ; }
	void LoadFTile();
	void ResetJTile();
	void StoreJTile();
	void GatherF(float F[6][N],const float w0[3][3][N],const float l0[3][3][N],const float qmdth);
	void ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti);
	void Gather(float F[6][N]);
	void Scatter();
};

struct BundleTilerPGC : BundleTiler3D
{
	Slice<float> lasx,chix;
	float las_tile[3][3][3][8];
	float chi_tile[5][5][5];

	BundleTilerPGC(Pusher *owner) : BundleTiler3D(owner) , ParticleBundle(owner) { ; }
	void LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void DepositSourceSlice(bool needsAtomic);
	void LoadLaserTile();
	void ResetChiTile();
	void StoreChiTile();
	void GatherLaser(float las[8][N],const float w0[3][3][N],const float q2m2);
	void ScatterChi(const float chi[N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N]);
	void Gather(float F[6][N],float las[8][N]);
	void Scatter(float chi[N]);
};

struct BundleTilerBohmian : virtual ParticleBundle
{
	Slice<float> Jx;
	float tile[3][3][3][4];
	BundleTilerBohmian(Pusher *owner) : ParticleBundle(owner) { ; }
	void LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void DepositSourceSlice(bool needsAtomic);
	void LoadTile();
	void GatherJ4(float J[4][N],const float w0[3][3][N]);
	void Gather();
};

struct ParticleBundle2D : BundleTiler2D,BundlePusherBoris
{
	ParticleBundle2D(Pusher *owner) : BundleTiler2D(owner),BundlePusherBoris(owner) , ParticleBundle(owner) { ; }
	void Advance();
};

struct ParticleBundle3D : BundleTiler3D,BundlePusherBoris
{
	ParticleBundle3D(Pusher *owner) : BundleTiler3D(owner),BundlePusherBoris(owner) , ParticleBundle(owner) { ; }
	void Advance();
};

struct ParticleBundlePGC : BundleTilerPGC,BundlePusherPGC
{
	ParticleBundlePGC(Pusher *owner) : BundleTilerPGC(owner),BundlePusherPGC(owner) , ParticleBundle(owner) { ; }
	void Advance();
};

struct ParticleBundleBohmian : BundleTilerBohmian,BundlePusherBohmian
{
	ParticleBundleBohmian(Pusher *owner) : BundleTilerBohmian(owner),BundlePusherBohmian(owner) , ParticleBundle(owner) { ; }
	void Advance();
};

struct ParticleBundleUnitary2D : BundleTiler2D,BundlePusherUnitary
{
	ParticleBundleUnitary2D(Pusher *owner) : BundleTiler2D(owner),BundlePusherUnitary(owner) , ParticleBundle(owner) { ; }
	void Advance();
};

struct ParticleBundleUnitary3D : BundleTiler3D,BundlePusherUnitary
{
	ParticleBundleUnitary3D(Pusher *owner) : BundleTiler3D(owner),BundlePusherUnitary(owner) , ParticleBundle(owner) { ; }
	void Advance();
};

//////////////////////////
//                      //
// INLINE PUSH ROUTINES //
//                      //
//////////////////////////


inline ParticleBundle::ParticleBundle(Pusher *owner)
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
inline void ParticleBundle::PadBundle()
{
	for (int ax=0;ax<4;ax++)
		for (int i=num;i<N;i++)
		{
			x[ax][i] = 0.0;
			u[ax][i] = 0.0;
		}
	for (int i=num;i<N;i++)
	{
		u[0][i] = 1.0; // u needs to be a good 4-vector
		cell[i] = cell[0];
		number[i] = 0.0;
	}
}
inline void ParticleBundle::Reset()
{
	num = 0;
}
inline bool ParticleBundle::Complete(const Particle& par)
{
	return (num && par.q.cell!=cell[0]) || num==N;
}
inline void ParticleBundle::Append(Particle& par)
{
	refs[num] = &par;
	cell[num] = par.q.cell;
	x[1][num] = par.q.x[0];
	x[2][num] = par.q.x[1];
	x[3][num] = par.q.x[2];
	u[1][num] = par.p[0]/m0;
	u[2][num] = par.p[1]/m0;
	u[3][num] = par.p[2]/m0;
	number[num] = par.number;
	num++;
}
inline void ParticleBundle::CopyBack()
{
	for (int i=0;i<num;i++)
	{
		refs[i]->q.cell = cell[i];
		refs[i]->q.x[0] = x[1][i];
		refs[i]->q.x[1] = x[2][i];
		refs[i]->q.x[2] = x[3][i];
		refs[i]->p[0] = u[1][i]*m0;
		refs[i]->p[1] = u[2][i]*m0;
		refs[i]->p[2] = u[3][i]*m0;
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
	}
}
inline void ParticleBundle::translate(float x[4][N],tw::Float vel[3][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(x,vel:AB)
		for (int i=0;i<N;i++)
			x[c+1][i] += vel[c][i]*dt*k[c];
}
inline void ParticleBundle::set_cell_mask(float cellMask[N],tw::Int cell0,tw::Int cell[N])
{
	for (int i=0;i<N;i++)
		cellMask[i] = float(cell[i]==cell0);
}
inline void ParticleBundle::get_cell_displ(tw::Int dc[3],tw::Int n)
{
	dc[0] = ijk[0][n] - ijk0[0];
	dc[1] = ijk[1][n] - ijk0[1];
	dc[2] = ijk[2][n] - ijk0[2];
}
inline void ParticleBundle::load_j4(float J[4][N],const float number[N],const tw::Float vel[3][N])
{
	#pragma omp simd aligned(J,number:AB)
	for (int i=0;i<N;i++)
		J[0][i] = q0*number[i];
	#pragma omp simd aligned(J,number,vel:AB)
	for (int i=0;i<N;i++)
		J[1][i] = q0*number[i]*vel[0][i]*k[0];
	#pragma omp simd aligned(J,number,vel:AB)
	for (int i=0;i<N;i++)
		J[2][i] = q0*number[i]*vel[1][i]*k[1];
	#pragma omp simd aligned(J,number,vel:AB)
	for (int i=0;i<N;i++)
		J[3][i] = q0*number[i]*vel[2][i]*k[2];
}
inline void ParticleBundle::ZeroArray(float q[][N],tw::Int s1,tw::Int s2)
{
	for (tw::Int s=s1;s<=s2;s++)
		for (tw::Int par=0;par<N;par++)
			q[s][par] = 0.0f;
}


///////////////////////////////////////////
// UNITARY PUSHER
//////////////////////////////////////////


inline void BundlePusherUnitary::estimate_ds(float ds[N],float F[6][N],tw::Float u[4][N])
{
	// ds is the ratio of proper time step to lab step
	// the dilation factor is folded into F (which already contains dt)
	#pragma omp simd aligned(ds,F,u:AB)
	for (int i=0;i<N;i++)
	{
		ds[i] = 1.0/u[0][i];
		ds[i] *= 1.0 - (F[0][i]*u[1][i] + F[1][i]*u[2][i] + F[2][i]*u[3][i])*dt/sqr(u[0][i]);
	}
	#pragma omp simd aligned(ds,F:AB)
	for (int i=0;i<N;i++)
	{
		F[0][i] *= ds[i];
		F[1][i] *= ds[i];
		F[2][i] *= ds[i];
		F[3][i] *= ds[i];
		F[4][i] *= ds[i];
		F[5][i] *= ds[i];
	}
}
inline void BundlePusherUnitary::copy_spinor(tw::Float dst[2][2][2][N],tw::Float src[2][2][2][N])
{
	#pragma omp simd aligned(dst,src:AB)
	for (int i=0;i<N;i++)
	{
		dst[0][0][0][i] = src[0][0][0][i];
		dst[0][0][1][i] = src[0][0][1][i];
		dst[0][1][0][i] = src[0][1][0][i];
		dst[0][1][1][i] = src[0][1][1][i];
		dst[1][0][0][i] = src[1][0][0][i];
		dst[1][0][1][i] = src[1][0][1][i];
		dst[1][1][0][i] = src[1][1][0][i];
		dst[1][1][1][i] = src[1][1][1][i];
	}
}
inline void BundlePusherUnitary::add_spinor(tw::Float dst[2][2][2][N],tw::Float src[2][2][2][N])
{
	#pragma omp simd aligned(dst,src:AB)
	for (int i=0;i<N;i++)
	{
		dst[0][0][0][i] += src[0][0][0][i];
		dst[0][0][1][i] += src[0][0][1][i];
		dst[0][1][0][i] += src[0][1][0][i];
		dst[0][1][1][i] += src[0][1][1][i];
		dst[1][0][0][i] += src[1][0][0][i];
		dst[1][0][1][i] += src[1][0][1][i];
		dst[1][1][0][i] += src[1][1][0][i];
		dst[1][1][1][i] += src[1][1][1][i];
	}
}
inline void BundlePusherUnitary::to_spinor(tw::Float u[4][N],tw::Float z[2][2][2][N])
{
	// Form the spinor, z, from the four velocity, u
	#pragma omp simd aligned(u,z:AB)
	for (int i=0;i<N;i++)
	{
		// need the following line because u[0] is not part of the global particle state
		u[0][i] = sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i]);
		// real part
		z[0][0][0][i] = u[0][i] + u[3][i];
		z[0][1][0][i] = u[1][i];
		z[1][0][0][i] = u[1][i];
		z[1][1][0][i] = u[0][i] - u[3][i];
		// imag part
		z[0][0][1][i] = 0.0f;
		z[0][1][1][i] = -u[2][i];
		z[1][0][1][i] = u[2][i];
		z[1][1][1][i] = 0.0f;
	}
}
inline void BundlePusherUnitary::to_vector(tw::Float u[4][N],tw::Float z[2][2][2][N],tw::Float a[3][N])
{
	// Restore the four velocity, u, from the spinor, z
	set_psi24(a,F);
	#pragma omp simd aligned(u,z:AB)
	for (int i=0;i<N;i++)
	{
		a[0][i] = 1.0/((1.0 - a[0][i])*(1.0 - a[0][i]) + a[1][i]*a[1][i]);
		u[0][i] = a[0][i]*0.5f*(z[0][0][0][i] + z[1][1][0][i]);
		u[1][i] = a[0][i]*z[0][1][0][i];
		u[2][i] = a[0][i]*z[1][0][1][i];
		u[3][i] = a[0][i]*0.5f*(z[0][0][0][i] - z[1][1][0][i]);
	}
}
inline void BundlePusherUnitary::dagger(tw::Float z[2][2][2][N])
{
	// Form the Hermitian conjugate of z in place
	#pragma omp simd aligned(z:AB)
	for (int i=0;i<N;i++)
	{
		std::swap(z[0][1][0][i],z[1][0][0][i]);
		std::swap(z[0][1][1][i],z[1][0][1][i]);
		z[0][0][1][i] *= -1.0f;
		z[0][1][1][i] *= -1.0f;
		z[1][0][1][i] *= -1.0f;
		z[1][1][1][i] *= -1.0f;
	}
}
inline void BundlePusherUnitary::scalar_mul(tw::Float a[3][N],tw::Float z[2][2][2][N])
{
	// multiply a spinor, z, by a complex scalar, a.
	// the extra element of the scalar is used as a temporary.
	#pragma omp simd aligned(a,z:AB)
	for (int i=0;i<N;i++)
	{
		a[2][i] = z[0][0][0][i];
		z[0][0][0][i] = z[0][0][0][i]*a[0][i] - z[0][0][1][i]*a[1][i];
		z[0][0][1][i] = a[2][i]*a[1][i] + z[0][0][1][i]*a[0][i];

		a[2][i] = z[0][1][0][i];
		z[0][1][0][i] = z[0][1][0][i]*a[0][i] - z[0][1][1][i]*a[1][i];
		z[0][1][1][i] = a[2][i]*a[1][i] + z[0][1][1][i]*a[0][i];

		a[2][i] = z[1][0][0][i];
		z[1][0][0][i] = z[1][0][0][i]*a[0][i] - z[1][0][1][i]*a[1][i];
		z[1][0][1][i] = a[2][i]*a[1][i] + z[1][0][1][i]*a[0][i];

		a[2][i] = z[1][1][0][i];
		z[1][1][0][i] = z[1][1][0][i]*a[0][i] - z[1][1][1][i]*a[1][i];
		z[1][1][1][i] = a[2][i]*a[1][i] + z[1][1][1][i]*a[0][i];
	}
}
inline void BundlePusherUnitary::left_mul_sig1(tw::Float z[2][2][2][N])
{
	// Left multiply by Pauli matrix sigma_1 (swap rows)
	#pragma omp simd aligned(z:AB)
	for (int i=0;i<N;i++)
	{
		std::swap(z[0][0][0][i],z[1][0][0][i]);
		std::swap(z[0][0][1][i],z[1][0][1][i]);
		std::swap(z[0][1][0][i],z[1][1][0][i]);
		std::swap(z[0][1][1][i],z[1][1][1][i]);
	}
}
inline void BundlePusherUnitary::left_mul_sig2(tw::Float z[2][2][2][N])
{
	// Left multiply by Pauli matrix sigma_2
	#pragma omp simd aligned(z:AB)
	for (int i=0;i<N;i++)
	{
		// First multiply top row by i and bottom row by -i
		std::swap(z[0][0][0][i],z[0][0][1][i]);
		z[0][0][0][i] *= -1.0f;
		std::swap(z[0][1][0][i],z[0][1][1][i]);
		z[0][1][0][i] *= -1.0f;
		std::swap(z[1][0][0][i],z[1][0][1][i]);
		z[1][0][1][i] *= -1.0f;
		std::swap(z[1][1][0][i],z[1][1][1][i]);
		z[1][1][1][i] *= -1.0f;
		// Now swap rows
		std::swap(z[0][0][0][i],z[1][0][0][i]);
		std::swap(z[0][0][1][i],z[1][0][1][i]);
		std::swap(z[0][1][0][i],z[1][1][0][i]);
		std::swap(z[0][1][1][i],z[1][1][1][i]);
	}
}
inline void BundlePusherUnitary::left_mul_sig3(tw::Float z[2][2][2][N])
{
	// Left multiply by Pauli matrix sigma_3 (change sign of bottom row)
	#pragma omp simd aligned(z:AB)
	for (int i=0;i<N;i++)
	{
		z[1][0][0][i] *= -1.0f;
		z[1][0][1][i] *= -1.0f;
		z[1][1][0][i] *= -1.0f;
		z[1][1][1][i] *= -1.0f;
	}
}
inline void BundlePusherUnitary::set_psi_1(tw::Float a[3][N],float F[6][N])
{
	#pragma omp simd aligned(a,F:AB)
	for (int i=0;i<N;i++)
	{
		a[0][i] = F[0][i];
		a[1][i] = F[3][i];
	}
}
inline void BundlePusherUnitary::set_psi_2(tw::Float a[3][N],float F[6][N])
{
	#pragma omp simd aligned(a,F:AB)
	for (int i=0;i<N;i++)
	{
		a[0][i] = F[1][i];
		a[1][i] = F[4][i];
	}
}
inline void BundlePusherUnitary::set_psi_3(tw::Float a[3][N],float F[6][N])
{
	#pragma omp simd aligned(a,F:AB)
	for (int i=0;i<N;i++)
	{
		a[0][i] = F[2][i];
		a[1][i] = F[5][i];
	}
}
inline void BundlePusherUnitary::set_psi24(tw::Float a[3][N],float F[6][N])
{
	#pragma omp simd aligned(a,F:AB)
	for (int i=0;i<N;i++)
	{
		a[0][i] = 0.25*(F[0][i]*F[0][i] - F[3][i]*F[3][i] + F[1][i]*F[1][i] - F[4][i]*F[4][i] + F[2][i]*F[2][i] - F[5][i]*F[5][i]);
		a[1][i] = 0.5*(F[0][i]*F[3][i] + F[1][i]*F[4][i] + F[2][i]*F[5][i]);
	}
}
inline void BundlePusherUnitary::velocity(tw::Float vel[3][N],tw::Float u[4][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,u:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = u[c+1][i]/u[0][i];
}

///////////////////////////////////////////
// BORIS PUSHER
//////////////////////////////////////////


inline void BundlePusherBoris::impulse(tw::Float u[4][N],float F[6][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(u,F:AB)
		for (int i=0;i<N;i++)
			u[c+1][i] += F[c][i];
}
inline void BundlePusherBoris::rotation1(tw::Float t[3][N],tw::Float u[4][N],float F[6][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(t,u,F:AB)
		for (int i=0;i<N;i++)
			t[c][i] = F[c+3][i]/sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i]);
}
inline void BundlePusherBoris::rotation2(tw::Float s[3][N],tw::Float t[3][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(s,t:AB)
		for (int i=0;i<N;i++)
			s[c][i] = 2.0*t[c][i]/(1.0 + t[0][i]*t[0][i] + t[1][i]*t[1][i] + t[2][i]*t[2][i]);
}
inline void BundlePusherBoris::rotation3(tw::Float s[3][N],tw::Float t[3][N],tw::Float vel[3][N],tw::Float u[4][N])
{
	// Use vel as temporary while adding rotational impulse, cross(u + cross(u,t),s)
	#pragma omp simd aligned(s,t,vel,u:AB)
	for (int i=0;i<N;i++)
	{
		vel[0][i] = u[1][i] + u[2][i]*t[2][i] - u[3][i]*t[1][i];
		vel[1][i] = u[2][i] + u[3][i]*t[0][i] - u[1][i]*t[2][i];
		vel[2][i] = u[3][i] + u[1][i]*t[1][i] - u[2][i]*t[0][i];
		u[1][i] += vel[1][i]*s[2][i] - vel[2][i]*s[1][i];
		u[2][i] += vel[2][i]*s[0][i] - vel[0][i]*s[2][i];
		u[3][i] += vel[0][i]*s[1][i] - vel[1][i]*s[0][i];
	}
}
inline void BundlePusherBoris::velocity(tw::Float vel[3][N],tw::Float u[4][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,u:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = u[c+1][i]/sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i]);
}


///////////////////////////////////////////
// PGC PUSHER
//////////////////////////////////////////


inline void BundlePusherPGC::avg_gam_1(tw::Float avgGam[N],tw::Float vel[3][N],tw::Float u[4][N],float F[6][N],float las[8][N])
{
	const tw::Float dth = 0.5*dt;
	alignas(AB) tw::Float g1[N];
	alignas(AB) tw::Float dA[N];
	// estimate avgGam at level n using solution of <g>^2 = g1<g> - dA/4
	#pragma omp simd aligned(g1,u,las:AB)
	for (int i=0;i<N;i++)
		g1[i] = sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i] + 0.5*las[7][i]);
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,u,g1:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = u[c+1][i]/g1[i];
	#pragma omp simd aligned(g1,F,vel,dA,las,avgGam:AB)
	for (int i=0;i<N;i++)
	{
		g1[i] += F[0][i]*vel[0][i] + F[1][i]*vel[1][i] + F[2][i]*vel[2][i];
		dA[i] = dth*(las[0][i]*vel[0][i] + las[1][i]*vel[1][i] + las[2][i]*vel[2][i]);
		avgGam[i] = 0.5*g1[i] + 0.5*sqrt(g1[i]*g1[i] - dA[i]);
	}
}
inline void BundlePusherPGC::avg_gam_2(tw::Float avgGam[N],tw::Float u[4][N],float las[8][N])
{
	// estimate avgGam at level n+1/2 using <g>^2 = g1^2 + (g1/<g>)dA/2
	// (solve cubic and expand in powers of dA/g1^2)
	const tw::Float dth = 0.5*dt;
	alignas(AB) tw::Float g1[N];
	alignas(AB) tw::Float dA[N];
	#pragma omp simd aligned(g1,dA,u,las,avgGam:AB)
	for (int i=0;i<N;i++)
	{
		g1[i] = sqrt(1.0 + u[1][i]*u[1][i] + u[2][i]*u[2][i] + u[3][i]*u[3][i] + 0.5*las[7][i]);
		dA[i] = dth*(las[3][i]*u[1][i] + las[4][i]*u[2][i] + las[5][i]*u[3][i])/g1[i];
		avgGam[i] = g1[i] + 0.25*dA[i]/g1[i] - 0.09375*dA[i]*dA[i]/(g1[i]*g1[i]*g1[i]);
	}
}
inline void BundlePusherPGC::impulse(tw::Float u[4][N],float F[6][N],float las[8][N],tw::Float avgGam[N])
{
	const tw::Float dth = 0.5*dt;
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(u,F,las,avgGam:AB)
		for (int i=0;i<N;i++)
			u[c+1][i] += F[c][i] - 0.25*las[c][i]*dth/avgGam[i];
}
inline void BundlePusherPGC::rotation1(tw::Float t[3][N],float F[6][N],tw::Float avgGam[N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(t,F,avgGam:AB)
		for (int i=0;i<N;i++)
			t[c][i] = F[c+3][i]/avgGam[i];
}
inline void BundlePusherPGC::velocity(tw::Float vel[3][N],tw::Float u[4][N],tw::Float avgGam[N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,u,avgGam:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = u[c+1][i]/avgGam[i];
}
inline void BundlePusherPGC::load_chi(float chi[N],float number[N],tw::Float avgGam[N])
{
	#pragma omp simd aligned(chi,number,avgGam:AB)
	for (int i=0;i<N;i++)
		chi[i] = -q0*q0*number[i]/(m0*avgGam[i]);
}


///////////////////////////////////////////
// BOHMIAN PUSHER
//////////////////////////////////////////


inline void BundlePusherBohmian::bohm_velocity(tw::Float vel[3][N],tw::Float u[4][N],float J[4][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,u,J:AB)
		for (int i=0;i<N;i++)
		{
			const tw::Float vn = J[c+1][i]/(tw::small_pos + J[0][i]);
			vel[c][i] = 1.5*vn - 0.5*u[c+1][i]; // extrapolate to n+1/2
			u[c+1][i] = vn;
		}
}


///////////////////////////////////////////
// BUNDLE TILERS - SLICING
//////////////////////////////////////////


inline void BundleTilerEM::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	owner->EM->LoadDataIntoImage<float>(&Fx);
}
inline void BundleTilerEM::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	Jx = 0.0f;
}
inline void BundleTilerEM::DepositSourceSlice(bool needsAtomic)
{
	if (needsAtomic)
		owner->sources->AddDataFromImageAtomic<float>(&Jx);
	else
		owner->sources->AddDataFromImage<float>(&Jx);
}

inline void BundleTilerPGC::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	owner->EM->LoadDataIntoImage<float>(&Fx);
	lasx.Resize(Element(0,7),low,high,ignorable);
	owner->laser->LoadDataIntoImage<float>(&lasx);
}
inline void BundleTilerPGC::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	Jx = 0.0f;
	chix.Resize(Element(0),low,high,ignorable);
	chix = 0.0f;
}
inline void BundleTilerPGC::DepositSourceSlice(bool needsAtomic)
{
	if (needsAtomic)
	{
		owner->sources->AddDataFromImageAtomic<float>(&Jx);
		owner->chi->AddDataFromImageAtomic<float>(&chix);
	}
	else
	{
		owner->sources->AddDataFromImage<float>(&Jx);
		owner->chi->AddDataFromImage<float>(&chix);
	}
}

inline void BundleTilerBohmian::LoadFieldSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	owner->qo_j4->LoadDataIntoImage<float>(&Jx);
}
inline void BundleTilerBohmian::InitSourceSlice(tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
}
inline void BundleTilerBohmian::DepositSourceSlice(bool needsAtomic)
{
}


///////////////////////////////////////////
// BUNDLE TILERS - TILING
//////////////////////////////////////////


inline void BundleTiler2D::LoadFTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int k=0;k<3;k++)
			for (tw::Int s=0;s<6;s++)
				F_tile[i][k][s] = Fx(ijk0[0]-1+i,0,ijk0[2]-1+k,s);
}
inline void BundleTiler2D::ResetJTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int k=0;k<5;k++)
			for (tw::Int s=0;s<4;s++)
				J_tile[i][k][s] = 0.0f;
}
inline void BundleTiler2D::StoreJTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int k=0;k<5;k++)
			for (tw::Int s=0;s<4;s++)
				Jx(ijk0[0]-2+i,0,ijk0[2]-2+k,s) += J_tile[i][k][s];
}
inline void BundleTiler3D::LoadFTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
				for (tw::Int s=0;s<6;s++)
					F_tile[i][j][k][s] = Fx(ijk0[0]-1+i,ijk0[1]-1+j,ijk0[2]-1+k,s);
}
inline void BundleTiler3D::ResetJTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
				for (tw::Int s=0;s<4;s++)
					J_tile[i][j][k][s] = 0.0f;
}
inline void BundleTiler3D::StoreJTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
				for (tw::Int s=0;s<4;s++)
					Jx(ijk0[0]-2+i,ijk0[1]-2+j,ijk0[2]-2+k,s) += J_tile[i][j][k][s];
}
inline void BundleTilerPGC::LoadLaserTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
				for (tw::Int s=0;s<8;s++)
					las_tile[i][j][k][s] = lasx(ijk0[0]-1+i,ijk0[1]-1+j,ijk0[2]-1+k,s);
}
inline void BundleTilerPGC::ResetChiTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
				chi_tile[i][j][k] = 0.0f;
}
inline void BundleTilerPGC::StoreChiTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
				chix(ijk0[0]-2+i,ijk0[1]-2+j,ijk0[2]-2+k,0) += chi_tile[i][j][k];
}
inline void BundleTilerBohmian::LoadTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
			{
				for (tw::Int s=0;s<4;s++)
					tile[i][j][k][s] = Jx(ijk0[0]-1+i,ijk0[1]-1+j,ijk0[2]-1+k,s);
			}
}
