struct Species;
struct Kinetics;

// Particle struct is defined in discreteSpace.h

struct TransferParticle
{
	// To avoid inconsistencies arising from FP comparisons on different nodes
	// the destination information is computed on the source domain and packaged with the particle
	// The position is kept as a double precision global coordinate until final call to AddParticle
	tw::Int dst[4];
	// dst[0] is rank of starting domain upon construction; gets set to destination domain later.
	// dst[1..3] are +-1, giving direction of movement; zero if no movement.
	tw::vec4 x,p; // can use x[0] or p[0] to pack extra info
	float number,aux1,aux2;
};

struct ParticleRef
{
	// Used to create sorting map within a thread for subsets of particle lists
	tw::Int idx,cell;
	ParticleRef() noexcept
	{
		idx = 0;
		cell = 0;
	}
	ParticleRef(tw::Int list_index,const Particle& par) noexcept
	{
		idx = list_index;
		cell = par.q.cell;
	}
	friend bool operator < (const ParticleRef& r1,const ParticleRef& r2)
	{
		return r1.cell < r2.cell;
	}
};

struct ParticleBundle
{
	// The purpose of the particle bundle is to allow for efficient vectorization of the pusher
	// Several inlined functions are needed due to OpenMP inability to work with member variables
	// Evidently this is related to uninstantiated variables being tagged as aligned
	// The workaround is to pass members back in as arguments
	tw::Int num,cell0,ijk0[3];
	static const tw::Int N = tw::max_bundle_size;
	static const tw::Int AB = tw::vec_align_bytes;
	// arrays storing particle state
	alignas(AB) tw::Int cell[N];
	alignas(AB) float x[3][N];
	alignas(AB) tw::Float p[3][N];
	alignas(AB) float number[N];
	// temporary arrays used during calculation
	alignas(AB) float domainMask[N];
	alignas(AB) float cellMask[N];
	alignas(AB) tw::Int ijk[3][N];
	alignas(AB) tw::Float t[3][N];
	alignas(AB) tw::Float s[3][N];
	alignas(AB) tw::Float vel[3][N];
	alignas(AB) float F[6][N];
	alignas(AB) float J[4][N];
	alignas(AB) float w0[3][3][N];
	alignas(AB) float w1[3][3][N];
	alignas(AB) float l0[3][3][N];

	std::valarray<Particle*> refs;

	ParticleBundle();
	void PadBundle();
	void Reset();
	bool Complete(const Particle& par);
	void Append(Particle& par);
	void CopyBack(Species *owner);
	void impulse(tw::Float p[3][N],float F[6][N],const tw::Float& q0, const tw::Float& dth);
	void rotation1(tw::Float t[3][N],tw::Float p[3][N],float F[6][N],const tw::Float& q0,const tw::Float& m0,const tw::Float& dth);
	void rotation2(tw::Float s[3][N],tw::Float t[3][N]);
	void rotation3(tw::Float s[3][N],tw::Float t[3][N],tw::Float vel[3][N],tw::Float p[3][N]);
	void velocity(tw::Float vel[3][N],tw::Float p[3][N],const tw::Float& m0);
	void bohm_velocity(tw::Float vel[3][N],tw::Float p[3][N],float J[4][N]);
	void translate(float x[3][N],tw::Float vel[3][N],const tw::Float k[3],const tw::Float& dt);
	void set_cell_mask(float cellMask[N],tw::Int cell0,tw::Int cell[N]);
	void get_cell_displ(tw::Int dc[3],tw::Int n);
	void load_j4(float J[4][N],const float number[N],const tw::Float vel[3][N],const tw::Float k[3],const tw::Float& q0);
	void ZeroArray(float q[][N],tw::Int s1,tw::Int s2);
};

struct ParticleBundleBohmian : ParticleBundle
{
	Slice<float> Jx;
	float tile[3][3][3][4];
	void LoadFieldSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void InitSourceSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void DepositSourceSlice(Species *owner,bool needsAtomic);
	void GatherJ4(float J[4][N],const float w0[3][3][N]);
	void Push(Species *owner);
	void LoadTile();
};

struct ParticleBundleEM : ParticleBundle
{
	Slice<float> Fx,Jx;
	void LoadFieldSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void InitSourceSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void DepositSourceSlice(Species *owner,bool needsAtomic);
};

struct ParticleBundle2D : ParticleBundleEM
{
	float tile[3][3][6];
	float xtile[5][5][4];

	void GatherF(float F[6][N],const float w0[3][3][N],const float l0[3][3][N]);
	void ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti);
	void Push(Species *owner);

	void LoadTile();
	void ResetXTile();
	void StoreXTile();
};

struct ParticleBundle3D : ParticleBundleEM
{
	float tile[3][3][3][6];
	float xtile[5][5][5][4];

	void GatherF(float F[6][N],const float w0[3][3][N],const float l0[3][3][N]);
	void ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti);
	void Push(Species *owner);

	void LoadTile();
	void ResetXTile();
	void StoreXTile();
};

struct ParticleBundleElectrostatic : ParticleBundle3D
{
	void LoadFieldSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void GatherF(float F[6][N],const float w0[3][3][N],const float l0[3][3][N]);
	void Push(Species *owner);
	void LoadTile();
};

struct ParticleBundlePGC : ParticleBundle3D
{
	alignas(AB) float las[8][N];
	alignas(AB) float chi[N];
	alignas(AB) tw::Float avgMass[N];
	Slice<float> lasx,chix;
	float las_tile[3][3][3][8];
	float chi_tile[5][5][5];

	void avg_mass_1(tw::Float avgMass[N],tw::Float vel[3][N],tw::Float p[3][N],float F[6][N],float las[8][N],const tw::Float& q0,const tw::Float& m0,const tw::Float& dth);
	void avg_mass_2(tw::Float avgMass[N],tw::Float p[3][N],float las[8][N],const tw::Float& q0,const tw::Float& m0,const tw::Float& dth);
	void impulse(tw::Float p[3][N],float F[6][N],float las[8][N],tw::Float avgMass[N],const tw::Float& q0,const tw::Float& dth);
	void rotation1(tw::Float t[3][N],float F[6][N],tw::Float avgMass[N],const tw::Float& q0,const tw::Float& dth);
	void velocity(tw::Float vel[3][N],tw::Float p[3][N],tw::Float avgMass[N]);
	void load_chi(float chi[N],float number[N],tw::Float avgMass[N],const tw::Float& q0);

	void LoadFieldSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void InitSourceSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4]);
	void DepositSourceSlice(Species *owner,bool needsAtomic);

	void GatherLaser(float las[8][N],const float w0[3][3][N]);
	void ScatterChi(const float chi[N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N]);
	void Push(Species *owner);

	void LoadTile();
	void ResetXTile();
	void StoreXTile();
};

struct LoadingData
{
	tw::Int i,j,k;
	tw::Float densToAdd,densNow,particleDensity;
	tw::vec3 thermalMomentum,driftMomentum;
	std::valarray<tw::vec3> subGrid;
	tw::Int pointsInSubGrid;
	bool neutralize;

	LoadingData(const tw::vec3& distributionInCell);
};

struct Species:Module
{
	tw::Float restMass,charge;
	tw::vec3 emissionTemp;
	tw::vec3 distributionInCell;
	tw::Float targetDensity;
	tw::Float minimumDensity;
	tw::Float accelerationTime,accelerationImpulse,accelerationForceNow;
	tw::bc::par bc0[4],bc1[4];
	bool mobile,radiationDamping;
	tw::Float meanFreePath;
	tw::Int count,sortPeriod;

	std::vector<Particle> particle;
	std::vector<TransferParticle> transfer;

	Ionizer* ionizer;
	Field* EM; // Ex,Ey,Ez,Bx,By,Bz
	Field* sources; // rho,Jx,Jy,Jz
	Field* laser; // F0x,F0y,F0z,F1x,F1y,F1z,aa0,aa1
	ComplexField* chi;
	tw::Float* carrierFrequency;
	tw_polarization_type* polarizationType;
	ScalarField* rho00;
	Vec3Field* ESField;

	Field *qo_j4; // 4-current from quantum optics modules

	Species(const std::string& name,Simulation* sim);
	virtual ~Species();
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void VerifyInput();
	virtual void Initialize();
	void AddParticle(const tw::vec3& p,const Primitive& q,const float& number);
	void AddParticle(const TransferParticle& newParticle);
	void AddTransferParticle(const Particle& src);
	void CleanParticleList();

	void ApplyGlobalBoundaryConditions();

	void ComputeTransferParticleDestinations();
	void PrepareTransfer(std::vector<TransferParticle>& accumulator,std::vector<tw::Int>& tally,tw::Int axis,tw::Int displ);
	void FinishTransfer(TransferParticle* inBuffer,tw::Int number,tw::Int axis,tw::Int displ);
	void CollectTransfers();

	virtual void MoveWindow();
	void BeginMoveWindow();
	void FinishMoveWindow();

	void GenerateParticles(bool init);
	tw::Float AddDensity(const LoadingData& theData);
	tw::Float AddDensityRandom(const LoadingData& theData);
	void DepositInitialCharge(const tw::vec3& pos,tw::Float macroCharge);
	void CalculateDensity(ScalarField& dens);

	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
	virtual void Report(Diagnostic&);
	virtual void WarningMessage(std::ostream *theStream);

	void GetSubarrayBounds(std::vector<ParticleRef>& sorted,tw::Int low[4],tw::Int high[4],tw::Int layers);
	void SpreadTasks(std::vector<tw::Int>& task_map);
	void BunchTasks(std::vector<tw::Int>& task_map);
	void DispatchPush();
	template <class BundleType>
	void Push();
	template <class BundleType>
	void PushSlice(tw::Int tasks,tw::Int tid,tw::Int bounds_data[][8]);
};

struct Kinetics:Module
{
	std::vector<Species*> species; // explicitly typed copy of submodule list

	ScalarField rho00;
	Field* sources;
	ComplexField* chi;
	ScalarField* ESRho;

	Kinetics(const std::string& name,Simulation* sim);
	virtual void Initialize();
	virtual bool ValidSubmodule(Module* sub);
	virtual void ExchangeResources();
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void Update();
	virtual void MoveWindow();
	void Ionize();

	void TransferParticles();
	tw::Float KineticEnergy(const Region& theRgn);
	virtual void Report(Diagnostic&);

	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};


//////////////////////////
//                      //
// INLINE PUSH ROUTINES //
//                      //
//////////////////////////


inline ParticleBundle::ParticleBundle()
{
	refs.resize(N);
	num = 0;
}
inline void ParticleBundle::PadBundle()
{
	for (int ax=0;ax<3;ax++)
		for (int i=num;i<N;i++)
		{
			x[ax][i] = 0.0;
			p[ax][i] = 0.0;
		}
	for (int i=num;i<N;i++)
	{
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
	x[0][num] = par.q.x[0];
	x[1][num] = par.q.x[1];
	x[2][num] = par.q.x[2];
	p[0][num] = par.p[0];
	p[1][num] = par.p[1];
	p[2][num] = par.p[2];
	number[num] = par.number;
	num++;
}
inline void ParticleBundle::CopyBack(Species *owner)
{
	for (int i=0;i<num;i++)
	{
		refs[i]->q.cell = cell[i];
		refs[i]->q.x[0] = x[0][i];
		refs[i]->q.x[1] = x[1][i];
		refs[i]->q.x[2] = x[2][i];
		refs[i]->p[0] = p[0][i];
		refs[i]->p[1] = p[1][i];
		refs[i]->p[2] = p[2][i];
		// If particle left MPI domain, put a copy on transfer list, and mark for disposal
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
inline void ParticleBundle::impulse(tw::Float p[3][N],float F[6][N],const tw::Float& q0, const tw::Float& dth)
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(p,F:AB)
		for (int i=0;i<N;i++)
			p[c][i] += q0*F[c][i]*dth;
}
inline void ParticleBundle::rotation1(tw::Float t[3][N],tw::Float p[3][N],float F[6][N],const tw::Float& q0,const tw::Float& m0,const tw::Float& dth)
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(t,p,F:AB)
		for (int i=0;i<N;i++)
			t[c][i] = q0*F[c+3][i]*dth/sqrt(m0*m0 + p[0][i]*p[0][i] + p[1][i]*p[1][i] + p[2][i]*p[2][i]);
}
inline void ParticleBundle::rotation2(tw::Float s[3][N],tw::Float t[3][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(s,t:AB)
		for (int i=0;i<N;i++)
			s[c][i] = 2.0*t[c][i]/(1.0 + t[0][i]*t[0][i] + t[1][i]*t[1][i] + t[2][i]*t[2][i]);
}
inline void ParticleBundle::rotation3(tw::Float s[3][N],tw::Float t[3][N],tw::Float vel[3][N],tw::Float p[3][N])
{
	// Use vel as temporary while adding rotational impulse, cross(p + cross(p,t),s)
	#pragma omp simd aligned(s,t,vel,p:AB)
	for (int i=0;i<N;i++)
	{
		vel[0][i] = p[0][i] + p[1][i]*t[2][i] - p[2][i]*t[1][i];
		vel[1][i] = p[1][i] + p[2][i]*t[0][i] - p[0][i]*t[2][i];
		vel[2][i] = p[2][i] + p[0][i]*t[1][i] - p[1][i]*t[0][i];
		p[0][i] += vel[1][i]*s[2][i] - vel[2][i]*s[1][i];
		p[1][i] += vel[2][i]*s[0][i] - vel[0][i]*s[2][i];
		p[2][i] += vel[0][i]*s[1][i] - vel[1][i]*s[0][i];
	}
}
inline void ParticleBundle::velocity(tw::Float vel[3][N],tw::Float p[3][N],const tw::Float& m0)
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,p:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = p[c][i]/sqrt(m0*m0 + p[0][i]*p[0][i] + p[1][i]*p[1][i] + p[2][i]*p[2][i]);
}
inline void ParticleBundle::bohm_velocity(tw::Float vel[3][N],tw::Float p[3][N],float J[4][N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,p,J:AB)
		for (int i=0;i<N;i++)
		{
			const tw::Float vn = J[c+1][i]/(tw::small_pos + J[0][i]);
			vel[c][i] = 1.5*vn - 0.5*p[c][i]; // extrapolate to n+1/2
			p[c][i] = vn;
		}
}
inline void ParticleBundle::translate(float x[3][N],tw::Float vel[3][N],const tw::Float k[3],const tw::Float& dt)
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(x,vel:AB)
		for (int i=0;i<N;i++)
			x[c][i] += vel[c][i]*dt*k[c];
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
inline void ParticleBundle::load_j4(float J[4][N],const float number[N],const tw::Float vel[3][N],const tw::Float k[3],const tw::Float& q0)
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
inline void ParticleBundleEM::LoadFieldSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	owner->EM->LoadDataIntoImage<float>(&Fx);
}
inline void ParticleBundleEM::InitSourceSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	Jx = 0.0f;
}
inline void ParticleBundleEM::DepositSourceSlice(Species *owner,bool needsAtomic)
{
	if (needsAtomic)
		owner->sources->AddDataFromImageAtomic<float>(&Jx);
	else
		owner->sources->AddDataFromImage<float>(&Jx);
}
inline void ParticleBundleElectrostatic::LoadFieldSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,2),low,high,ignorable);
	owner->ESField->LoadDataIntoImage<float>(&Fx);
}
inline void ParticleBundleElectrostatic::LoadTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
				for (tw::Int s=0;s<3;s++)
					tile[i][j][k][s] = Fx(ijk0[0]-1+i,ijk0[1]-1+j,ijk0[2]-1+k,s);
}
inline void ParticleBundlePGC::LoadFieldSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Fx.Resize(Element(0,5),low,high,ignorable);
	owner->EM->LoadDataIntoImage<float>(&Fx);
	lasx.Resize(Element(0,7),low,high,ignorable);
	owner->laser->LoadDataIntoImage<float>(&lasx);
}
inline void ParticleBundlePGC::InitSourceSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	Jx = 0.0f;
	chix.Resize(Element(0),low,high,ignorable);
	chix = 0.0f;
}
inline void ParticleBundlePGC::DepositSourceSlice(Species *owner,bool needsAtomic)
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
inline void ParticleBundle2D::LoadTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int k=0;k<3;k++)
			for (tw::Int s=0;s<6;s++)
				tile[i][k][s] = Fx(ijk0[0]-1+i,0,ijk0[2]-1+k,s);
}
inline void ParticleBundle2D::ResetXTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int k=0;k<5;k++)
			for (tw::Int s=0;s<4;s++)
				xtile[i][k][s] = 0.0f;
}
inline void ParticleBundle2D::StoreXTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int k=0;k<5;k++)
			for (tw::Int s=0;s<4;s++)
				Jx(ijk0[0]-2+i,0,ijk0[2]-2+k,s) += xtile[i][k][s];
}
inline void ParticleBundle3D::LoadTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
				for (tw::Int s=0;s<6;s++)
					tile[i][j][k][s] = Fx(ijk0[0]-1+i,ijk0[1]-1+j,ijk0[2]-1+k,s);
}
inline void ParticleBundle3D::ResetXTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
				for (tw::Int s=0;s<4;s++)
					xtile[i][j][k][s] = 0.0f;
}
inline void ParticleBundle3D::StoreXTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
				for (tw::Int s=0;s<4;s++)
					Jx(ijk0[0]-2+i,ijk0[1]-2+j,ijk0[2]-2+k,s) += xtile[i][j][k][s];
}

inline void ParticleBundlePGC::avg_mass_1(tw::Float avgMass[N],tw::Float vel[3][N],tw::Float p[3][N],float F[6][N],float las[8][N],const tw::Float& q0,const tw::Float& m0,const tw::Float& dth)
{
	alignas(AB) tw::Float m1[N];
	alignas(AB) tw::Float dA[N];
	// estimate avgMass at level n using solution of <m>^2 = m1<m> - dA/4
	#pragma omp simd aligned(m1,p,las:AB)
	for (int i=0;i<N;i++)
		m1[i] = sqrt(m0*m0 + p[0][i]*p[0][i] + p[1][i]*p[1][i] + p[2][i]*p[2][i] + 0.5*q0*q0*las[7][i]);
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,p,m1:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = p[c][i]/m1[i];
	#pragma omp simd aligned(m1,F,vel,dA,las,avgMass:AB)
	for (int i=0;i<N;i++)
	{
		m1[i] += dth*q0*(F[0][i]*vel[0][i] + F[1][i]*vel[1][i] + F[2][i]*vel[2][i]);
		dA[i] = dth*q0*q0*(las[0][i]*vel[0][i] + las[1][i]*vel[1][i] + las[2][i]*vel[2][i]);
		avgMass[i] = 0.5*m1[i] + 0.5*sqrt(m1[i]*m1[i] - dA[i]);
	}
}
inline void ParticleBundlePGC::avg_mass_2(tw::Float avgMass[N],tw::Float p[3][N],float las[8][N],const tw::Float& q0,const tw::Float& m0,const tw::Float& dth)
{
	// estimate avgMass at level n+1/2 using <m>^2 = m1^2 + (m1/<m>)dA/2
	// (solve cubic and expand in powers of dA/m1^2)
	alignas(AB) tw::Float m1[N];
	alignas(AB) tw::Float dA[N];
	#pragma omp simd aligned(m1,dA,p,las,avgMass:AB)
	for (int i=0;i<N;i++)
	{
		m1[i] = sqrt(m0*m0 + p[0][i]*p[0][i] + p[1][i]*p[1][i] + p[2][i]*p[2][i] + 0.5*q0*q0*las[7][i]);
		dA[i] = dth*q0*q0*(las[3][i]*p[0][i] + las[4][i]*p[1][i] + las[5][i]*p[2][i])/m1[i];
		avgMass[i] = m1[i] + 0.25*dA[i]/m1[i] - 0.09375*dA[i]*dA[i]/(m1[i]*m1[i]*m1[i]);
	}
}
inline void ParticleBundlePGC::impulse(tw::Float p[3][N],float F[6][N],float las[8][N],tw::Float avgMass[N],const tw::Float& q0,const tw::Float& dth)
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(p,F,las,avgMass:AB)
		for (int i=0;i<N;i++)
			p[c][i] += q0*F[c][i]*dth - 0.25*q0*q0*las[c][i]*dth/avgMass[i];
}
inline void ParticleBundlePGC::rotation1(tw::Float t[3][N],float F[6][N],tw::Float avgMass[N],const tw::Float& q0,const tw::Float& dth)
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(t,F,avgMass:AB)
		for (int i=0;i<N;i++)
			t[c][i] = q0*F[c+3][i]*dth/avgMass[i];
}
inline void ParticleBundlePGC::velocity(tw::Float vel[3][N],tw::Float p[3][N],tw::Float avgMass[N])
{
	for (int c=0;c<3;c++)
		#pragma omp simd aligned(vel,p,avgMass:AB)
		for (int i=0;i<N;i++)
			vel[c][i] = p[c][i]/avgMass[i];
}
inline void ParticleBundlePGC::load_chi(float chi[N],float number[N],tw::Float avgMass[N],const tw::Float& q0)
{
	#pragma omp simd aligned(chi,number,avgMass:AB)
	for (int i=0;i<N;i++)
		chi[i] = -q0*q0*number[i]/avgMass[i];
}

inline void ParticleBundlePGC::LoadTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
			{
				for (tw::Int s=0;s<6;s++)
					tile[i][j][k][s] = Fx(ijk0[0]-1+i,ijk0[1]-1+j,ijk0[2]-1+k,s);
				for (tw::Int s=0;s<8;s++)
					las_tile[i][j][k][s] = lasx(ijk0[0]-1+i,ijk0[1]-1+j,ijk0[2]-1+k,s);
			}
}
inline void ParticleBundlePGC::ResetXTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
			{
				for (tw::Int s=0;s<4;s++)
					xtile[i][j][k][s] = 0.0f;
				chi_tile[i][j][k] = 0.0f;
			}
}
inline void ParticleBundlePGC::StoreXTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
			{
				for (tw::Int s=0;s<4;s++)
					Jx(ijk0[0]-2+i,ijk0[1]-2+j,ijk0[2]-2+k,s) += xtile[i][j][k][s];
				chix(ijk0[0]-2+i,ijk0[1]-2+j,ijk0[2]-2+k,0) += chi_tile[i][j][k];
			}
}

inline void ParticleBundleBohmian::LoadFieldSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
	Jx.Resize(Element(0,3),low,high,ignorable);
	owner->qo_j4->LoadDataIntoImage<float>(&Jx);
}
inline void ParticleBundleBohmian::InitSourceSlice(Species *owner,tw::Int low[4],tw::Int high[4],tw::Int ignorable[4])
{
}
inline void ParticleBundleBohmian::DepositSourceSlice(Species *owner,bool needsAtomic)
{
}
inline void ParticleBundleBohmian::LoadTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
			{
				for (tw::Int s=0;s<4;s++)
					tile[i][j][k][s] = Jx(ijk0[0]-1+i,ijk0[1]-1+j,ijk0[2]-1+k,s);
			}
}
