struct Species;
struct Kinetics;

// Particle struct is defined in discreteSpace.h
// Heavy lifting is done by Bundle classes in pusher.h

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

struct LoadingData
{
	tw::cell cell;
	tw::Float C0,C1,C2,densToAdd,densNow,particleDensity;
	tw::vec3 thermalMomentum,driftMomentum;
	std::valarray<tw::vec3> subGrid;
	tw::Int pointsInSubGrid;
	bool neutralize;

	LoadingData(const Simulation& sim,const tw::vec3& distributionInCell,const tw::cell& c);
	tw::Float GeometryFactor(tw::Float r,tw::Float r0) const
	{
		return fabs(C0 + SafeDiv(C1*r,r0) + sqr(SafeDiv(sqrt(C2)*r,r0)));
	}
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
	std::vector<Species*> species; // explicitly typed submodules

	ScalarField rho00;
	Field* sources;
	ComplexField* chi;
	ScalarField* ESRho;

	Kinetics(const std::string& name,Simulation* sim);
	virtual void Initialize();
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
