struct Species;
struct Kinetics;

// Particle struct is defined in discreteSpace.h
// Heavy lifting is done by classes in particles_*.h

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

	Mover* mover;
	Ionizer* ionizer;

	Field* EM; // Ex,Ey,Ez,Bx,By,Bz
	Field* sources; // rho,Jx,Jy,Jz
	Field* laser; // F0x,F0y,F0z,F1x,F1y,F1z,aa0,aa1
	ComplexField* chi;
	tw::Float* carrierFrequency;
	tw_polarization_type* polarizationType;
	ScalarField* rho00;

	Field *qo_j4; // 4-current from quantum optics modules

	Species(const std::string& name,Simulation* sim);
	virtual ~Species();
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void VerifyInput();
	virtual void Initialize();
	void AddParticle(const float& number,const Primitive& q,const tw::vec4& p,const tw::vec4& s);
	void AddParticle(const TransferParticle& newParticle);
	void CleanParticleList();

	void ApplyGlobalBoundaryConditions();

	void ComputeTransferParticleDestinations();
	void PrepareTransfer(std::vector<TransferParticle>& accumulator,std::vector<tw::Int>& tally,tw::Int axis,tw::Int displ);
	void FinishTransfer(TransferParticle* inBuffer,tw::Int number,tw::Int axis,tw::Int displ);
	void CollectTransfers();

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

	virtual bool Test(tw::Int& id);
	void ReflectionTest();
	void MoveWindowTest();
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
