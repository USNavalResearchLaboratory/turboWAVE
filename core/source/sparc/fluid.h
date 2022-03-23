namespace sparc
{
	enum laserModel { vacuum, isotropic };
	enum radiationModel { noRadiation, thin, thick };
	enum plasmaModel { neutral, quasineutral };
}

struct Fluid:Module
{
	tw::Float charge,mass,thermalMomentum,enCrossSection,initialIonizationFraction;
	Field state0,state1; // density,p1,p2,p3 (unlike SPARC state, p is not a momentum density)
	ScalarField fixed,gas;
	bool coulombCollisions;

	// temporaries used in Update
	Field vel; // gammaAvg,v1,v2,v3

	Ionizer *ionizer;

	Field *EM,*J4;
	Field *laser;
	ComplexField *chi;
	tw::Float *carrierFrequency;

	Fluid(const std::string& name,Simulation* sim);
	~Fluid();
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void VerifyInput();
	virtual void Initialize();
	virtual void Update();
	virtual void MoveWindow();
	virtual void AddDensity(tw::Float densityToAdd,tw::Int i,tw::Int j,tw::Int k);

	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void Report(Diagnostic&);
	virtual bool Test(tw::Int& id);
	void AdvectionTest();
	void ConservationTest();
};

struct EquilibriumGroup;

struct Chemical:Module
{
	EquilibriumGroup *group; // explictly typed super
	EOSComponent *eosData;
	Ionizer *ionizer;
	UniformProfile *background;
	sparc::material mat;
	tw::Int indexInState;

	Chemical(const std::string& name,Simulation* sim);
	virtual ~Chemical();
	void SetupIndexing();
	virtual void VerifyInput();

	bool LoadFluid(Field& hydro);
	void LoadInternalEnergy(Field& hydro,Field& eos);
};

struct EquilibriumGroup:Module
{
	std::vector<Chemical*> chemical; // explicitly typed submodule list
	EOSMixture *eosMixData;
	bool mobile;

	// The hydro set contains indices into the state vector for this group.
	// The mass density index corresponds to the first chemical in the group.
	sparc::hydro_set hidx;
	// The eos set contains indices into the eos vector for this group.
	sparc::eos_set eidx;
	// The material data is packed and processed for optimization
	sparc::material_set matset;
	// Following is used to limit motion by zeroing forces on this group
	tw::Float forceFilter;

	tw::Float DensitySum(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=hidx.first;s<hidx.first+hidx.num;s++)
			ans += f(cell,s);
		return ans;
	}
	tw::Float DensityWeightedSum(const Field& f,std::valarray<tw::Float>& qty,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=hidx.first;s<hidx.first+hidx.num;s++)
			ans += f(cell,s)*qty[s-hidx.first];
		return ans;
	}
	void LoadMassDensity(ScalarField& nm,const Field& f)
	{
		for (auto cell : EntireCellRange(*this))
			nm(cell) = DensityWeightedSum(f,matset.mass,cell);
	}
	void LoadMassDensityCv(ScalarField& nmcv,const Field& f)
	{
		for (auto cell : EntireCellRange(*this))
			nmcv(cell) = DensityWeightedSum(f,matset.cvm,cell);
	}
	tw::vec3 Velocity(const Field& f,const tw::cell& cell)
	{
		tw::Float nm = DensityWeightedSum(f,matset.mass,cell);
		return tw::vec3(f(cell,hidx.npx),f(cell,hidx.npy),f(cell,hidx.npz))/(tw::small_pos + nm);
	}
	void LoadVelocity(ScalarField& vel,const Field& f,tw::Int ax)
	{
		// assumes velocity components appear in order in state vector
		tw::Float nm;
		for (auto cell : EntireCellRange(*this))
		{
			nm = DensityWeightedSum(f,matset.mass,cell);
			vel(cell) = f(cell,hidx.npx+ax-1)/(tw::small_pos + nm);
		}
	}

	EquilibriumGroup(const std::string& name,Simulation* sim);
	virtual ~EquilibriumGroup(); // ASHER_MOD
	void SetupIndexing();
	virtual void VerifyInput();
	bool GenerateFluid(Field& hydro,Field& eos);
};

namespace sparc
{
struct HydroManager:Module
{
	std::vector<EquilibriumGroup*> group; // explicitly type submodule list
	std::vector<Reaction*> reaction;
	std::vector<Excitation*> excitation;
	std::vector<Collision*> collision;
	Field state0,state1,creationRate,destructionRate,eos0,eos1;

	ParabolicSolver *parabolicSolver;
	EllipticSolver *ellipticSolver;
	IsotropicPropagator *laserPropagator;
	tw::vec3 dipoleCenter;
	tw::Float laserFrequency; // derived from pulse list during initialization
	tw::Float backgroundDensity,backgroundTemperature;

	sparc::radiationModel radModel;
	sparc::laserModel lasModel;
	sparc::plasmaModel plasModel;
	bool electrostaticHeating;
	Chemical *electrons; // pointer to electron species, if present
	tw::Int ie; // index to electron density in state vector
	ScalarField scratch,scratch2,fluxMask;
	ScalarField rho0,rho,phi,nu_e,me_eff,radiativeLosses,radiationIntensity;
	ComplexField laserAmplitude,refractiveIndex;

	// items needed for step size control and user feedback
	tw::Float epsilonFactor;
	std::stringstream statusMessage;

	// The distinction between creation and destruction is used for adaptive timestep adjustment.
	// In the case of constant time step, there is no difference between create(x) and destroy(-x).
	// For signed quantities like momentum, the distinction is never used at all.
	void CreateMass(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ creationRate(cell,h.ni) += val; }
	void DestroyMass(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ destructionRate(cell,h.ni) += val; }
	void CreateMomentum(const tw::cell& cell,const tw::Int& ax,const tw::Float& val,const sparc::hydro_set& h)
		{ creationRate(cell,h.npx+ax-1) += val; }
	void DestroyMomentum(const tw::cell& cell,const tw::Int& ax,const tw::Float& val,const sparc::hydro_set& h)
		{ destructionRate(cell,h.npx+ax-1) += val; }
	void CreateTotalEnergy(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ creationRate(cell,h.u) += val; }
	void DestroyTotalEnergy(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ destructionRate(cell,h.u) += val; }
	void CreateVibrations(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ creationRate(cell,h.x) += val; }
	void DestroyVibrations(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ destructionRate(cell,h.x) += val; }
	void CreateTotalAndVibrational(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
	{
		creationRate(cell,h.u) += val;
		creationRate(cell,h.x) += val;
	}
	void DestroyTotalAndVibrational(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
	{
		destructionRate(cell,h.u) += val;
		destructionRate(cell,h.x) += val;
	}

	HydroManager(const std::string& name,Simulation* sim);
	virtual ~HydroManager();
	void SetupIndexing();
	virtual void Initialize();
	virtual void Reset();
	void LoadCollisionRate(Collision *coll,ScalarField& R);
	void ComputeElectronCollisionFrequency();
	void ComputeCollisionalSources();
	void ComputeRadiativeSources();
	void ComputeHydroSources();
	tw::vec3 ComputeForceOnBody(tw::Int i,tw::Int j,tw::Int k);
	void ComputeSources();
	tw::Float EstimateTimeStep();

	void HydroAdvance(const tw::grid::axis& axis,tw::Float dt);
	void LaserAdvance(tw::Float dt);
	void ChemAdvance(tw::Float dt);
	void DiffusionAdvance(tw::Float dt);
	void FieldAdvance(tw::Float dt);
	void EOSAdvance(tw::Float dt);
	void FirstOrderAdvance(tw::Float dt,bool computeSources);
	virtual void Update();

	virtual void VerifyInput();
	virtual bool ReadQuasitoolBlock(const tw::input::Preamble& preamble,std::stringstream& inputString);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void Report(Diagnostic&);
	virtual void StatusMessage(std::ostream *theStream);
};
}
