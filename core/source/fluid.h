namespace sparc
{
	enum collisionType { hard_sphere, coulomb, metallic };
	enum laserModel { vacuum, isotropic };
	enum radiationModel { noRadiation, thin, thick };
	enum plasmaModel { neutral, quasineutral };
}

inline tw::Float clippingFunction(tw::Float x)
{
	if (x>0.0 && x<2.0)
		return 1.0 + 2.0*cub(0.5*x-1) + sqr(0.5*x-1)*sqr(0.5*x-1);
	if (x>=2.0)
		return 1.0;
	return 0.0;
}

struct Fluid:Module
{
	tw::Float charge,mass,thermalMomentum,enCrossSection,initialIonizationFraction;
	Field state0,state1; // density,p1,p2,p3 (unlike SPARC state, p is not a momentum density)
	ScalarField fixed,gas;
	bool coulombCollisions;

	// temporaries used in Update
	Field vel; // gammaAvg,v1,v2,v3

	IonizationData ionization;

	Field *EM,*J4;
	Field *laser;
	ComplexField *chi;
	tw::Float *carrierFrequency;

	Fluid(Grid* theGrid);
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void Initialize();
	virtual void Update();
	virtual void MoveWindow();
	virtual void AddDensity(tw::Float densityToAdd,tw::Int i,tw::Int j,tw::Int k);

	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);

	virtual void BoxDiagnosticHeader(GridDataDescriptor*);
	virtual void BoxDiagnose(GridDataDescriptor*);
	virtual void PointDiagnosticHeader(std::ofstream& outFile);
	virtual void PointDiagnose(std::ofstream& outFile,const weights_3D& w);
};

struct Chemical;
struct Chemistry;

struct SubReaction
{
	std::vector<sparc::hydro_set> reactants,products;
	tw::Float heat,vheat;

	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct Reaction
{
	std::vector<SubReaction*> sub;
	tw::Int catalyst,numBodies;
	tw::Float T0,T1; // temperature range
	tw::Float c1,c2,c3; // arrhenius form
	tw::Float b[9]; // janev coefficients

	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct Excitation
{
	sparc::hydro_set exciter,excitee;
	tw::Float T0,T1; // temperature range (not used)
	tw::Float c1,c2,c3; // arrhenius form
	tw::Float b[9]; // janev coefficients
	tw::Float level;

	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct Collision
{
	sparc::hydro_set body1,body2;
	sparc::collisionType type;
	tw::Float crossSection;
	tw::Float ks,T_ref,n_ref;

	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct EquilibriumGroup:Module
{
	std::vector<Chemical*> chemical; // explicitly typed submodule list

	// Following contains indices in state array corresponding to density of each chemical in the group
	Element e;
	// The hydro set contains indices into the state vector for this group.
	// The mass density index corresponds to the first chemical in the group.
	sparc::hydro_set hidx;
	// The eos set contains indices into the eos vector for this group.
	sparc::eos_set eidx;
	// Following is used to limit motion by zeroing forces on this group
	tw::Float forceFilter;

	EOSMixture *eosMixData;
	// duplicate data from chemical list for efficiency
	std::valarray<tw::Float> mass,charge,cvm,excitationEnergy,thermo_cond_cvm,k_visc_m;

	tw::Float DensitySum(const Field& f,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=e.low;s<=e.high;s++)
			ans += f(cell,s);
		return ans;
	}
	tw::Float DensityWeightedSum(const Field& f,std::valarray<tw::Float>& qty,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=e.low;s<=e.high;s++)
			ans += f(cell,s)*qty[s-e.low];
		return ans;
	}
	tw::Float ConditionalDensitySum(const Field& f,std::valarray<tw::Float>& qty,const CellIterator& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=e.low;s<=e.high;s++)
			ans += qty[s-e.low] > 0.0 ? f(cell,s) : 0.0;
		return ans;
	}
	void LoadMassDensity(ScalarField& nm,const Field& f)
	{
		for (CellIterator cell(*this,true);cell<cell.end();++cell)
			nm(cell) = DensityWeightedSum(f,mass,cell);
	}
	void LoadMassDensityCv(ScalarField& nmcv,const Field& f)
	{
		for (CellIterator cell(*this,true);cell<cell.end();++cell)
			nmcv(cell) = DensityWeightedSum(f,cvm,cell);
	}
	tw::vec3 Velocity(const Field& f,const CellIterator& cell)
	{
		tw::Float nm = DensityWeightedSum(f,mass,cell);
		return tw::vec3(f(cell,npx),f(cell,npy),f(cell,npz))/(tw::small_pos + nm);
	}
	void LoadVelocity(ScalarField& vel,const Field& f,tw::Int ax)
	{
		// assumes velocity components appear in order in state vector
		tw::Float nm;
		for (CellIterator cell(*this,true);cell<cell.end();++cell)
		{
			nm = DensityWeightedSum(f,mass,cell);
			vel(cell) = f(cell,npx+ax-1)/(tw::small_pos + nm);
		}
	}

	EquilibriumGroup(Grid* theGrid);
	virtual ~EquilibriumGroup(); // ASHER_MOD
	virtual void Initialize();

	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct Chemical:Module
{
	EquilibriumGroup *group; // explictly typed super
	tw::Int indexInGroup,indexInState;
	tw::Float mass,charge,cvm; // cvm = cv*mass/kB = 3/2 for point particles
	IonizationData ionization;
	tw::Float excitationEnergy,thermometricConductivity,kinematicViscosity;

	EOSDataTool *eosData;

	// electron motion in lattice
	tw::Float effectiveMass,transitionDensity;
	tw::Complex permittivity;

	Chemical(Grid* theGrid);
	virtual ~Chemical();
	virtual void Initialize();

	bool GenerateFluid(Field& f,bool init);
	void DefaultEOS();

	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct Chemistry:Module
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

	sparc::radiationModel radModel;
	sparc::laserModel lasModel;
	sparc::plasmaModel plasModel;
	Chemical *electrons; // pointer to electron species, if present
	tw::Int ie; // index to electron density in state vector
	ScalarField scratch,scratch2,fluxMask;
	ScalarField rho0,rho,phi,nu_e,me_eff,radiativeLosses,radiationIntensity;
	tw::Float phi_lbc,phi_rbc;
	ComplexField laserAmplitude,refractiveIndex;
	std::valarray<tw::Float> gres;

	// items needed for step size control and relaxation
	tw::Float epsilonFactor;
	tw::Float tolerance,overrelaxation;
	tw::Int maxIterations;
	std::stringstream statusMessage;

	// The distinction between creation and destruction is used for adaptive timestep adjustment.
	// In the case of constant time step, there is no difference between create(x) and destroy(-x).
	// For signed quantities like momentum, the distinction is never used at all.
	void CreateMass(const CellIterator& cell,const tw::Int& c,const tw::Float& val,EquilibriumGroup *g) { creationRate(cell,g->e.low+c) += val; }
	void DestroyMass(const CellIterator& cell,const tw::Int& c,const tw::Float& val,EquilibriumGroup *g) { destructionRate(cell,g->e.low+c) += val; }
	void CreateMomentum(const CellIterator& cell,const tw::Int& ax,const tw::Float& val,EquilibriumGroup *g) { creationRate(cell,g->npx+ax-1) += val * g->forceFilter; }
	void DestroyMomentum(const CellIterator& cell,const tw::Int& ax,const tw::Float& val,EquilibriumGroup *g) { destructionRate(cell,g->npx+ax-1) += val * g->forceFilter; }
	void CreateTotalEnergy(const CellIterator& cell,const tw::Float& val,EquilibriumGroup *g) { creationRate(cell,g->U) += val; }
	void DestroyTotalEnergy(const CellIterator& cell,const tw::Float& val,EquilibriumGroup *g) { destructionRate(cell,g->U) += val; }
	void CreateVibrations(const CellIterator& cell,const tw::Float& val,EquilibriumGroup *g) { creationRate(cell,g->Xi) += val; }
	void DestroyVibrations(const CellIterator& cell,const tw::Float& val,EquilibriumGroup *g) { destructionRate(cell,g->Xi) += val; }
	void CreateTotalAndVibrational(const CellIterator& cell,const tw::Float& val,EquilibriumGroup *g)
	{
		creationRate(cell,g->U) += val;
		creationRate(cell,g->Xi) += val;
	}
	void DestroyTotalAndVibrational(const CellIterator& cell,const tw::Float& val,EquilibriumGroup *g)
	{
		destructionRate(cell,g->U) += val;
		destructionRate(cell,g->Xi) += val;
	}

	Chemistry(Grid* theGrid);
	virtual ~Chemistry();
	virtual void Initialize();
	virtual void Reset();

	tw::Float CoulombCrossSection(tw::Float q1,tw::Float q2,tw::Float m12,tw::Float v12,tw::Float N1,tw::Float N2,tw::Float T1,tw::Float T2) const;
	tw::Float ElectronPhononRateCoeff(tw::Float Ti,tw::Float EFermi,tw::Float ks,tw::Float nref) const;
	void ComputeElectronCollisionFrequency();
	void ComputeCollisionalSources();
	void ComputeRadiativeSources();
	void ComputeHydroSources();
	tw::vec3 ComputeForceOnBody(tw::Int i,tw::Int j,tw::Int k);
	void ComputeSources();
	tw::Float EstimateTimeStep();

	void HydroAdvance(const axisSpec& axis,tw::Float dt);
	void LaserAdvance(tw::Float dt);
	void ChemAdvance(tw::Float dt);
	void DiffusionAdvance(tw::Float dt);
	void FieldAdvance(tw::Float dt);
	void ApplyEOS(Field& hydro,Field& eos);
	void FirstOrderAdvance(tw::Float dt,bool computeSources);
	virtual void Update();

	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
	void ParseReaction(std::stringstream& inputString);
	void ParseExcitation(std::stringstream& inputString);
	void ParseCollision(std::stringstream& inputString);

	virtual void EnergyHeadings(std::ofstream& outFile);
	virtual void EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn);
	virtual void BoxDiagnosticHeader(GridDataDescriptor*);
	virtual void BoxDiagnose(GridDataDescriptor*);
	virtual void CustomDiagnose();
	virtual void StatusMessage(std::ostream *theStream);
};
