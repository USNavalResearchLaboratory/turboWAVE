struct LaserSolver:Module
{
	tw::Float laserFreq;
	tw_polarization_type polarizationType;
	ComplexField a0,a1; // vector potential
	ComplexField chi; // defined by j = chi*a

	LaserPropagator *propagator;

	LaserSolver(const std::string& name,Simulation* sim);
	virtual ~LaserSolver();
	virtual void ExchangeResources();
	virtual void Initialize();
	virtual void Reset();

	virtual void VerifyInput();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void Update();
	tw::vec3 GetIonizationKick(const tw::Float& a2,const tw::Float& q0,const tw::Float& m0);
};

struct QSSolver:LaserSolver
{
	QSSolver(const std::string& name,Simulation* sim);
};

struct PGCSolver:LaserSolver
{
	Field F;

	PGCSolver(const std::string& name,Simulation* sim);
	virtual void ExchangeResources();
	virtual void Initialize();

	virtual void MoveWindow();
	virtual void AntiMoveWindow();

	virtual void Update();
	virtual void ComputeFinalFields();

	virtual void Report(Diagnostic&);
};
