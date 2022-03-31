struct Electrostatic:FieldSolver
{
	Field F,J4;
	ScalarField phi,source;

	Electrostatic(const std::string& name,Simulation* sim);
	virtual void Initialize();
	virtual void ExchangeResources();
	virtual void Reset();
	virtual void Update();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void SetupInitialPotential();
	virtual void ComputeFinalFields();

	virtual void Report(Diagnostic&);
};
