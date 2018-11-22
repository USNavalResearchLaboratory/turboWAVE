struct Electrostatic:FieldSolver
{
	Field sources;
	ScalarField phi,source;
	Vec3Field Ef;
	std::valarray<tw::Float> lbc,rbc; // 2-ghost cell layers will break lbc and rbc access
	tw::Float electrodeRadius,electrodePotential,slewRate;
	
	Electrostatic(const std::string& name,Grid* theGrid);
	virtual void Initialize();
	virtual void ExchangeResources();
	virtual void Reset();
	virtual void Update();
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);

	virtual void SetupInitialPotential();
	virtual void SetupElectrodePotential();
	virtual void ComputeFinalFields();

	virtual void EnergyHeadings(std::ofstream& outFile);
	virtual void EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn);
	virtual void BoxDiagnosticHeader(GridDataDescriptor*);
	virtual void BoxDiagnose(GridDataDescriptor*);
	virtual void PointDiagnosticHeader(std::ofstream& outFile);
	virtual void PointDiagnose(std::ofstream& outFile,const weights_3D& w);
};
