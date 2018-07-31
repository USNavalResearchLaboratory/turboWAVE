enum tw_module {	nullModule,
					electrostatic,
					coulombSolver,directSolver,curvilinearDirectSolver,
					qsLaser,pgcLaser,
					kinetics,species,fluidFields,equilibriumGroup,chemical,chemistry,
					boundElectrons,atomicPhysics,schroedingerModule,pauliModule,kleinGordon,dirac};

struct Module:DiscreteSpace
{
	tw::Float dt,dth,dti;
	std::string name;
	std::vector<Profile*> profile;
	Grid* owner;

	tw::Int updateOrderingIndex;
	tw_module typeCode;
	bool suppressNextUpdate;

	// OpenCL Support
	private:
	std::string programFilename;
	std::string buildLog;
	public:
	#ifdef USE_OPENCL
	cl_program program;
	#endif

	Module(Grid* theGrid);
	~Module();
	virtual void PublishResource(void* resource,const std::string& description);
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void ExchangeResources() {;}
	virtual void Initialize();
	virtual void Reset() {;}
	virtual void Update() {;}
	virtual void MoveWindow() {;}
	virtual void AntiMoveWindow() {;}
	virtual void AdaptGrid() {;}
	virtual tw::Float AdaptTimestep() { return 0.0;}

	void InitializeCLProgram(const std::string& filename);

	virtual void ReadInputFileBlock(std::stringstream& inputString);
	virtual void ReadInputFileTerm(std::stringstream& inputString,std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
	static Module* CreateObjectFromFile(std::ifstream& inFile,Grid* thGrid);

	virtual void StartDiagnostics();
	virtual void EnergyHeadings(std::ofstream& outFile) {;}
	virtual void EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn) {;}
	virtual void BoxDiagnosticHeader(GridDataDescriptor*) {;}
	virtual void BoxDiagnose(GridDataDescriptor*) {;}
	virtual void PointDiagnosticHeader(std::ofstream& outFile) {;}
	virtual void PointDiagnose(std::ofstream& outFile,const weights_3D& w) {;}
	virtual void CustomDiagnose() {;}
	virtual void WarningMessage(std::ostream *theStream);
	virtual void StatusMessage(std::ostream *theStream) {;}
};

struct ModuleComparator
{
	bool operator() (Module* const& m1,Module* const& m2)
	{
		return m1->updateOrderingIndex < m2->updateOrderingIndex;
	}
};
