namespace tw
{
	enum class module_type {	nullModule,
					electrostatic,
					coulombSolver,directSolver,curvilinearDirectSolver,
					qsLaser,pgcLaser,
					kinetics,species,fluidFields,equilibriumGroup,chemical,chemistry,
					boundElectrons,schroedinger,pauli,kleinGordon,dirac};
}

struct Module:DiscreteSpace
{
	tw::Float dt,dth,dti;
	std::string name;
	std::vector<Profile*> profile;
	Grid* owner;
	Module* super;
	std::vector<Module*> submodule;
	std::vector<ComputeTool*> moduleTool;

	tw::Int updateSequencePriority;
	tw::module_type typeCode;
	bool suppressNextUpdate;

	// OpenCL Support
	private:
	std::string programFilename;
	std::string buildLog;
	public:
	#ifdef USE_OPENCL
	cl_program program;
	#endif

	Module(const std::string& name,Grid* theGrid);
	~Module();
	bool AddSubmodule(Module* sub);
	virtual bool ValidSubmodule(Module* sub);
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

	virtual void VerifyInput();
	virtual void ReadInputFileBlock(std::stringstream& inputString);
	virtual bool ReadQuasitoolBlock(const std::vector<std::string>& preamble,std::stringstream& inputString);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);

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

	static bool SingularType(tw::module_type theType);
	static tw::module_type CreateSupermoduleTypeFromSubmoduleKey(const std::string& key);
	static bool QuasitoolNeedsModule(const std::vector<std::string>& preamble);
	static tw::module_type CreateTypeFromInput(const std::vector<std::string>& preamble);
	static Module* CreateObjectFromType(const std::string& name,tw::module_type theType,Grid* theGrid);
	static Module* CreateObjectFromFile(std::ifstream& inFile,Grid* thGrid);
};

struct ModuleComparator
{
	bool operator() (Module* const& m1,Module* const& m2)
	{
		return m1->updateSequencePriority < m2->updateSequencePriority;
	}
};
