namespace tw
{
	enum class module_type {	none,
					electrostatic,
					coulombSolver,directSolver,curvilinearDirectSolver,
					qsLaser,pgcLaser,
					kinetics,species,fluidFields,equilibriumGroup,chemical,sparcHydroManager,
					boundElectrons,schroedinger,pauli,kleinGordon,dirac};
	enum class priority { diagnostic, source, field };
	inline std::map<tw::priority,tw::Int> priority_sort_map()
	{
		std::map<tw::priority,tw::Int> ans = {{priority::diagnostic,100},{priority::source,200},{priority::field,300}};
		return ans;
	}
}

struct Module:DiscreteSpace
{
	std::string name;
	Simulation* owner;
	Module* super;
	std::vector<Module*> submodule;
	std::vector<ComputeTool*> moduleTool;
	tw::input::DirectiveReader directives;

	tw::priority updateSequencePriority;
	tw::Int subSequencePriority;
	tw::Int smoothing[4],compensation[4];
	tw::module_type typeCode;
	bool suppressNextUpdate;

	// Strongly typed ComputeTool lists that are provided for free
	std::vector<Profile*> profile;
	std::vector<Wave*> wave;
	std::vector<Conductor*> conductor;

	// OpenCL Support
	private:
	std::string programFilename;
	std::string buildLog;
	public:
	#ifdef USE_OPENCL
	cl_program program;
	#endif

	Module(const std::string& name,Simulation* sim);
	virtual ~Module();
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
	virtual bool ReadQuasitoolBlock(const tw::input::Preamble& preamble,std::stringstream& inputString);
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void StartDiagnostics();
	virtual void Report(Diagnostic&);
	virtual void WarningMessage(std::ostream *theStream);
	virtual void StatusMessage(std::ostream *theStream) {;}

	static std::map<std::string,tw::module_type> Map();
	static bool SingularType(tw::module_type theType);
	static tw::module_type CreateSupermoduleTypeFromSubmoduleKey(const std::string& key);
	static bool QuasitoolNeedsModule(const tw::input::Preamble& preamble);
	static tw::module_type CreateTypeFromInput(const tw::input::Preamble& preamble);
	static Module* CreateObjectFromType(const std::string& name,tw::module_type theType,Simulation* sim);
	static Module* CreateObjectFromFile(std::ifstream& inFile,Simulation* sim);
};

struct ModuleComparator
{
	bool operator() (Module* const& m1,Module* const& m2)
	{
		tw::Int idx1 = tw::priority_sort_map()[m1->updateSequencePriority] + m1->subSequencePriority;
		tw::Int idx2 = tw::priority_sort_map()[m2->updateSequencePriority] + m2->subSequencePriority;
		return idx1 < idx2;
	}
};
