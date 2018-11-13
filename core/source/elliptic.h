struct EllipticSolver:ComputeTool
{
	ScalarField *coeff;
	boundarySpec x0,x1,y0,y1,z0,z1;
	boundarySpec x0s,x1s,y0s,y1s,z0s,z1s; // saved BC's
	std::valarray<tw::Float> lbc,rbc,lbc_t,rbc_t;
	tw::Float gammaBeam;

	EllipticSolver(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void Initialize();
	virtual void SetCoefficients(ScalarField *coefficients);
	virtual void SetBoundaryConditions(ScalarField& phi);
	virtual void SetBoundaryConditions(boundarySpec x0,boundarySpec x1,boundarySpec y0,boundarySpec y1,boundarySpec z0,boundarySpec z1);
	virtual void SaveBoundaryConditions();
	virtual void RestoreBoundaryConditions();
	virtual void FixPotential(ScalarField& phi,Region* theRegion,const tw::Float& thePotential);
	virtual void TransformBoundaryValues() {;}
	virtual void ZeroModeGhostCellValues(tw::Float *phi0,tw::Float *phiN1,ScalarField& rho,tw::Float mul);
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul) = 0;
	void FormOperatorStencil(std::valarray<tw::Float>& D,tw::Int i,tw::Int j,tw::Int k);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct IterativePoissonSolver:EllipticSolver
{
	char *mask1,*mask2;
	tw::Int iterationsPerformed;
	tw::Float normResidualAchieved,normSource;
	tw::Float overrelaxationChange;
	tw::Int maxIterations;
	tw::Float tolerance,overrelaxation,minimumNorm;

	IterativePoissonSolver(const std::string& name,MetricSpace *m,Task *tsk);
	~IterativePoissonSolver();
	virtual void Initialize();
	virtual void FixPotential(ScalarField& phi,Region* theRegion,const tw::Float& thePotential);
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
	virtual void StatusMessage(std::ostream *theStream);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct EllipticSolver1D:EllipticSolver
{
	GlobalIntegrator<tw::Float> *globalIntegrator;

	EllipticSolver1D(const std::string& name,MetricSpace *m,Task *tsk);
	~EllipticSolver1D();
	virtual void Initialize();
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
};

struct PoissonSolver:EllipticSolver
{
	GlobalIntegrator<tw::Float> *globalIntegrator;

	PoissonSolver(const std::string& name,MetricSpace *m,Task *tsk);
	~PoissonSolver();
	virtual void Initialize();
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
};

struct EigenmodePoissonSolver:EllipticSolver
{
	std::valarray<tw::Float> eigenvalue,hankel,inverseHankel;
	GlobalIntegrator<tw::Float> *globalIntegrator;

	EigenmodePoissonSolver(const std::string& name,MetricSpace *m,Task *tsk);
	~EigenmodePoissonSolver();
	virtual void Initialize();
	virtual void TransformBoundaryValues();
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
};
