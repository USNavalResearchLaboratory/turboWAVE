struct EllipticSolver:ComputeTool
{
	ScalarField *coeff;
	boundarySpec x0,x1,y0,y1,z0,z1;
	boundarySpec x0s,x1s,y0s,y1s,z0s,z1s; // saved BC's
	std::valarray<tw::Float> lbc,rbc,lbc_t,rbc_t;
	tw::Int maxIterations;
	tw::Float tolerance,gammaBeam,overrelaxation,minimumNorm;
	
	EllipticSolver(MetricSpace *m,Task *tsk,bool shared);
	virtual void Initialize();
	virtual void SetCoefficients(ScalarField *coefficients);
	virtual void SetBoundaryConditions(ScalarField& phi,boundarySpec x0,boundarySpec x1,boundarySpec y0,boundarySpec y1,boundarySpec z0,boundarySpec z1);
	virtual void SaveBoundaryConditions();
	virtual void RestoreBoundaryConditions();
	virtual void FixPotential(ScalarField& phi,Region* theRegion,const tw::Float& thePotential);
	virtual void TransformBoundaryValues() {;}
	virtual void ZeroModeGhostCellValues(tw::Float *phi0,tw::Float *phiN1,ScalarField& rho,tw::Float mul);
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul) = 0;
	void FormOperatorStencil(std::valarray<tw::Float>& D,tw::Int i,tw::Int j,tw::Int k);
};

struct IterativePoissonSolver:EllipticSolver
{
	char *mask1,*mask2;
	tw::Int iterationsPerformed;
	tw::Float normResidualAchieved,normSource;
	tw::Float overrelaxationChange;
	
	IterativePoissonSolver(MetricSpace *m,Task *tsk,bool shared);
	~IterativePoissonSolver();
	virtual void Initialize();
	virtual void FixPotential(ScalarField& phi,Region* theRegion,const tw::Float& thePotential);
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
	
	//void Coarser(ScalarField& RFine,ScalarField& R,ScalarField& eps);
	//void Finer(ScalarField& correction);

	virtual void StatusMessage(std::ostream *theStream);
};

struct EllipticSolver1D:EllipticSolver
{
	GlobalIntegrator<tw::Float> *globalIntegrator;

	EllipticSolver1D(MetricSpace *m,Task *tsk,bool shared);
	~EllipticSolver1D();
	virtual void Initialize();
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
};

struct PoissonSolver:EllipticSolver
{
	GlobalIntegrator<tw::Float> *globalIntegrator;
	
	PoissonSolver(MetricSpace *m,Task *tsk,bool shared);
	~PoissonSolver();
	virtual void Initialize();
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
};

struct EigenmodePoissonSolver:EllipticSolver
{
	std::valarray<tw::Float> eigenvalue,hankel,inverseHankel;
	GlobalIntegrator<tw::Float> *globalIntegrator;
	
	EigenmodePoissonSolver(MetricSpace *m,Task *tsk,bool shared);
	~EigenmodePoissonSolver();
	virtual void Initialize();
	virtual void TransformBoundaryValues();
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
};

