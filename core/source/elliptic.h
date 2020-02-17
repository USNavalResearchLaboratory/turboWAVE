// IterativePoissonSolver and EllipticSolver1D handle equations of the form div(coeff*grad(phi)) = mul*source
// PoissonSolver and EigenmodePoissonSolver only handle div(phi) = mul*source
// For PoissonSolver and EigenmodePoissonSolver, only z-boundary conditions are programmable

// Boundary Conditions:
// These tools expect boundary values to be stored in the outer ghost cells
// E.g., if phi(i,j,-1) = V0 then:
// BC Type          Explicit BC
// -------------    ------------
// dirichletCell    phi(i,j,0) = V0
// neumannWall      phi(i,j,1) - phi(i,j,0) = V0
// dirichletWall    phi(i,j,0)/2 + phi(i,j,1)/2 = V0

// N.b. the BC handling assumes at least 2 ghost cell layers in the field

struct EllipticSolver:ComputeTool
{
	ScalarField *coeff;
	boundarySpec x0,x1,y0,y1,z0,z1;
	boundarySpec x0s,x1s,y0s,y1s,z0s,z1s; // saved BC's
	tw::Float gammaBeam;

	EllipticSolver(const std::string& name,MetricSpace *m,Task *tsk);
	virtual void SetCoefficients(ScalarField *coefficients);
	virtual void SetBoundaryConditions(ScalarField& phi);
	virtual void SetBoundaryConditions(boundarySpec x0,boundarySpec x1,boundarySpec y0,boundarySpec y1,boundarySpec z0,boundarySpec z1);
	virtual void SaveBoundaryConditions();
	virtual void RestoreBoundaryConditions();
	virtual void FixPotential(ScalarField& phi,Region* theRegion,const tw::Float& thePotential);
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
	tw::Int maxIterations;
	tw::Float tolerance,overrelaxation,minimumNorm;

	IterativePoissonSolver(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~IterativePoissonSolver();
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
	virtual ~EllipticSolver1D();
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
};

struct PoissonSolver:EllipticSolver
{
	GlobalIntegrator<tw::Float> *globalIntegrator;

	PoissonSolver(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~PoissonSolver();
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
};

struct EigenmodePoissonSolver:EllipticSolver
{
	std::valarray<tw::Float> eigenvalue,hankel,inverseHankel;
	GlobalIntegrator<tw::Float> *globalIntegrator;

	EigenmodePoissonSolver(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~EigenmodePoissonSolver();
	virtual void Initialize();
	virtual void Solve(ScalarField& phi,ScalarField& source,tw::Float mul);
};
