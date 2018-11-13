struct FieldSolver:Module
{
	EllipticSolver *ellipticSolver;

	FieldSolver(Grid* theGrid);
	~FieldSolver();
	virtual void Initialize();
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);
};

struct Electromagnetic:FieldSolver
{
	Field F,sources;
	ScalarField scratch1,scratch2,conductorMask;

	tw::vec3 dipoleCenter;
	tw::Float gammaBeam;
	std::vector<FarFieldDetectorDescriptor*> farFieldDetector;

	Electromagnetic(Grid* theGrid);
	virtual void ExchangeResources();
	virtual void Initialize();
	virtual void Reset();
	virtual void Update();
	virtual void MoveWindow();
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);

	template <tw::Int X,tw::Int Y,tw::Int Z>
	void LoadVectorPotential(Field& A,tw::Float t);
	template <tw::Int X,tw::Int Y,tw::Int Z>
	void CleanDivergence(Field& A,tw::Float charge_multiplier);
	void InitializeConductors();
	void SetExteriorBoundaryConditionsE(Field& A,const Element& ex,const Element& ey,const Element& ez);
	void SetExteriorBoundaryConditionsB(Field& A,const Element& bx,const Element& by,const Element& bz);
	void ForceQuasistaticVectorPotential(Field& A,ScalarField& DtPhi);

	virtual void StartDiagnostics();
	virtual void EnergyHeadings(std::ofstream& outFile);
	virtual void EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn);
	virtual void BoxDiagnosticHeader(GridDataDescriptor*);
	virtual void BoxDiagnose(GridDataDescriptor*);
	virtual void PointDiagnosticHeader(std::ofstream& outFile);
	virtual void PointDiagnose(std::ofstream& outFile,const weights_3D& w);
	virtual void CustomDiagnose();
};

struct CoulombSolver:Electromagnetic
{
	Field A4;
	LindmanBoundary L1,L2,R1,R2;

	CoulombSolver(Grid* theGrid);
	virtual void ExchangeResources();
	virtual void Initialize();
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);

	virtual void Update();
	virtual void ComputeFinalFields();

	virtual void BoxDiagnosticHeader(GridDataDescriptor*);
	virtual void BoxDiagnose(GridDataDescriptor*);
	virtual void PointDiagnosticHeader(std::ofstream& outFile);
	virtual void PointDiagnose(std::ofstream& outFile,const weights_3D& w);
};

struct DirectSolver:Electromagnetic
{
	Field A;
	Field PMLx,PMLy,PMLz;
	tw::Int layerThicknessX0,layerThicknessY0,layerThicknessZ0;
	tw::Int layerThicknessX1,layerThicknessY1,layerThicknessZ1;
	tw::Float reflectionCoefficient;
	bool enforceChargeConservation;

	YeePropagatorPML *yeeTool;

	DirectSolver(Grid* theGrid);
	~DirectSolver();
	virtual void Initialize();
	virtual void Update();
	virtual void MoveWindow();
	virtual void ReadInputFileDirective(std::stringstream& inputString,const std::string& command);
	virtual void ReadData(std::ifstream& inFile);
	virtual void WriteData(std::ofstream& outFile);

	void SetupPML(Field& pml,tw::Int g0,tw::Int gN,tw::Int L0,tw::Int L1,tw::Float R,tw::Float ds);
};

struct CurvilinearDirectSolver:DirectSolver
{
	CurvilinearDirectSolver(Grid* theGrid);
	virtual void Initialize();
	virtual void Update();
	virtual void SetSingularPointsE();
	virtual void SetSingularPointsB();
};

template <tw::Int X,tw::Int Y,tw::Int Z>
void Electromagnetic::LoadVectorPotential(Field& A,tw::Float t)
{
	#pragma omp parallel
	{
		tw::Int i,j,k,s;
		tw::vec3 r1,r2,r3;
		for (CellIterator cell(A,true);cell<cell.end();++cell)
		{
			cell.Decode(&i,&j,&k);
			r1.x = owner->X(i,1) - 0.5*owner->dX(i,1);
			r1.y = owner->X(j,2);
			r1.z = owner->X(k,3);
			r2.x = owner->X(i,1);
			r2.y = owner->X(j,2) - 0.5*owner->dX(j,2);
			r2.z = owner->X(k,3);
			r3.x = owner->X(i,1);
			r3.y = owner->X(j,2);
			r3.z = owner->X(k,3) - 0.5*owner->dX(k,3);
			owner->CurvilinearToCartesian(&r1);
			owner->CurvilinearToCartesian(&r2);
			owner->CurvilinearToCartesian(&r3);
			A(i,j,k,X) = 0.0;
			A(i,j,k,Y) = 0.0;
			A(i,j,k,Z) = 0.0;
			for (s=0;s<owner->wave.size();s++)
			{
				A(i,j,k,X) += owner->wave[s]->VectorPotential(t,r1).x;
				A(i,j,k,Y) += owner->wave[s]->VectorPotential(t,r2).y;
				A(i,j,k,Z) += owner->wave[s]->VectorPotential(t,r3).z;
			}
		}
	}
	CleanDivergence<X,Y,Z>(A,0.0);
}

template <tw::Int X,tw::Int Y,tw::Int Z>
void Electromagnetic::CleanDivergence(Field& A,tw::Float charge_multiplier)
{
	tw::Int i,j,k;

	#pragma omp parallel for private(i,j,k) collapse(3) schedule(static)
	for (i=1;i<=dim[1];i++)
		for (j=1;j<=dim[2];j++)
			for (k=1;k<=dim[3];k++)
			{
				scratch1(i,j,k) = charge_multiplier * sources(i,j,k,0);
				scratch1(i,j,k) -= divE<X,Y,Z>(A,i,j,k,*owner);
			}

	ellipticSolver->SaveBoundaryConditions();
	ellipticSolver->SetBoundaryConditions(scratch2,neumannWall,neumannWall,neumannWall,neumannWall,natural,natural);
	ellipticSolver->Solve(scratch2,scratch1,1.0);
	scratch2.ApplyBoundaryCondition();
	ellipticSolver->RestoreBoundaryConditions();

	add_grad<0,X,Y,Z>(scratch2,A,*owner,1.0);

	A.CopyFromNeighbors(Element(X,Z));
	A.ApplyBoundaryCondition(Element(X,Z));
}
