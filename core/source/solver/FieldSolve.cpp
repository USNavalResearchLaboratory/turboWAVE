module;

#include "tw_includes.h"
#include "tw_logger.h"

export module field_solve;
import input;
import compute_tool;
import twmodule;
import fields;
import elliptic;
import hyperbolic;
import diagnostics;
import injection;
import functions;
import logger;

using namespace tw::bc;

export struct FieldSolver:Module
{
	std::shared_ptr<EllipticSolver> ellipticSolver;
	tw::Waves waves;
	tw::Conductors conductors;

	FieldSolver(const std::string& name,Simulation* sim);
	virtual void VerifyInput();
};

export struct Electromagnetic:FieldSolver
{
	Field F,sources;
	ScalarField scratch1,scratch2,conductorMask;

	tw::vec3 dipoleCenter,externalE,externalB;
	tw::Float gammaBeam;

	Electromagnetic(const std::string& name,Simulation* sim);
	virtual void ExchangeResources();
	virtual void Initialize();
	virtual void Reset();
	virtual void Update();
	virtual void MoveWindow();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	template <tw::Int X,tw::Int Y,tw::Int Z>
	void LoadVectorPotential(Field& A,tw::Float t);
	template <tw::Int X,tw::Int Y,tw::Int Z>
	void CleanDivergence(Field& A,tw::Float charge_multiplier);
	void InitializeConductors();
	void SetExteriorBoundaryConditionsE(Field& A,const Rng& ex,const Rng& ey,const Rng& ez);
	void SetExteriorBoundaryConditionsB(Field& A,const Rng& bx,const Rng& by,const Rng& bz);
	void ForceQuasistaticVectorPotential(Field& A,ScalarField& DtPhi);

	virtual void StartDiagnostics();
	virtual void Report(Diagnostic&);
};

export struct CoulombSolver:Electromagnetic
{
	Field A4;
	LindmanBoundary L1,L2,R1,R2;

	CoulombSolver(const std::string& name,Simulation* sim);
	virtual void ExchangeResources();
	virtual void Initialize();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void Update();
	virtual void ComputeFinalFields();
	virtual void Report(Diagnostic&);
};

export struct DirectSolver:Electromagnetic
{
	Field A;
	Field PMLx,PMLy,PMLz;
	tw::Int layerThickness[6];
	tw::Float reflectionCoefficient[6];
	bool enforceChargeConservation;

	std::shared_ptr<YeePropagatorPML> yeeTool;

	DirectSolver(const std::string& name,Simulation* sim);
	virtual void Initialize();
	virtual void Update();
	virtual void MoveWindow();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	void SetupPML(Field& pml,tw::Int g0,tw::Int gN,tw::Int L0,tw::Int L1,tw::Float R0,tw::Float R1,tw::Float ds);
};

export struct CurvilinearDirectSolver:DirectSolver
{
	CurvilinearDirectSolver(const std::string& name,Simulation* sim);
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
		for (auto cell : EntireCellRange(A,1))
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
			A(1,i,j,k,X) = 0.0;
			A(1,i,j,k,Y) = 0.0;
			A(1,i,j,k,Z) = 0.0;
			for (auto wave : waves)
			{
				A(1,i,j,k,X) += wave->VectorPotential(t,r1).x;
				A(1,i,j,k,Y) += wave->VectorPotential(t,r2).y;
				A(1,i,j,k,Z) += wave->VectorPotential(t,r3).z;
			}
		}
	}
	CleanDivergence<X,Y,Z>(A,0.0);
}

template <tw::Int X,tw::Int Y,tw::Int Z>
void Electromagnetic::CleanDivergence(Field& A,tw::Float charge_multiplier)
{
	using namespace tw::bc;

	// TODO: #pragma omp parallel for collapse(3) schedule(static)
	for (tw::Int i=1;i<=dim[1];i++)
		for (tw::Int j=1;j<=dim[2];j++)
			for (tw::Int k=1;k<=dim[3];k++)
			{
				scratch1(i,j,k) = charge_multiplier * sources(1,i,j,k,0);
				scratch1(i,j,k) -= divE<X,Y,Z>(A,1,i,j,k,*owner);
			}

	ellipticSolver->SaveBoundaryConditions();
	ellipticSolver->SetBoundaryConditions(fld::neumannWall,fld::neumannWall,fld::neumannWall,fld::neumannWall,fld::natural,fld::natural);
	ellipticSolver->SetFieldsBoundaryConditions(scratch2,Rng(0));
	ellipticSolver->Solve(scratch2,scratch1,1.0);
	scratch2.ApplyBoundaryCondition();
	ellipticSolver->RestoreBoundaryConditions();

	add_grad<0,X,Y,Z>(1,scratch2,A,*owner,1.0);

	A.CopyFromNeighbors(Rng(X,Z+1));
	A.ApplyBoundaryCondition(Rng(X,Z+1));
}

export struct FarFieldDiagnostic : Module
{
	// In the far field space the index space goes like (t,theta,phi)
	tw::Float radius,bounds[6];
	tw::Int dims[4];
	tw::Int period;

	Vec3Field A;
	Field *J4;

	FarFieldDiagnostic(const std::string& name,Simulation *sim);
	virtual void Initialize();
	virtual bool InspectResource(void *resource,const std::string& description);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
	virtual void Update();
	void FinalReport(); // TODO: need a way to trigger, add module method?
};

/////////////////////////////
// FIELD SOLVER BASE CLASS //
/////////////////////////////

// The base class is essentially for managing the ellipticSolver tool.
// This serves as a basic template for how to manage a ComputeTool.

FieldSolver::FieldSolver(const std::string& name,Simulation* sim):Module(name,sim)
{
	if (native.native!=tw::units::plasma)
		throw tw::FatalError("FieldSolver module requires <native units = plasma>");

	updateSequencePriority = tw::priority::field;
}

void FieldSolver::VerifyInput()
{
	Module::VerifyInput();
	// Populate strongly typed tools
	for (auto tool : tools)
	{
		if (std::dynamic_pointer_cast<EllipticSolver>(tool)) {
			ellipticSolver = std::dynamic_pointer_cast<EllipticSolver>(tool);
		} else if (std::dynamic_pointer_cast<Wave>(tool)) {
			waves.push_back(std::dynamic_pointer_cast<Wave>(tool));
		} else if (std::dynamic_pointer_cast<Conductor>(tool)) {
			conductors.push_back(std::dynamic_pointer_cast<Conductor>(tool));
		}
	}
	// If no elliptic solver, create one automatically
	if (!ellipticSolver) {
		auto name = owner->CreateTool("default_facr_poisson_solver",tw::tool_type::facrPoissonSolver);
		ellipticSolver = std::dynamic_pointer_cast<EllipticSolver>(owner->UseTool(name));
	}
}


////////////////////////////////
// ELECTROMAGNETIC BASE CLASS //
////////////////////////////////


Electromagnetic::Electromagnetic(const std::string& name,Simulation* sim):FieldSolver(name,sim)
{
	if (dk(0) <= std::sqrt(
		(sim->GlobalDim(1)==1 ? 0.0 : 1.0)*sqr(dk(1)) +
		(sim->GlobalDim(2)==1 ? 0.0 : 1.0)*sqr(dk(2)) +
		(sim->GlobalDim(3)==1 ? 0.0 : 1.0)*sqr(dk(3))))
		throw tw::FatalError("Courant condition is violated " + 
			std::to_string(dk(0)) + " " + std::to_string(dk(3)));

	F.Initialize(6,*this,owner);
	sources.Initialize(4,*this,owner);
	scratch1.Initialize(*this,owner);
	scratch2.Initialize(*this,owner);
	conductorMask.Initialize(*this,owner);

	gammaBeam = 1.0;

	#ifdef USE_OPENCL
	F.InitializeComputeBuffer();
	sources.InitializeComputeBuffer();
	#endif

	directives.Add("dipole center",new tw::input::Vec3(&dipoleCenter),false);
	directives.Add("external E",new tw::input::Vec3(&externalE),false);
	directives.Add("external B",new tw::input::Vec3(&externalB),false);
	directives.Add("gamma beam",new tw::input::Float(&gammaBeam),false);
}

void Electromagnetic::ExchangeResources()
{
	PublishResource(&F,"electromagnetic:F");
	PublishResource(&sources,"electromagnetic:sources");
}

void Electromagnetic::Initialize()
{
	logger::TRACE("initialize EM base");

	FieldSolver::Initialize();
	InitializeConductors();

	SetExteriorBoundaryConditionsE(F,Rng(0),Rng(1),Rng(2));
	SetExteriorBoundaryConditionsB(F,Rng(3),Rng(4),Rng(5));

	sources.SetBoundaryConditions(Rng(0,4),tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	sources.SetBoundaryConditions(Rng(1),tw::grid::x,fld::normalFluxFixed,fld::normalFluxFixed);
	if (owner->gridGeometry==tw::grid::cylindrical)
		sources.SetBoundaryConditions(Rng(0),tw::grid::x,fld::neumannWall,fld::dirichletCell);

	sources.SetBoundaryConditions(Rng(0,4),tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	sources.SetBoundaryConditions(Rng(2),tw::grid::y,fld::normalFluxFixed,fld::normalFluxFixed);

	sources.SetBoundaryConditions(Rng(0,4),tw::grid::z,fld::dirichletCell,fld::dirichletCell);
	sources.SetBoundaryConditions(Rng(3),tw::grid::z,fld::normalFluxFixed,fld::normalFluxFixed);

	Electromagnetic::Update(); // globalize charge density (local deposition in source modules)
}

#ifdef USE_OPENCL
void Electromagnetic::Reset()
{
	sources.MADDComputeBuffer(0.0,0.0);
}
#else
void Electromagnetic::Reset()
{
	sources = 0.0;
}
#endif

#ifdef USE_OPENCL
void Electromagnetic::Update()
{
	sources.WeightComputeBufferByVolume(*owner,1.0);
	sources.UpdateGhostCellsInComputeBuffer();
}
#else
void Electromagnetic::Update()
{
	sources.DepositFromNeighbors(Rng(0,4));
	sources.ApplyFoldingCondition(Rng(0,4)); // only apply before conversion to density
	conserved_current_to_dens<0,1,2,3>(1,sources,*owner);
	sources.ApplyBoundaryCondition(Rng(0,4)); // only apply after conversion to density
	sources.Smooth(Rng(0,4),*owner,smoothing,compensation);
}
#endif

void Electromagnetic::MoveWindow()
{
	Module::MoveWindow();
	for (auto s : StripRange(*this,3,0,1,strongbool::yes))
		F.Shift(Rng(0,6),s,-1,0.0);
	F.DownwardCopy(Rng(0,6),tw::grid::z,1);
}

void Electromagnetic::InitializeConductors()
{
	// Conductor resides in cells shifted back by 1/2
	// These have E known along upper edges and B normal to upper walls
	// The conductor fills the whole cell or none of the cell
	tw::vec3 shiftedCenter;
	//#pragma omp parallel for private(i,j,k,s,shiftedCenter) collapse(3) schedule(static)
	for (auto cell : EntireCellRange(*this,1))
	{
		conductorMask(cell) = 1.0;
		shiftedCenter = owner->Pos(cell) - 0.5*owner->dPos(cell);
		for (auto c : conductors)
			if (c->affectsA && c->theRgn->Inside(shiftedCenter,*owner))
				conductorMask(cell) = 0.0;
	}
}

void Electromagnetic::SetExteriorBoundaryConditionsE(Field& A,const Rng& ex,const Rng& ey,const Rng& ez)
{
	// in the following, it doesn't hurt to zero out low-side cell walls in low-side ghost cells
	A.SetBoundaryConditions(ex,tw::grid::x,fld::dirichletCell,fld::none);
	A.SetBoundaryConditions(ex,tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	A.SetBoundaryConditions(ex,tw::grid::z,fld::dirichletCell,fld::dirichletCell);

	A.SetBoundaryConditions(ey,tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	A.SetBoundaryConditions(ey,tw::grid::y,fld::dirichletCell,fld::none);
	A.SetBoundaryConditions(ey,tw::grid::z,fld::dirichletCell,fld::dirichletCell);

	A.SetBoundaryConditions(ez,tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	A.SetBoundaryConditions(ez,tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	A.SetBoundaryConditions(ez,tw::grid::z,fld::dirichletCell,fld::none);
}

void Electromagnetic::SetExteriorBoundaryConditionsB(Field& A,const Rng& bx,const Rng& by,const Rng& bz)
{
	// in the following, it doesn't hurt to zero out low-side cell walls in low-side ghost cells
	A.SetBoundaryConditions(bx,tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	A.SetBoundaryConditions(bx,tw::grid::y,fld::dirichletCell,fld::none);
	A.SetBoundaryConditions(bx,tw::grid::z,fld::dirichletCell,fld::none);

	A.SetBoundaryConditions(by,tw::grid::x,fld::dirichletCell,fld::none);
	A.SetBoundaryConditions(by,tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	A.SetBoundaryConditions(by,tw::grid::z,fld::dirichletCell,fld::none);

	A.SetBoundaryConditions(bz,tw::grid::x,fld::dirichletCell,fld::none);
	A.SetBoundaryConditions(bz,tw::grid::y,fld::dirichletCell,fld::none);
	A.SetBoundaryConditions(bz,tw::grid::z,fld::dirichletCell,fld::dirichletCell);
}

void Electromagnetic::ForceQuasistaticVectorPotential(Field& A4,ScalarField& DtPhi)
{
	// only works in cartesian, assumes coulomb gauge
	tw::Int i,j,k,ax;
	tw::Float beta = std::sqrt(1.0 - 1.0/sqr(gammaBeam));
	tw::Float w = beta*dx(0)*dk(3);
	ScalarField ans,rhs;
	ans.Initialize(*this,owner);
	rhs.Initialize(*this,owner);

	ellipticSolver->gammaBeam = gammaBeam;

	for (ax=1;ax<=3;ax++)
	{
		// Find Solution at current time

		for (i=1;i<=dim[1];i++)
			for (j=1;j<=dim[2];j++)
			{
				tw::xstrip<3> v(*this,0,tw::node4{1,i,j,0});
				for (k=1;k<=dim[3];k++)
					rhs(v,k) = DtPhi.d1(v,k,0,ax) - sources(v,k,ax);
			}
		rhs.CopyFromNeighbors();

		ellipticSolver->Solve(ans,rhs,1.0);
		for (auto cell : EntireCellRange(*this,1))
			A4(cell,ax+4) = ans(cell);

		// Find Solution at last time

		for (i=1;i<=dim[1];i++)
			for (j=1;j<=dim[2];j++)
			{
				tw::xstrip<3> v(*this,0,tw::node4{1,i,j,0});
				for (k=1;k<=dim[3];k++)
					rhs(v,k) = DtPhi.d1(v,k,0,ax) - sources(v,k,ax);
			}
		rhs.CopyFromNeighbors();

		for (k=lfg[3];k<=dim[3];k++)
			for (j=lfg[2];j<=ufg[2];j++)
				for (i=lfg[1];i<=ufg[1];i++)
					rhs(i,j,k) = (1.0 - w)*rhs(i,j,k) + w*rhs(i,j,k+1);
		rhs.DownwardCopy(tw::grid::z,1);

		ellipticSolver->Solve(ans,rhs,1.0);
		for (auto cell : EntireCellRange(*this,1))
			A4(cell,ax) = ans(cell);
	}

	ellipticSolver->gammaBeam = 1.0;
}

void Electromagnetic::ReadCheckpoint(std::ifstream& inFile)
{
	FieldSolver::ReadCheckpoint(inFile);
	inFile.read((char*)&dipoleCenter,sizeof(tw::vec3));
	inFile.read((char*)&gammaBeam,sizeof(tw::Float));
	F.ReadCheckpoint(inFile);
}

void Electromagnetic::WriteCheckpoint(std::ofstream& outFile)
{
	FieldSolver::WriteCheckpoint(outFile);
	outFile.write((char*)&dipoleCenter,sizeof(tw::vec3));
	outFile.write((char*)&gammaBeam,sizeof(tw::Float));
	F.WriteCheckpoint(outFile);
}

void Electromagnetic::StartDiagnostics()
{
	#ifdef USE_OPENCL
	F.ReceiveFromComputeBuffer();
	sources.ReceiveFromComputeBuffer();
	#endif
}

void Electromagnetic::Report(Diagnostic& diagnostic)
{
	FieldSolver::Report(diagnostic);

	ScalarField temp;
	temp.Initialize(*this,owner);

	for (auto cell : InteriorCellRange(*this,1))
		temp(cell) = 0.5*(Norm(F.Vec3(cell,0))+Norm(F.Vec3(cell,3)));
	diagnostic.VolumeIntegral("FieldEnergy",temp,1,0);

	diagnostic.VolumeIntegral("TotalCharge",sources,1,0);
	diagnostic.FirstMoment("Dx",sources,1,0,dipoleCenter,tw::grid::x);
	diagnostic.FirstMoment("Dy",sources,1,0,dipoleCenter,tw::grid::y);
	diagnostic.FirstMoment("Dz",sources,1,0,dipoleCenter,tw::grid::z);
	diagnostic.ReportField("Ex",F,1,0,tw::dims::electric_field,"$E_x$");
	diagnostic.ReportField("Ey",F,1,1,tw::dims::electric_field,"$E_y$");
	diagnostic.ReportField("Ez",F,1,2,tw::dims::electric_field,"$E_z$");
	diagnostic.ReportField("Bx",F,1,3,tw::dims::magnetic_field,"$B_x$");
	diagnostic.ReportField("By",F,1,4,tw::dims::magnetic_field,"$B_y$");
	diagnostic.ReportField("Bz",F,1,5,tw::dims::magnetic_field,"$B_z$");
	diagnostic.ReportField("rho",sources,1,0,tw::dims::charge_density,"$\\rho$");
	diagnostic.ReportField("Jx",sources,1,1,tw::dims::current_density,"$j_x$");
	diagnostic.ReportField("Jy",sources,1,2,tw::dims::current_density,"$j_y$");
	diagnostic.ReportField("Jz",sources,1,3,tw::dims::current_density,"$j_z$");
}


////////////////////////////////////////////////
// COULOMB GAUGE FIELD SOLVER (original WAVE) //
////////////////////////////////////////////////


CoulombSolver::CoulombSolver(const std::string& name,Simulation* sim):Electromagnetic(name,sim)
{
	A4.Initialize(8,*this,owner);

	L1.Initialize(sim,sim,tw::grid::z,tw::grid::low,1,5);
	L2.Initialize(sim,sim,tw::grid::z,tw::grid::low,2,6);
	R1.Initialize(sim,sim,tw::grid::z,tw::grid::high,1,5);
	R2.Initialize(sim,sim,tw::grid::z,tw::grid::high,2,6);
}

void CoulombSolver::ExchangeResources()
{
	Electromagnetic::ExchangeResources();
}

void CoulombSolver::Initialize()
{
	Electromagnetic::Initialize();

	// Overrides the boundary conditions set by the tool
	ellipticSolver->SetBoundaryConditions(fld::periodic,fld::periodic,fld::periodic,fld::periodic,fld::natural,fld::natural);
	ellipticSolver->SetFieldsBoundaryConditions(scratch2,Rng(0));

	SetExteriorBoundaryConditionsE(A4,Rng(1),Rng(2),Rng(3));
	SetExteriorBoundaryConditionsE(A4,Rng(5),Rng(6),Rng(7));

	// Initialize radiation fields

	LoadVectorPotential<1,2,3>(A4,-0.5*dx(0));
	LoadVectorPotential<5,6,7>(A4,0.5*dx(0));

	// Initialize static fields

	CopyFieldData(scratch1,Rng(0),sources,Rng(0));
	ellipticSolver->Solve(scratch2,scratch1,-1.0);
	CopyFieldData(A4,Rng(0),scratch2,Rng(0));
	CopyFieldData(A4,Rng(4),scratch2,Rng(0));

	ComputeFinalFields();
}

void CoulombSolver::ReadCheckpoint(std::ifstream& inFile)
{
	Electromagnetic::ReadCheckpoint(inFile);

	A4.ReadCheckpoint(inFile);
	L1.ReadCheckpoint(inFile);
	L2.ReadCheckpoint(inFile);
	R1.ReadCheckpoint(inFile);
	R2.ReadCheckpoint(inFile);
}

void CoulombSolver::WriteCheckpoint(std::ofstream& outFile)
{
	Electromagnetic::WriteCheckpoint(outFile);

	A4.WriteCheckpoint(outFile);
	L1.WriteCheckpoint(outFile);
	L2.WriteCheckpoint(outFile);
	R1.WriteCheckpoint(outFile);
	R2.WriteCheckpoint(outFile);
}

void CoulombSolver::Update()
{
	Electromagnetic::Update();
	const tw::Float dt = dx(0);
	const tw::Float dth = 0.5*dt;
	const tw::Float dti = dk(0);

	// upon entry : A0(-1/2) , A1(1/2) , phi0(-1) , phi1(0) , J(1/2) , rho(1)
	// upon exit : A0(1/2) , A1(3/2) , phi0(0) , phi1(1)

	// ADVANCE THE POTENTIALS

	// Update scalar potential and compute dphi/dt

	CopyFieldData(A4,Rng(0),A4,Rng(4));
	CopyFieldData(scratch1,Rng(0),sources,Rng(0));
	ellipticSolver->Solve(scratch2,scratch1,-1.0);
	CopyFieldData(A4,Rng(4),scratch2,Rng(0));
	AddMulFieldData(scratch2,Rng(0),A4,Rng(0),-1.0);
	scratch2 *= dti;

	// Deal with beam initialization

	//if (gammaBeam!=1.0 && owner->WindowPos(0)==0.0)
	//	ForceQuasistaticVectorPotential(A4,scratch2);

	// Must update boundary memory before vector potential

	L1.UpdateBoundaryMemory(A4,dt);
	L2.UpdateBoundaryMemory(A4,dt);
	R1.UpdateBoundaryMemory(A4,dt);
	R2.UpdateBoundaryMemory(A4,dt);

	// Update the interior : do this in place by putting new data in A0

	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,0,1,false))
		{
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				A4(v,k,1) = 2.0*A4(v,k,5) - A4(v,k,1) + dt*dt*(sources(v,k,1) + A4.d2(v,k,5,1) + A4.d2(v,k,5,2) + A4.d2(v,k,5,3) - scratch2.dbak(v,k,0,1));
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				A4(v,k,2) = 2.0*A4(v,k,6) - A4(v,k,2) + dt*dt*(sources(v,k,2) + A4.d2(v,k,6,1) + A4.d2(v,k,6,2) + A4.d2(v,k,6,3) - scratch2.dbak(v,k,0,2));
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				A4(v,k,3) = 2.0*A4(v,k,7) - A4(v,k,3) + dt*dt*(sources(v,k,3) + A4.d2(v,k,7,1) + A4.d2(v,k,7,2) + A4.d2(v,k,7,3) - scratch2.dbak(v,k,0,3));
		}
	}

	// Boundary Conditions
	// x,y components from Lindman absorbing boundary conditions
	// z component from gauge condition

	A4.ApplyBoundaryCondition(Rng(1,4));
	A4.CopyFromNeighbors(Rng(1,4));

	if (owner->n0[3]==MPI_PROC_NULL)
	{
		L1.Set(A4,&waves,owner->WindowPos(0)+dth,dt);
		L2.Set(A4,&waves,owner->WindowPos(0)+dth,dt);
		for (auto strip : StripRange(A4,3,0,1,strongbool::no))
			A4(strip,0,3) = A4(strip,1,3) + dx(3)*(A4.dfwd(strip,0,1,1) + A4.dfwd(strip,0,2,2));
	}

	if (owner->n1[3]==MPI_PROC_NULL)
	{
		R1.Set(A4,&waves,owner->WindowPos(0)+dth,dt);
		R2.Set(A4,&waves,owner->WindowPos(0)+dth,dt);
		for (auto strip : StripRange(A4,3,0,1,strongbool::no))
			A4(strip,dim[3]+1,3) = A4(strip,dim[3],3) - dx(3)*(A4.dfwd(strip,dim[3],1,1) + A4.dfwd(strip,dim[3],2,2));
	}

	A4.DownwardCopy(Rng(1,4),tw::grid::x,1);
	A4.UpwardCopy(Rng(1,4),tw::grid::x,1);
	A4.DownwardCopy(Rng(1,4),tw::grid::y,1);
	A4.UpwardCopy(Rng(1,4),tw::grid::y,1);

	// Swap A0 and A1 so A1 has most recent data

	A4.Swap(Rng(1,4),Rng(5,8));

	// COMPUTE E AND B FIELDS FOR PUSHER

	ComputeFinalFields();
}

void CoulombSolver::ComputeFinalFields()
{
	const tw::Float dti = this->dk(0);
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,0,1,false))
		{
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				F(v,k,0) = dti*(A4(v,k,1)-A4(v,k,5)) - A4.dbak(v,k,4,1);
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				F(v,k,1) = dti*(A4(v,k,2)-A4(v,k,6)) - A4.dbak(v,k,4,2);
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				F(v,k,2) = dti*(A4(v,k,3)-A4(v,k,7)) - A4.dbak(v,k,4,3);

			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				F(v,k,3) = 0.5*(A4.dbak(v,k,3,2) - A4.dbak(v,k,2,3) + A4.dbak(v,k,7,2) - A4.dbak(v,k,6,3));
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				F(v,k,4) = 0.5*(A4.dbak(v,k,1,3) - A4.dbak(v,k,3,1) + A4.dbak(v,k,5,3) - A4.dbak(v,k,7,1));
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				F(v,k,5) = 0.5*(A4.dbak(v,k,2,1) - A4.dbak(v,k,1,2) + A4.dbak(v,k,6,1) - A4.dbak(v,k,5,2));
		}
	}

	F.CopyFromNeighbors(All(F));
	F.ApplyBoundaryCondition(All(F));
}

void CoulombSolver::Report(Diagnostic& diagnostic)
{
	// Upon entry : A0(-1/2) , A1(1/2) , phi0(-1) , phi1(0)
	Electromagnetic::Report(diagnostic);

	diagnostic.ReportField("phi",A4,1,4,tw::dims::scalar_potential,"$\\phi$");

	std::map<tw::Int,std::string> xyz = {{1,"x"},{2,"y"},{3,"z"}};

	for (tw::Int ax=1;ax<=3;ax++)
	{
		scratch1 = 0.0;
		AddMulFieldData(scratch1,Rng(0),A4,Rng(ax),0.5);
		AddMulFieldData(scratch1,Rng(0),A4,Rng(ax+4),0.5);
		diagnostic.ReportField("A"+xyz[ax],scratch1,1,0,tw::dims::vector_potential,"$A_"+xyz[ax]+"$");

		scratch1 = 0.0;
		AddMulFieldData(scratch1,Rng(0),A4,Rng(ax),-dk(0));
		AddMulFieldData(scratch1,Rng(0),A4,Rng(ax+4),dk(0));
		diagnostic.ReportField("A"+xyz[ax]+"Dot",scratch1,1,0,tw::dims::electric_field,"$\\dot{A}_"+xyz[ax]+"$");
	}
}


///////////////////////
//  DIRECT SOLVER    //
//  ( with PML's )   //
///////////////////////


DirectSolver::DirectSolver(const std::string& name,Simulation* sim):Electromagnetic(name,sim)
{
	auto yee_name = owner->CreateTool("yee",tw::tool_type::yeePropagatorPML);
	yeeTool = std::dynamic_pointer_cast<YeePropagatorPML>(owner->UseTool(yee_name));
	enforceChargeConservation = true;
	for (tw::Int i=0;i<6;i++) layerThickness[i] = 0;
	for (tw::Int i=0;i<6;i++) reflectionCoefficient[i] = 0.01;
	A.Initialize(12,*this,owner);
	const tw::node5 xdims {1,dim[1],1,1,6};
	const tw::node5 ydims {1,dim[2],1,1,6};
	const tw::node5 zdims {1,dim[3],1,1,6};
	PMLx.Initialize(StaticSpace(xdims,sim->PhysicalSize(),std_packing,std_layers),owner); // sx,tx,jx,sxstar,txstar,jxstar
	PMLy.Initialize(StaticSpace(ydims,sim->PhysicalSize(),std_packing,std_layers),owner);
	PMLz.Initialize(StaticSpace(zdims,sim->PhysicalSize(),std_packing,std_layers),owner);

	#ifdef USE_OPENCL
	A.InitializeComputeBuffer();
	PMLx.InitializeComputeBuffer();
	PMLy.InitializeComputeBuffer();
	PMLz.InitializeComputeBuffer();
	yeeTool->SetupComputeKernels(F,A,PMLx,PMLy,PMLz,sources);
	#endif

	directives.Add("enforce charge conservation",new tw::input::Bool(&enforceChargeConservation),false);
	directives.Add("layers",new tw::input::Numbers<tw::Int>(layerThickness,6),false);
	directives.Add("reflection coefficient",new tw::input::Numbers<tw::Float>(reflectionCoefficient,6),false);
}

void DirectSolver::SetupPML(Field& pml,tw::Int g0,tw::Int gN,tw::Int L0,tw::Int L1,tw::Float R0,tw::Float R1,tw::Float ds)
{
	// works for uniform grid only
	// Each PML straddles a cell wall, occupying half of each adjacent cell.  Half of a ghost cell is occupied.
	// sigma known at cell centers for E_T, sigma-star at cell walls for B_T
	// therefore, sigma is specified at L+1 points, sigma-star at L points, given L layers

	// Example (L0=L1=2):

	//     o     o     o     o     o     o     o     o
	// ||  |  x     x  |  x     x     x  |  x     x  |  ||
	//     0     1     L               N1-L  N1-1    N1

	// Legend:  o=E_T, x=B_T , |=PML boundary, ||=outer ghost cell boundary

	tw::Int i,ig,L,N1=gN+1;
	tw::Float delta,maxConductivity,p,p0,sigma;
	const tw::Int bufferZone = 8;

	for (i=pml.LFG(1);i<=pml.UFG(1);i++)
	{
		pml(1,i,0,0,0) = 1.0;
		pml(1,i,0,0,1) = dx(0);
		pml(1,i,0,0,2) = 1.0;
		pml(1,i,0,0,3) = 1.0;
		pml(1,i,0,0,4) = dx(0);
		pml(1,i,0,0,5) = 1.0;
	}

	if (L0 && gN>1)
	{
		L = L0;
		delta = L * ds;
		maxConductivity = -1.5*std::log(R0)/delta;
		p0 = delta - 0.5*ds;
		for (i=pml.LFG(1);i<=pml.UFG(1);i++)
		{
			ig = g0+i;
			p = tw::Float(ig-1)*ds + 0.5*ds;
			if (ig >= 0 && ig <= L)
			{
				sigma = maxConductivity * (sqr((p - p0)/delta) + sqr(ds/delta)/12.0);
				sigma *= ig==0 || ig==L ? 0.5 : 1.0;
				pml(1,i,0,0,0) = std::exp(-sigma*dx(0));
				pml(1,i,0,0,1) = (1.0 - pml(1,i,0,0,0))/sigma;
				pml(1,i,0,0,2) = 0.0;
			}
			p = tw::Float(ig-1)*ds;
			if (ig > 0 && ig <= L)
			{
				sigma = maxConductivity * (sqr((p - p0)/delta) + sqr(ds/delta)/12.0);
				pml(1,i,0,0,3) = std::exp(-sigma*dx(0));
				pml(1,i,0,0,4) = (1.0 - pml(1,i,0,0,3))/sigma;
				pml(1,i,0,0,5) = 0.0;
			}
			if (ig > L && ig <= L+bufferZone)
			{
				pml(1,i,0,0,2) = QuinticRise(tw::Float(ig-L)/tw::Float(bufferZone));
				pml(1,i,0,0,5) = QuinticRise(tw::Float(ig-L-0.5)/tw::Float(bufferZone));
			}
		}
	}

	if (L1 && gN>1)
	{
		L = L1;
		delta = L * ds;
		maxConductivity = -1.5*std::log(R1)/delta;
		p0 = ds*N1 - delta - 0.5*ds;
		for (i=pml.LFG(1);i<=pml.UFG(1);i++)
		{
			ig = g0+i;
			p = tw::Float(ig-1)*ds + 0.5*ds;
			if (ig >= N1-L && ig <= N1)
			{
				sigma = maxConductivity * (sqr((p - p0)/delta) + sqr(ds/delta)/12.0);
				sigma *= ig==N1 || ig==N1-L ? 0.5 : 1.0;
				pml(1,i,0,0,0) = std::exp(-sigma*dx(0));
				pml(1,i,0,0,1) = (1.0 - pml(1,i,0,0,0))/sigma;
				pml(1,i,0,0,2) = 0.0;
			}
			p = tw::Float(ig-1)*ds;
			if (ig > N1-L && ig <= N1)
			{
				sigma = maxConductivity * (sqr((p - p0)/delta) + sqr(ds/delta)/12.0);
				pml(1,i,0,0,3) = std::exp(-sigma*dx(0));
				pml(1,i,0,0,4) = (1.0 - pml(1,i,0,0,3))/sigma;
				pml(1,i,0,0,5) = 0.0;
			}
			if (ig < N1-L && ig >= N1-L-bufferZone)
			{
				pml(1,i,0,0,2) = QuinticRise(tw::Float(N1-L-ig)/tw::Float(bufferZone));
				pml(1,i+1,0,0,5) = QuinticRise(tw::Float(N1-L-ig-0.5)/tw::Float(bufferZone));
			}
		}
	}
}

void DirectSolver::Initialize()
{
	logger::TRACE("initialize maxwell");

	Field A0;

	Electromagnetic::Initialize();
	A0.Initialize(3,*this,owner);

	SetExteriorBoundaryConditionsE(A,Rng(0,2),Rng(2,4),Rng(4,6));
	SetExteriorBoundaryConditionsB(A,Rng(6,8),Rng(8,10),Rng(10,12));

	// Setup PML media

	SetupPML(PMLx,owner->GlobalCellIndex(0,1),owner->GlobalDim(1),layerThickness[0],layerThickness[1],reflectionCoefficient[0],reflectionCoefficient[1],spacing[1]);
	SetupPML(PMLy,owner->GlobalCellIndex(0,2),owner->GlobalDim(2),layerThickness[2],layerThickness[3],reflectionCoefficient[2],reflectionCoefficient[3],spacing[2]);
	SetupPML(PMLz,owner->GlobalCellIndex(0,3),owner->GlobalDim(3),layerThickness[4],layerThickness[5],reflectionCoefficient[4],reflectionCoefficient[5],spacing[3]);

	#ifdef USE_OPENCL
	A.SendToComputeBuffer();
	F.SendToComputeBuffer();
	PMLx.SendToComputeBuffer();
	PMLy.SendToComputeBuffer();
	PMLz.SendToComputeBuffer();
	#endif

	// Initialize radiation fields

	LoadVectorPotential<0,1,2>(A0,-0.5*dx(0));
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,0,1,true))
		{
			#pragma omp simd
			for (tw::Int k=lfg[3];k<=ufg[3];k++)
			{
				A(v,k,0) = 0.5*dk(0)*A0(v,k,0);
				A(v,k,1) = 0.5*dk(0)*A0(v,k,0);
				A(v,k,2) = 0.5*dk(0)*A0(v,k,1);
				A(v,k,3) = 0.5*dk(0)*A0(v,k,1);
				A(v,k,4) = 0.5*dk(0)*A0(v,k,2);
				A(v,k,5) = 0.5*dk(0)*A0(v,k,2);
				A(v,k,6) = 0.0;
				A(v,k,7) = 0.0;
				A(v,k,8) = 0.0;
				A(v,k,9) = 0.0;
				A(v,k,10) = 0.0;
				A(v,k,11) = 0.0;
			}
		}
	}
	// Need B(-dt/2) to prep centered fields
	add_curlE<0,1,2,6,8,10>(1,A0,A,*owner,0.5);
	add_curlE<0,1,2,7,9,11>(1,A0,A,*owner,0.5);
	A.CopyFromNeighbors(All(A));
	A.ApplyBoundaryCondition(All(A));
	yeeTool->PrepCenteredFields(F,A,externalE,externalB);

	LoadVectorPotential<0,1,2>(A0,0.5*dx(0));
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,0,1,true))
		{
			#pragma omp simd
			for (tw::Int k=lfg[3];k<=ufg[3];k++)
			{
				A(v,k,0) -= 0.5*dk(0)*A0(v,k,0);
				A(v,k,1) -= 0.5*dk(0)*A0(v,k,0);
				A(v,k,2) -= 0.5*dk(0)*A0(v,k,1);
				A(v,k,3) -= 0.5*dk(0)*A0(v,k,1);
				A(v,k,4) -= 0.5*dk(0)*A0(v,k,2);
				A(v,k,5) -= 0.5*dk(0)*A0(v,k,2);
				A(v,k,6) = 0.0;
				A(v,k,7) = 0.0;
				A(v,k,8) = 0.0;
				A(v,k,9) = 0.0;
				A(v,k,10) = 0.0;
				A(v,k,11) = 0.0;
			}
		}
	}
	add_curlE<0,1,2,6,8,10>(1,A0,A,*owner,0.5);
	add_curlE<0,1,2,7,9,11>(1,A0,A,*owner,0.5);

	// Initialize static fields

	CopyFieldData(scratch1,Rng(0),sources,Rng(0));
	ellipticSolver->Solve(scratch2,scratch1,-1.0);
	add_grad<0,0,2,4>(1,scratch2,A,*owner,-0.5);
	add_grad<0,1,3,5>(1,scratch2,A,*owner,-0.5);

	// Clean up ghost cells

	A.CopyFromNeighbors(All(A));
	A.ApplyBoundaryCondition(All(A));
	yeeTool->CenteredFields(F,A);
}

void DirectSolver::MoveWindow()
{
	Electromagnetic::MoveWindow();
	for (auto s : StripRange(*this,3,0,1,strongbool::yes))
		A.Shift(All(A),s,-1,0.0);
	A.DownwardCopy(All(A),tw::grid::z,1);
}

void DirectSolver::ReadCheckpoint(std::ifstream& inFile)
{
	Electromagnetic::ReadCheckpoint(inFile);
	A.ReadCheckpoint(inFile);
	PMLx.ReadCheckpoint(inFile);
	PMLy.ReadCheckpoint(inFile);
	PMLz.ReadCheckpoint(inFile);
	inFile.read((char *)&enforceChargeConservation,sizeof(bool));
	inFile.read((char *)&layerThickness[0],sizeof(layerThickness));
	inFile.read((char *)&reflectionCoefficient[0],sizeof(reflectionCoefficient));
}

void DirectSolver::WriteCheckpoint(std::ofstream& outFile)
{
	Electromagnetic::WriteCheckpoint(outFile);
	A.WriteCheckpoint(outFile);
	PMLx.WriteCheckpoint(outFile);
	PMLy.WriteCheckpoint(outFile);
	PMLz.WriteCheckpoint(outFile);
	outFile.write((char *)&enforceChargeConservation,sizeof(bool));
	outFile.write((char *)&layerThickness[0],sizeof(layerThickness));
	outFile.write((char *)&reflectionCoefficient[0],sizeof(reflectionCoefficient));
}

void DirectSolver::Update()
{
	logger::TRACE("start Maxwell update");

	// Add electric antenna currents

	if (conductors.size())
	{
		#ifdef USE_OPENCL
		sources.ReceiveFromComputeBuffer();
		#endif
		for (auto c : conductors)
			if (c->currentType==EM::current::electric)
				c->DepositSources(sources, owner->WindowPos(0), dx(0));
		#ifdef USE_OPENCL
		sources.SendToComputeBuffer();
		#endif
	}

	Electromagnetic::Update();

	// Advance the electric field

	yeeTool->AdvanceE(A,PMLx,PMLy,PMLz,sources);
	if (conductors.size())
		yeeTool->UpdateInteriorBoundaryE(A,conductorMask);

	// Save the old magnetic field so final fields can be centered

	yeeTool->PrepCenteredFields(F,A,externalE,externalB);

	// Advance the magnetic field

	yeeTool->AdvanceB(A,PMLx,PMLy,PMLz);
	if (conductors.size())
		yeeTool->UpdateInteriorBoundaryB(A,conductorMask);

	// Setup the final field for the particle push

	yeeTool->CenteredFields(F,A);
}



///////////////////////////////////
//   CURVILINEAR DIRECT SOLVER   //
///////////////////////////////////


CurvilinearDirectSolver::CurvilinearDirectSolver(const std::string& name,Simulation* sim):DirectSolver(name,sim)
{
	A.Initialize(6,*this,owner);
}

void CurvilinearDirectSolver::Initialize()
{
	Field A0;
	const tw::Float dt = spacing[0];
	const tw::Float dti = freq[0];
	const tw::Float dth = 0.5*dt;

	Electromagnetic::Initialize();
	A0.Initialize(3,*this,owner);

	SetExteriorBoundaryConditionsE(A,Rng(0),Rng(1),Rng(2));
	SetExteriorBoundaryConditionsE(A,Rng(3),Rng(4),Rng(5));

	if (owner->gridGeometry==tw::grid::cylindrical)
	{
		A.SetBoundaryConditions(Rng(0),tw::grid::x,fld::none,fld::none);
		A.SetBoundaryConditions(Rng(1),tw::grid::x,fld::dirichletWall,fld::dirichletCell);
		A.SetBoundaryConditions(Rng(2),tw::grid::x,fld::neumannWall,fld::dirichletCell);
		A.SetBoundaryConditions(Rng(3),tw::grid::x,fld::dirichletWall,fld::dirichletCell);
		A.SetBoundaryConditions(Rng(4,6),tw::grid::x,fld::none,fld::none);
	}

	// Initialize radiation fields

	LoadVectorPotential<0,1,2>(A0,-dth);
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,0,1,true))
		{
			#pragma omp simd
			for (tw::Int k=lfg[3];k<=ufg[3];k++)
			{
				A(v,k,0) = dti*A0(v,k,0);
				A(v,k,1) = dti*A0(v,k,1);
				A(v,k,2) = dti*A0(v,k,2);
			}
		}
	}
	// Save the B field at -dt/2 for computing centered fields
	add_curlE<0,1,2,3,4,5>(1,A0,F,*owner,1.0);
	F.CopyFromNeighbors(All(F));
	F.ApplyBoundaryCondition(All(F));

	LoadVectorPotential<0,1,2>(A0,dth);
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,0,1,true))
		{
			#pragma omp simd
			for (tw::Int k=lfg[3];k<=ufg[3];k++)
			{
				A(v,k,0) -= dti*A0(v,k,0);
				A(v,k,1) -= dti*A0(v,k,1);
				A(v,k,2) -= dti*A0(v,k,2);
			}
		}
	}
	add_curlE<0,1,2,3,4,5>(1,A0,A,*owner,1.0);

	// Initialize static fields

	CopyFieldData(scratch1,Rng(0),sources,Rng(0));
	ellipticSolver->Solve(scratch2,scratch1,-1.0);
	add_grad<0,0,1,2>(1,scratch2,A,*owner,-1.0);

	// Clean up ghost cells

	SetSingularPointsE();
	SetSingularPointsB();
	A.CopyFromNeighbors(All(A));
	A.ApplyBoundaryCondition(All(A));

	// Time centered fields for particle pusher

	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,0,1,true))
		{
			#pragma omp simd
			for (tw::Int k=lfg[3];k<=ufg[3];k++)
			{
				F(v,k,0) = A(v,k,0);
				F(v,k,1) = A(v,k,1);
				F(v,k,2) = A(v,k,2);
				F(v,k,3) = 0.5*(A(v,k,3) + F(v,k,3));
				F(v,k,4) = 0.5*(A(v,k,4) + F(v,k,4));
				F(v,k,5) = 0.5*(A(v,k,5) + F(v,k,5));
			}
		}
	}
}

void CurvilinearDirectSolver::SetSingularPointsE()
{
	tw::Int i,j,k;
	if (owner->gridGeometry==tw::grid::cylindrical && owner->X(0,1)<0.0)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (j=lfg[2];j<=ufg[2];j++)
				A(1,1,j,k,0) = 0.0;
	}
	if (owner->gridGeometry==tw::grid::spherical && owner->X(0,1)<0.0)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (j=lfg[2];j<=ufg[2];j++)
				A(1,1,j,k,0) = 0.0;
	}
	if (owner->gridGeometry==tw::grid::spherical && owner->X(0,2)<0.0)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (i=lfg[1];i<=ufg[1];i++)
				A(1,i,1,k,1) = 0.0;
	}
	if (owner->gridGeometry==tw::grid::spherical && owner->X(ufg[2],2)>pi)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (i=lfg[1];i<=ufg[1];i++)
				A(1,i,ufg[2],k,1) = 0.0;
	}
}

void CurvilinearDirectSolver::SetSingularPointsB()
{
	// Assume that F holds the old staggered B-field at this point
	// Then use curlE = -dB/dt to update B
	tw::Int i,j,k;
	const tw::Float dt = dx(0);
	if (owner->gridGeometry==tw::grid::cylindrical && owner->X(0,1)<0.0)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (j=lfg[2];j<=ufg[2];j++)
			{
				A(1,1,j,k,4) = 0.0;
				A(1,1,j,k,5) = F(1,1,j,k,5) - 2.0*dt*A(1,1,j,k,1)/owner->X(1,1);
			}
	}
	if (owner->gridGeometry==tw::grid::spherical && owner->X(0,2)<0.0)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (i=lfg[1];i<=ufg[1];i++)
			{
				A(1,i,1,k,5) = 0.0;
				A(1,i,1,k,3) = F(1,i,1,k,3) - 2.0*dt*A(1,i,1,k,2)/(owner->X(i,1)*std::sin(owner->X(1,2)));
			}
	}
	if (owner->gridGeometry==tw::grid::spherical && owner->X(ufg[2],2)>pi)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (i=lfg[1];i<=ufg[1];i++)
			{
				A(1,i,ufg[2],k,5) = 0.0;
				A(1,i,ufg[2],k,3) = F(1,i,ufg[2],k,3) + 2.0*dt*A(1,i,dim[2],k,2)/(owner->X(i,1)*std::sin(owner->X(dim[2],2)));
			}
	}
}

void CurvilinearDirectSolver::Update()
{
	tw::vec3 S;
	const tw::Float dt = dx(0);

	// Add electric antenna currents

	if (conductors.size())
	{
		#ifdef USE_OPENCL
		sources.ReceiveFromComputeBuffer();
		#endif
		for (auto c : conductors)
			if (c->currentType==EM::current::electric)
				c->DepositSources(sources, owner->WindowPos(0), dt);
		#ifdef USE_OPENCL
		sources.SendToComputeBuffer();
		#endif
	}

	Electromagnetic::Update();

	// Advance the electric field

	add_curlB<3,4,5,0,1,2>(1,A,A,*owner,dt);
	AddMulFieldData(A,Rng(0,3),sources,Rng(1,4),-dt);
	SetSingularPointsE();
	A.CopyFromNeighbors(Rng(0,3));

	// Save the old magnetic field so final fields can be centered

	CopyFieldData(F,Rng(3,6),A,Rng(3,6));

	// Advance the magnetic field

	add_curlE<0,1,2,3,4,5>(1,A,A,*owner,-dt);
	SetSingularPointsB();
	A.UpwardCopy(Rng(3,6),tw::grid::x,1);
	A.UpwardCopy(Rng(3,6),tw::grid::y,1);
	A.UpwardCopy(Rng(3,6),tw::grid::z,1);

	// Setup the final field for the particle push

	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,0,1,true))
		{
			#pragma omp simd
			for (tw::Int k=lfg[3];k<=ufg[3];k++)
			{
				F(v,k,0) = A(v,k,0);
				F(v,k,1) = A(v,k,1);
				F(v,k,2) = A(v,k,2);
				F(v,k,3) = 0.5*(A(v,k,3) + F(v,k,3));
				F(v,k,4) = 0.5*(A(v,k,4) + F(v,k,4));
				F(v,k,5) = 0.5*(A(v,k,5) + F(v,k,5));
			}
		}
	}
}


FarFieldDiagnostic::FarFieldDiagnostic(const std::string& name,Simulation *sim) : Module(name,sim)
{
	J4 = NULL;
	directives.Add("radius",new tw::input::Float(&radius));
	directives.Add("bounds",new tw::input::Numbers<tw::Float>(&bounds[0],6));
	directives.Add("dims",new tw::input::Numbers<tw::Int>(&dims[1],3));
	directives.Add("period",new tw::input::Int(&period));
}

void FarFieldDiagnostic::Initialize()
{
	if (J4==NULL)
		throw tw::FatalError("Far field diagnostic did not find a source field.");
	tw::vec4 size(dx(0),bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]);
	A.Initialize(StaticSpace(tw::node5{1,dims[1],dims[2],dims[3],3},size,std_packing,std_layers),owner);
}

bool FarFieldDiagnostic::InspectResource(void *resource,const std::string& description)
{
	if (description=="electromagnetic:sources")
	{
		J4 = (Field*)resource;
		return true;
	}
	return false;
}

void FarFieldDiagnostic::ReadCheckpoint(std::ifstream& inFile)
{
	Module::ReadCheckpoint(inFile);
	A.ReadCheckpoint(inFile);
}

void FarFieldDiagnostic::WriteCheckpoint(std::ofstream& outFile)
{
	Module::WriteCheckpoint(outFile);
	A.WriteCheckpoint(outFile);
}

void FarFieldDiagnostic::Update()
{
	// Integrate over a plane in the source region for each interrogation cell in the far field
	const tw::Float zmin = owner->X(1,3) - 0.5*owner->dX(1,3);
	const tw::Float zmax = owner->X(Dim(3),3) + 0.5*owner->dX(Dim(3),3);
	const tw::Float tp = owner->WindowPos(0);
	const tw::Float dtau = dx(0) * tw::Float(period);

	if (owner->StepNow() % period==0)
		for (auto farCell : InteriorCellRange(A,1))
		{
			const tw::Float tNow = bounds[0] + A.dx(1)*tw::Float(farCell.dcd1());
			const tw::Float thetaNow = bounds[2] + A.dx(2)*tw::Float(farCell.dcd2());
			const tw::Float phiNow = bounds[4] + A.dx(3)*tw::Float(farCell.dcd3());
			const tw::vec3 n(std::sin(thetaNow)*std::cos(phiNow),std::sin(thetaNow)*std::sin(phiNow),std::cos(thetaNow));
			#pragma omp parallel
			{
				for (auto s : StripRange(*this,3,0,1,strongbool::no))
				{
					std::valarray<tw::Float> j4(4);
					tw::vec3 rp = owner->Pos(s,1);
					rp.z = ((tp - tNow) - rp.x*n.x - rp.y*n.y)/n.z;
					const tw::Float dS = owner->dS(s,1,3)/n.z;
					if (rp.z > zmin && rp.z < zmax)
					{
						weights_3D w;
						owner->GetWeights(&w,tw::vec4(tp,rp));
						J4->Interpolate(All(*J4),j4,w);
						A.Pack(farCell,A(farCell) + tw::vec3(j4[1],j4[2],j4[3]) * dS * dtau / radius);
					}
				}
			}
		}
}

void FarFieldDiagnostic::FinalReport()
{
	std::string fileName;
	tw::Int master = 0;
	tw::Int curr = owner->strip[0].Get_rank();

	// Each node has different contributions to the same far field grid.
	// We merely need to form the superposition of all these grids.
	Vec3Field accum;
	accum.Initialize(A,owner);
	owner->strip[0].Sum(&A(1,0,0,0,0),&accum(1,0,0,0,0),sizeof(tw::vec3)*accum.TotalCells(),master);

	if (curr==master)
	{
		// put vector potential in coulomb gauge and spherical coordinates
		for (auto farCell : InteriorCellRange(accum,1))
		{
			const tw::Float theta = bounds[2] + A.dx(2)*tw::Float(farCell.dcd2());
			const tw::Float phi = bounds[4] + A.dx(3)*tw::Float(farCell.dcd3());
			const tw::vec3 nr( std::sin(theta)*std::cos(phi) , std::sin(theta)*std::sin(phi) , std::cos(theta) );
			const tw::vec3 nq( std::cos(theta)*std::cos(phi) , std::cos(theta)*std::sin(phi) , -std::sin(theta) );
			const tw::vec3 nf( -std::sin(phi) , std::cos(phi) , 0.0 );
			const tw::vec3 ACG = nr | (accum(farCell) | nr); // form coulomb gauge A
			accum.Pack(farCell, tw::vec3( ACG^nr , ACG^nq , ACG^nf )); // put in spherical coordinates
		}

		npy_writer writer;
		std::valarray<float> gData(dims[1]*dims[2]*dims[3]);

		fileName = name + "-Atheta.npy";
		writer.write_header(fileName,dims);
		for (tw::Int i=1;i<=dims[1];i++)
			for (tw::Int j=1;j<=dims[2];j++)
				for (tw::Int k=1;k<=dims[3];k++)
					gData[(i-1)*dims[2]*dims[3] + (j-1)*dims[3] + (k-1)] = accum(i,j,k).y;
		writer.add_frame(fileName,&gData[0],dims);

		fileName = name + "-Aphi.npy";
		writer.write_header(fileName,dims);
		for (tw::Int i=1;i<=dims[1];i++)
			for (tw::Int j=1;j<=dims[2];j++)
				for (tw::Int k=1;k<=dims[3];k++)
					gData[(i-1)*dims[2]*dims[3] + (j-1)*dims[3] + (k-1)] = accum(i,j,k).z;
		writer.add_frame(fileName,&gData[0],dims);
	}
}
