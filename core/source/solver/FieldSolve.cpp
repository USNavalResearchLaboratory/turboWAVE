#include "simulation.h"
#include "fieldSolve.h"
using namespace tw::bc;

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
	ellipticSolver = NULL;
}

FieldSolver::~FieldSolver()
{
	if (ellipticSolver!=NULL)
		owner->RemoveTool(ellipticSolver);
}

void FieldSolver::VerifyInput()
{
	Module::VerifyInput();
	// Find an elliptic solver on the list of tools associated with this module
	for (auto tool : moduleTool)
	{
		ellipticSolver = dynamic_cast<EllipticSolver*>(tool);
		if (ellipticSolver!=NULL)
			break;
	}
	// If the tool could not be found, create one automatically
	if (ellipticSolver==NULL)
		ellipticSolver = (EllipticSolver*)owner->CreateTool("default_facr_poisson_solver",tw::tool_type::facrPoissonSolver);
}


////////////////////////////////
// ELECTROMAGNETIC BASE CLASS //
////////////////////////////////


Electromagnetic::Electromagnetic(const std::string& name,Simulation* sim):FieldSolver(name,sim)
{
	if (1.0/dt <= sqrt(
		(sim->globalCells[1]==1 ? 0.0 : 1.0)/sqr(dx(*sim)) +
		(sim->globalCells[2]==1 ? 0.0 : 1.0)/sqr(dy(*sim)) +
		(sim->globalCells[3]==1 ? 0.0 : 1.0)/sqr(dz(*sim))))
		throw tw::FatalError("Courant condition is violated.");

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
	directives.Add("gamma beam",new tw::input::Float(&gammaBeam),false);
}

void Electromagnetic::ExchangeResources()
{
	PublishResource(&F,"electromagnetic:F");
	PublishResource(&sources,"electromagnetic:sources");
}

void Electromagnetic::Initialize()
{
	FieldSolver::Initialize();
	InitializeConductors();

	SetExteriorBoundaryConditionsE(F,Element(0),Element(1),Element(2));
	SetExteriorBoundaryConditionsB(F,Element(3),Element(4),Element(5));

	sources.SetBoundaryConditions(tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	sources.SetBoundaryConditions(Element(1),tw::grid::x,fld::normalFluxFixed,fld::normalFluxFixed);
	if (owner->gridGeometry==tw::grid::cylindrical)
		sources.SetBoundaryConditions(Element(0),tw::grid::x,fld::neumannWall,fld::dirichletCell);

	sources.SetBoundaryConditions(tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	sources.SetBoundaryConditions(Element(2),tw::grid::y,fld::normalFluxFixed,fld::normalFluxFixed);

	sources.SetBoundaryConditions(tw::grid::z,fld::dirichletCell,fld::dirichletCell);
	sources.SetBoundaryConditions(Element(3),tw::grid::z,fld::normalFluxFixed,fld::normalFluxFixed);

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
	sources.DepositFromNeighbors();
	sources.ApplyFoldingCondition(); // only apply before conversion to density
	conserved_current_to_dens<0,1,2,3>(sources,*owner);
	sources.ApplyBoundaryCondition(); // only apply after conversion to density
	sources.Smooth(*owner,smoothing,compensation);
}
#endif

void Electromagnetic::MoveWindow()
{
	Module::MoveWindow();
	for (auto s : StripRange(*this,3,strongbool::yes))
		F.Shift(s,-1,0.0);
	F.DownwardCopy(tw::grid::z,1);
}

void Electromagnetic::InitializeConductors()
{
	// Conductor resides in cells shifted back by 1/2
	// These have E known along upper edges and B normal to upper walls
	// The conductor fills the whole cell or none of the cell
	tw::vec3 shiftedCenter;
	//#pragma omp parallel for private(i,j,k,s,shiftedCenter) collapse(3) schedule(static)
	for (auto cell : EntireCellRange(*this))
	{
		conductorMask(cell) = 1.0;
		shiftedCenter = owner->Pos(cell) - 0.5*owner->dPos(cell);
		for (auto c : conductor)
			if (c->affectsA && c->theRgn->Inside(shiftedCenter,*owner))
				conductorMask(cell) = 0.0;
	}
}

void Electromagnetic::SetExteriorBoundaryConditionsE(Field& A,const Element& ex,const Element& ey,const Element& ez)
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

void Electromagnetic::SetExteriorBoundaryConditionsB(Field& A,const Element& bx,const Element& by,const Element& bz)
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
	tw::Float beta = sqrt(1.0 - 1.0/sqr(gammaBeam));
	tw::Float w = beta*dt*freq.z;
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
				tw::xstrip<3> v(*this,i,j,0);
				for (k=1;k<=dim[3];k++)
					rhs(v,k,0) = DtPhi(v,k,0,ax) - sources(v,k,ax);
			}
		rhs.CopyFromNeighbors();

		ellipticSolver->Solve(ans,rhs,1.0);
		for (auto cell : EntireCellRange(*this))
			A4(cell,ax+4) = ans(cell);

		// Find Solution at last time

		for (i=1;i<=dim[1];i++)
			for (j=1;j<=dim[2];j++)
			{
				tw::xstrip<3> v(*this,i,j,0);
				for (k=1;k<=dim[3];k++)
					rhs(v,k,0) = DtPhi(v,k,0,ax) - sources(v,k,ax);
			}
		rhs.CopyFromNeighbors();

		for (k=lfg[3];k<=dim[3];k++)
			for (j=lfg[2];j<=ufg[2];j++)
				for (i=lfg[1];i<=ufg[1];i++)
					rhs(i,j,k) = (1.0 - w)*rhs(i,j,k) + w*rhs(i,j,k+1);
		rhs.DownwardCopy(tw::grid::z,1);

		ellipticSolver->Solve(ans,rhs,1.0);
		for (auto cell : EntireCellRange(*this))
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

	for (auto cell : InteriorCellRange(*this))
		temp(cell) = 0.5*(Norm(F.Vec3(cell,0))+Norm(F.Vec3(cell,3)));
	diagnostic.VolumeIntegral("FieldEnergy",temp,0);

	diagnostic.VolumeIntegral("TotalCharge",sources,0);
	diagnostic.FirstMoment("Dx",sources,0,dipoleCenter,tw::grid::x);
	diagnostic.FirstMoment("Dy",sources,0,dipoleCenter,tw::grid::y);
	diagnostic.FirstMoment("Dz",sources,0,dipoleCenter,tw::grid::z);
	diagnostic.Field("Ex",F,0,tw::dims::electric_field,"$E_x$");
	diagnostic.Field("Ey",F,1,tw::dims::electric_field,"$E_y$");
	diagnostic.Field("Ez",F,2,tw::dims::electric_field,"$E_z$");
	diagnostic.Field("Bx",F,3,tw::dims::magnetic_field,"$B_x$");
	diagnostic.Field("By",F,4,tw::dims::magnetic_field,"$B_y$");
	diagnostic.Field("Bz",F,5,tw::dims::magnetic_field,"$B_z$");
	diagnostic.Field("rho",sources,0,tw::dims::charge_density,"$\\rho$");
	diagnostic.Field("Jx",sources,1,tw::dims::current_density,"$j_x$");
	diagnostic.Field("Jy",sources,2,tw::dims::current_density,"$j_y$");
	diagnostic.Field("Jz",sources,3,tw::dims::current_density,"$j_z$");
}


////////////////////////////////////////////////
// COULOMB GAUGE FIELD SOLVER (original WAVE) //
////////////////////////////////////////////////


CoulombSolver::CoulombSolver(const std::string& name,Simulation* sim):Electromagnetic(name,sim)
{
	A4.Initialize(8,*this,owner);

	L1.Initialize(sim,sim,&wave,tw::grid::z,tw::grid::low,1,5);
	L2.Initialize(sim,sim,&wave,tw::grid::z,tw::grid::low,2,6);
	R1.Initialize(sim,sim,&wave,tw::grid::z,tw::grid::high,1,5);
	R2.Initialize(sim,sim,&wave,tw::grid::z,tw::grid::high,2,6);
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
	ellipticSolver->SetFieldsBoundaryConditions(scratch2,Element(0));

	SetExteriorBoundaryConditionsE(A4,Element(1),Element(2),Element(3));
	SetExteriorBoundaryConditionsE(A4,Element(5),Element(6),Element(7));

	// Initialize radiation fields

	LoadVectorPotential<1,2,3>(A4,-dth);
	LoadVectorPotential<5,6,7>(A4,dth);

	// Initialize static fields

	CopyFieldData(scratch1,Element(0),sources,Element(0));
	ellipticSolver->Solve(scratch2,scratch1,-1.0);
	CopyFieldData(A4,Element(0),scratch2,Element(0));
	CopyFieldData(A4,Element(4),scratch2,Element(0));

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

	// upon entry : A0(-1/2) , A1(1/2) , phi0(-1) , phi1(0) , J(1/2) , rho(1)
	// upon exit : A0(1/2) , A1(3/2) , phi0(0) , phi1(1)

	// ADVANCE THE POTENTIALS

	// Update scalar potential and compute dphi/dt

	CopyFieldData(A4,Element(0),A4,Element(4));
	CopyFieldData(scratch1,Element(0),sources,Element(0));
	ellipticSolver->Solve(scratch2,scratch1,-1.0);
	CopyFieldData(A4,Element(4),scratch2,Element(0));
	AddMulFieldData(scratch2,Element(0),A4,Element(0),-1.0);
	scratch2 *= dti;

	// Deal with beam initialization

	//if (gammaBeam!=1.0 && owner->elapsedTime==0.0)
	//	ForceQuasistaticVectorPotential(A4,scratch2);

	// Must update boundary memory before vector potential

	L1.UpdateBoundaryMemory(A4,dt);
	L2.UpdateBoundaryMemory(A4,dt);
	R1.UpdateBoundaryMemory(A4,dt);
	R2.UpdateBoundaryMemory(A4,dt);

	// Update the interior : do this in place by putting new data in A0

	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,false))
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

	A4.ApplyBoundaryCondition(Element(1,3));
	A4.CopyFromNeighbors(Element(1,3));

	if (owner->n0[3]==MPI_PROC_NULL)
	{
		L1.Set(A4,owner->elapsedTime+dth,dt);
		L2.Set(A4,owner->elapsedTime+dth,dt);
		for (auto strip : StripRange(A4,3,strongbool::no))
			A4(strip,0,3) = A4(strip,1,3) + spacing.z*(A4.dfwd(strip,0,1,1) + A4.dfwd(strip,0,2,2));
	}

	if (owner->n1[3]==MPI_PROC_NULL)
	{
		R1.Set(A4,owner->elapsedTime+dth,dt);
		R2.Set(A4,owner->elapsedTime+dth,dt);
		for (auto strip : StripRange(A4,3,strongbool::no))
			A4(strip,dim[3]+1,3) = A4(strip,dim[3],3) - spacing.z*(A4.dfwd(strip,dim[3],1,1) + A4.dfwd(strip,dim[3],2,2));
	}

	A4.DownwardCopy(tw::grid::x,Element(1,3),1);
	A4.UpwardCopy(tw::grid::x,Element(1,3),1);
	A4.DownwardCopy(tw::grid::y,Element(1,3),1);
	A4.UpwardCopy(tw::grid::y,Element(1,3),1);

	// Swap A0 and A1 so A1 has most recent data

	A4.Swap(Element(1,3),Element(5,7));

	// COMPUTE E AND B FIELDS FOR PUSHER

	ComputeFinalFields();
}

void CoulombSolver::ComputeFinalFields()
{
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,false))
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

	F.CopyFromNeighbors();
	F.ApplyBoundaryCondition();
}

void CoulombSolver::Report(Diagnostic& diagnostic)
{
	// Upon entry : A0(-1/2) , A1(1/2) , phi0(-1) , phi1(0)
	Electromagnetic::Report(diagnostic);

	diagnostic.Field("phi",A4,4,tw::dims::scalar_potential,"$\\phi$");

	std::map<tw::Int,std::string> xyz = {{1,"x"},{2,"y"},{3,"z"}};

	for (tw::Int ax=1;ax<=3;ax++)
	{
		scratch1 = 0.0;
		AddMulFieldData(scratch1,Element(0),A4,Element(ax),0.5);
		AddMulFieldData(scratch1,Element(0),A4,Element(ax+4),0.5);
		diagnostic.Field("A"+xyz[ax],scratch1,0,tw::dims::vector_potential,"$A_"+xyz[ax]+"$");

		scratch1 = 0.0;
		AddMulFieldData(scratch1,Element(0),A4,Element(ax),-1.0/dt);
		AddMulFieldData(scratch1,Element(0),A4,Element(ax+4),1.0/dt);
		diagnostic.Field("A"+xyz[ax]+"Dot",scratch1,0,tw::dims::electric_field,"$\\dot{A}_"+xyz[ax]+"$");
	}
}


///////////////////////
//  DIRECT SOLVER    //
//  ( with PML's )   //
///////////////////////


DirectSolver::DirectSolver(const std::string& name,Simulation* sim):Electromagnetic(name,sim)
{
	yeeTool = (YeePropagatorPML*)owner->CreateTool("yee",tw::tool_type::yeePropagatorPML);
	enforceChargeConservation = true;
	for (tw::Int i=0;i<6;i++) layerThickness[i] = 0;
	for (tw::Int i=0;i<6;i++) reflectionCoefficient[i] = 0.01;
	A.Initialize(12,*this,owner);
	PMLx.Initialize(6,DiscreteSpace(dim[1],1,1,corner,size),owner); // sx,tx,jx,sxstar,txstar,jxstar
	PMLy.Initialize(6,DiscreteSpace(dim[2],1,1,corner,size),owner);
	PMLz.Initialize(6,DiscreteSpace(dim[3],1,1,corner,size),owner);

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

DirectSolver::~DirectSolver()
{
	owner->RemoveTool(yeeTool);
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
		pml(i,0,0,0) = 1.0;
		pml(i,0,0,1) = dt;
		pml(i,0,0,2) = 1.0;
		pml(i,0,0,3) = 1.0;
		pml(i,0,0,4) = dt;
		pml(i,0,0,5) = 1.0;
	}

	if (L0 && gN>1)
	{
		L = L0;
		delta = L * ds;
		maxConductivity = -1.5*log(R0)/delta;
		p0 = delta - 0.5*ds;
		for (i=pml.LFG(1);i<=pml.UFG(1);i++)
		{
			ig = g0+i;
			p = tw::Float(ig-1)*ds + 0.5*ds;
			if (ig >= 0 && ig <= L)
			{
				sigma = maxConductivity * (sqr((p - p0)/delta) + sqr(ds/delta)/12.0);
				sigma *= ig==0 || ig==L ? 0.5 : 1.0;
				pml(i,0,0,0) = exp(-sigma*dt);
				pml(i,0,0,1) = (1.0 - pml(i,0,0,0))/sigma;
				pml(i,0,0,2) = 0.0;
			}
			p = tw::Float(ig-1)*ds;
			if (ig > 0 && ig <= L)
			{
				sigma = maxConductivity * (sqr((p - p0)/delta) + sqr(ds/delta)/12.0);
				pml(i,0,0,3) = exp(-sigma*dt);
				pml(i,0,0,4) = (1.0 - pml(i,0,0,3))/sigma;
				pml(i,0,0,5) = 0.0;
			}
			if (ig > L && ig <= L+bufferZone)
			{
				pml(i,0,0,2) = QuinticRise(tw::Float(ig-L)/tw::Float(bufferZone));
				pml(i,0,0,5) = QuinticRise(tw::Float(ig-L-0.5)/tw::Float(bufferZone));
			}
		}
	}

	if (L1 && gN>1)
	{
		L = L1;
		delta = L * ds;
		maxConductivity = -1.5*log(R1)/delta;
		p0 = ds*N1 - delta - 0.5*ds;
		for (i=pml.LFG(1);i<=pml.UFG(1);i++)
		{
			ig = g0+i;
			p = tw::Float(ig-1)*ds + 0.5*ds;
			if (ig >= N1-L && ig <= N1)
			{
				sigma = maxConductivity * (sqr((p - p0)/delta) + sqr(ds/delta)/12.0);
				sigma *= ig==N1 || ig==N1-L ? 0.5 : 1.0;
				pml(i,0,0,0) = exp(-sigma*dt);
				pml(i,0,0,1) = (1.0 - pml(i,0,0,0))/sigma;
				pml(i,0,0,2) = 0.0;
			}
			p = tw::Float(ig-1)*ds;
			if (ig > N1-L && ig <= N1)
			{
				sigma = maxConductivity * (sqr((p - p0)/delta) + sqr(ds/delta)/12.0);
				pml(i,0,0,3) = exp(-sigma*dt);
				pml(i,0,0,4) = (1.0 - pml(i,0,0,3))/sigma;
				pml(i,0,0,5) = 0.0;
			}
			if (ig < N1-L && ig >= N1-L-bufferZone)
			{
				pml(i,0,0,2) = QuinticRise(tw::Float(N1-L-ig)/tw::Float(bufferZone));
				pml(i+1,0,0,5) = QuinticRise(tw::Float(N1-L-ig-0.5)/tw::Float(bufferZone));
			}
		}
	}
}

void DirectSolver::Initialize()
{
	Field A0;

	Electromagnetic::Initialize();
	A0.Initialize(3,*this,owner);

	SetExteriorBoundaryConditionsE(A,Element(0,1),Element(2,3),Element(4,5));
	SetExteriorBoundaryConditionsB(A,Element(6,7),Element(8,9),Element(10,11));

	// Setup PML media

	SetupPML(PMLx,owner->GlobalCellIndex(0,1),owner->globalCells[1],layerThickness[0],layerThickness[1],reflectionCoefficient[0],reflectionCoefficient[1],spacing.x);
	SetupPML(PMLy,owner->GlobalCellIndex(0,2),owner->globalCells[2],layerThickness[2],layerThickness[3],reflectionCoefficient[2],reflectionCoefficient[3],spacing.y);
	SetupPML(PMLz,owner->GlobalCellIndex(0,3),owner->globalCells[3],layerThickness[4],layerThickness[5],reflectionCoefficient[4],reflectionCoefficient[5],spacing.z);

	#ifdef USE_OPENCL
	A.SendToComputeBuffer();
	F.SendToComputeBuffer();
	PMLx.SendToComputeBuffer();
	PMLy.SendToComputeBuffer();
	PMLz.SendToComputeBuffer();
	#endif

	// Initialize radiation fields

	LoadVectorPotential<0,1,2>(A0,-dth);
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,true))
		{
			#pragma omp simd
			for (tw::Int k=lfg[3];k<=ufg[3];k++)
			{
				A(v,k,0) = 0.5*dti*A0(v,k,0);
				A(v,k,1) = 0.5*dti*A0(v,k,0);
				A(v,k,2) = 0.5*dti*A0(v,k,1);
				A(v,k,3) = 0.5*dti*A0(v,k,1);
				A(v,k,4) = 0.5*dti*A0(v,k,2);
				A(v,k,5) = 0.5*dti*A0(v,k,2);
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
	add_curlE<0,1,2,6,8,10>(A0,A,*owner,0.5);
	add_curlE<0,1,2,7,9,11>(A0,A,*owner,0.5);
	A.CopyFromNeighbors();
	A.ApplyBoundaryCondition();
	yeeTool->PrepCenteredFields(F,A);

	LoadVectorPotential<0,1,2>(A0,dth);
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,true))
		{
			#pragma omp simd
			for (tw::Int k=lfg[3];k<=ufg[3];k++)
			{
				A(v,k,0) -= 0.5*dti*A0(v,k,0);
				A(v,k,1) -= 0.5*dti*A0(v,k,0);
				A(v,k,2) -= 0.5*dti*A0(v,k,1);
				A(v,k,3) -= 0.5*dti*A0(v,k,1);
				A(v,k,4) -= 0.5*dti*A0(v,k,2);
				A(v,k,5) -= 0.5*dti*A0(v,k,2);
				A(v,k,6) = 0.0;
				A(v,k,7) = 0.0;
				A(v,k,8) = 0.0;
				A(v,k,9) = 0.0;
				A(v,k,10) = 0.0;
				A(v,k,11) = 0.0;
			}
		}
	}
	add_curlE<0,1,2,6,8,10>(A0,A,*owner,0.5);
	add_curlE<0,1,2,7,9,11>(A0,A,*owner,0.5);

	// Initialize static fields

	CopyFieldData(scratch1,Element(0),sources,Element(0));
	ellipticSolver->Solve(scratch2,scratch1,-1.0);
	add_grad<0,0,2,4>(scratch2,A,*owner,-0.5);
	add_grad<0,1,3,5>(scratch2,A,*owner,-0.5);

	// Clean up ghost cells

	A.CopyFromNeighbors();
	A.ApplyBoundaryCondition();
	yeeTool->CenteredFields(F,A);
}

void DirectSolver::MoveWindow()
{
	Electromagnetic::MoveWindow();
	for (auto s : StripRange(*this,3,strongbool::yes))
		A.Shift(s,-1,0.0);
	A.DownwardCopy(tw::grid::z,1);
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
	// Add electric antenna currents

	if (conductor.size())
	{
		#ifdef USE_OPENCL
		sources.ReceiveFromComputeBuffer();
		#endif
		for (auto c : conductor)
			if (c->currentType==EM::current::electric)
				c->DepositSources(sources, owner->elapsedTime, dt);
		#ifdef USE_OPENCL
		sources.SendToComputeBuffer();
		#endif
	}

	Electromagnetic::Update();

	// Advance the electric field

	yeeTool->AdvanceE(A,PMLx,PMLy,PMLz,sources);
	if (conductor.size())
		yeeTool->UpdateInteriorBoundaryE(A,conductorMask);

	// Save the old magnetic field so final fields can be centered

	yeeTool->PrepCenteredFields(F,A);

	// Advance the magnetic field

	yeeTool->AdvanceB(A,PMLx,PMLy,PMLz);
	if (conductor.size())
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

	Electromagnetic::Initialize();
	A0.Initialize(3,*this,owner);

	SetExteriorBoundaryConditionsE(A,Element(0),Element(1),Element(2));
	SetExteriorBoundaryConditionsE(A,Element(3),Element(4),Element(5));

	if (owner->gridGeometry==tw::grid::cylindrical)
	{
		A.SetBoundaryConditions(Element(0),tw::grid::x,fld::none,fld::none);
		A.SetBoundaryConditions(Element(1),tw::grid::x,fld::dirichletWall,fld::dirichletCell);
		A.SetBoundaryConditions(Element(2),tw::grid::x,fld::neumannWall,fld::dirichletCell);
		A.SetBoundaryConditions(Element(3),tw::grid::x,fld::dirichletWall,fld::dirichletCell);
		A.SetBoundaryConditions(Element(4,5),tw::grid::x,fld::none,fld::none);
	}

	// Initialize radiation fields

	LoadVectorPotential<0,1,2>(A0,-dth);
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,true))
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
	add_curlE<0,1,2,3,4,5>(A0,F,*owner,1.0);
	F.CopyFromNeighbors();
	F.ApplyBoundaryCondition();

	LoadVectorPotential<0,1,2>(A0,dth);
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,true))
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
	add_curlE<0,1,2,3,4,5>(A0,A,*owner,1.0);

	// Initialize static fields

	CopyFieldData(scratch1,Element(0),sources,Element(0));
	ellipticSolver->Solve(scratch2,scratch1,-1.0);
	add_grad<0,0,1,2>(scratch2,A,*owner,-1.0);

	// Clean up ghost cells

	SetSingularPointsE();
	SetSingularPointsB();
	A.CopyFromNeighbors();
	A.ApplyBoundaryCondition();

	// Time centered fields for particle pusher

	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,true))
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
				A(1,j,k,0) = 0.0;
	}
	if (owner->gridGeometry==tw::grid::spherical && owner->X(0,1)<0.0)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (j=lfg[2];j<=ufg[2];j++)
				A(1,j,k,0) = 0.0;
	}
	if (owner->gridGeometry==tw::grid::spherical && owner->X(0,2)<0.0)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (i=lfg[1];i<=ufg[1];i++)
				A(i,1,k,1) = 0.0;
	}
	if (owner->gridGeometry==tw::grid::spherical && owner->X(ufg[2],2)>pi)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (i=lfg[1];i<=ufg[1];i++)
				A(i,ufg[2],k,1) = 0.0;
	}
}

void CurvilinearDirectSolver::SetSingularPointsB()
{
	// Assume that F holds the old staggered B-field at this point
	// Then use curlE = -dB/dt to update B
	tw::Int i,j,k;
	if (owner->gridGeometry==tw::grid::cylindrical && owner->X(0,1)<0.0)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (j=lfg[2];j<=ufg[2];j++)
			{
				A(1,j,k,4) = 0.0;
				A(1,j,k,5) = F(1,j,k,5) - 2.0*dt*A(1,j,k,1)/owner->X(1,1);
			}
	}
	if (owner->gridGeometry==tw::grid::spherical && owner->X(0,2)<0.0)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (i=lfg[1];i<=ufg[1];i++)
			{
				A(i,1,k,5) = 0.0;
				A(i,1,k,3) = F(i,1,k,3) - 2.0*dt*A(i,1,k,2)/(owner->X(i,1)*sin(owner->X(1,2)));
			}
	}
	if (owner->gridGeometry==tw::grid::spherical && owner->X(ufg[2],2)>pi)
	{
		for (k=lfg[3];k<=ufg[3];k++)
			for (i=lfg[1];i<=ufg[1];i++)
			{
				A(i,ufg[2],k,5) = 0.0;
				A(i,ufg[2],k,3) = F(i,ufg[2],k,3) + 2.0*dt*A(i,dim[2],k,2)/(owner->X(i,1)*sin(owner->X(dim[2],2)));
			}
	}
}

void CurvilinearDirectSolver::Update()
{
	tw::vec3 S;

	// Add electric antenna currents

	if (conductor.size())
	{
		#ifdef USE_OPENCL
		sources.ReceiveFromComputeBuffer();
		#endif
		for (auto c : conductor)
			if (c->currentType==EM::current::electric)
				c->DepositSources(sources, owner->elapsedTime, dt);
		#ifdef USE_OPENCL
		sources.SendToComputeBuffer();
		#endif
	}

	Electromagnetic::Update();

	// Advance the electric field

	add_curlB<3,4,5,0,1,2>(A,A,*owner,dt);
	AddMulFieldData(A,Element(0,2),sources,Element(1,3),-dt);
	SetSingularPointsE();
	A.CopyFromNeighbors(Element(0,2));

	// Save the old magnetic field so final fields can be centered

	CopyFieldData(F,Element(3,5),A,Element(3,5));

	// Advance the magnetic field

	add_curlE<0,1,2,3,4,5>(A,A,*owner,-dt);
	SetSingularPointsB();
	A.UpwardCopy(tw::grid::x,Element(3,5),1);
	A.UpwardCopy(tw::grid::y,Element(3,5),1);
	A.UpwardCopy(tw::grid::z,Element(3,5),1);

	// Setup the final field for the particle push

	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,true))
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
	directives.Add("dims",new tw::input::Numbers<tw::Int>(&dims[0],3));
	directives.Add("period",new tw::input::Int(&period));
}

void FarFieldDiagnostic::Initialize()
{
	if (J4==NULL)
		throw tw::FatalError("Far field diagnostic did not find a source field.");
	tw::vec3 corner(bounds[0],bounds[2],bounds[4]);
	tw::vec3 size(bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]);
	A.Initialize(DiscreteSpace(dims[0],dims[1],dims[2],corner,size,1),owner);
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
	const tw::Float tp = owner->elapsedTime;
	const tw::Float dtau = dt * tw::Float(period);

	if (owner->stepNow % period==0)
		for (auto farCell : InteriorCellRange(A))
		{
			const tw::Float tNow = bounds[0] + dx(A)*tw::Float(farCell.dcd1());
			const tw::Float thetaNow = bounds[2] + dy(A)*tw::Float(farCell.dcd2());
			const tw::Float phiNow = bounds[4] + dz(A)*tw::Float(farCell.dcd3());
			const tw::vec3 n(sin(thetaNow)*cos(phiNow),sin(thetaNow)*sin(phiNow),cos(thetaNow));
			#pragma omp parallel
			{
				for (auto s : StripRange(*this,3,strongbool::no))
				{
					std::valarray<tw::Float> j4(4);
					tw::vec3 rp = owner->Pos(s,1);
					rp.z = ((tp - tNow) - rp.x*n.x - rp.y*n.y)/n.z;
					const tw::Float dS = owner->dS(s,1,3)/n.z;
					if (rp.z > zmin && rp.z < zmax)
					{
						weights_3D w;
						owner->GetWeights(&w,rp);
						J4->Interpolate(j4,w);
						A(farCell) += tw::vec3(j4[1],j4[2],j4[3]) * dS * dtau / radius;
					}
				}
			}
		}

	if (owner->IsLastStep())
		SpecialReport();
}

void FarFieldDiagnostic::SpecialReport()
{
	std::string fileName;
	tw::Int master = 0;
	tw::Int curr = owner->strip[0].Get_rank();

	// Each node has different contributions to the same far field grid.
	// We merely need to form the superposition of all these grids.
	Vec3Field accum;
	accum.Initialize(A,owner);
	owner->strip[0].Sum(&A(0,0,0),&accum(0,0,0),sizeof(tw::vec3)*accum.TotalCells(),master);

	if (curr==master)
	{
		// put vector potential in coulomb gauge and spherical coordinates
		for (auto farCell : InteriorCellRange(accum))
		{
			const tw::Float theta = bounds[2] + dy(A)*tw::Float(farCell.dcd2());
			const tw::Float phi = bounds[4] + dz(A)*tw::Float(farCell.dcd3());
			const tw::vec3 nr( sin(theta)*cos(phi) , sin(theta)*sin(phi) , cos(theta) );
			const tw::vec3 nq( cos(theta)*cos(phi) , cos(theta)*sin(phi) , -sin(theta) );
			const tw::vec3 nf( -sin(phi) , cos(phi) , 0.0 );
			const tw::vec3 ACG = nr | (accum(farCell) | nr); // form coulomb gauge A
			accum(farCell) = tw::vec3( ACG^nr , ACG^nq , ACG^nf ); // put in spherical coordinates
		}

		npy_writer writer;
		tw::Int shape[4] = { 1 , dims[0] , dims[1] , dims[2] };
		std::valarray<float> gData(shape[1]*shape[2]*shape[3]);

		fileName = name + "-Atheta.npy";
		writer.write_header(fileName,shape);
		for (tw::Int i=1;i<=shape[1];i++)
			for (tw::Int j=1;j<=shape[2];j++)
				for (tw::Int k=1;k<=shape[3];k++)
					gData[(i-1)*shape[2]*shape[3] + (j-1)*shape[3] + (k-1)] = accum(i,j,k).y;
		writer.add_frame(fileName,(char*)&gData[0],shape);

		fileName = name + "-Aphi.npy";
		writer.write_header(fileName,shape);
		for (tw::Int i=1;i<=shape[1];i++)
			for (tw::Int j=1;j<=shape[2];j++)
				for (tw::Int k=1;k<=shape[3];k++)
					gData[(i-1)*shape[2]*shape[3] + (j-1)*shape[3] + (k-1)] = accum(i,j,k).z;
		writer.add_frame(fileName,(char*)&gData[0],shape);
	}
}
