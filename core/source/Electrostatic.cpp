#include "simulation.h"
#include "fieldSolve.h"
#include "electrostatic.h"
using namespace tw::bc;


///////////////////////////////
//                           //
// ELECTROSTATIC FIELD GROUP //
// (General Curvilinear)     //
//                           //
///////////////////////////////


Electrostatic::Electrostatic(const std::string& name,Simulation* sim):FieldSolver(name,sim)
{
	typeCode = tw::module_type::electrostatic;
	phi.Initialize(*this,owner);
	source.Initialize(*this,owner);
	Ef.Initialize(*this,owner);
	J4.Initialize(4,*this,owner);

	Ef.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	Ef.SetBoundaryConditions(Element(0),tw::grid::x,fld::none,fld::normalFluxFixed);

	Ef.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	Ef.SetBoundaryConditions(Element(1),tw::grid::y,fld::none,fld::normalFluxFixed);

	Ef.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	Ef.SetBoundaryConditions(Element(2),tw::grid::z,fld::none,fld::normalFluxFixed);

	J4.SetBoundaryConditions(tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Element(1),tw::grid::x,fld::normalFluxFixed,fld::normalFluxFixed);

	J4.SetBoundaryConditions(tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Element(2),tw::grid::y,fld::normalFluxFixed,fld::normalFluxFixed);

	J4.SetBoundaryConditions(tw::grid::z,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Element(3),tw::grid::z,fld::normalFluxFixed,fld::normalFluxFixed);
}

void Electrostatic::Initialize()
{
	FieldSolver::Initialize();
	ellipticSolver->SetBoundaryConditions(phi);
	for (auto c : conductor)
		ellipticSolver->FixPotential(phi,c->theRgn,c->Voltage(owner->elapsedTime));
	SetupInitialPotential();

	Update();
}

void Electrostatic::ExchangeResources()
{
	PublishResource(&Ef,"electrostatic:E");
	PublishResource(&phi,"electrostatic:phi");
	PublishResource(&J4,"electromagnetic:sources"); // use electromagnetic format for sources
	PublishResource(this,"electrostatic:module");
}

void Electrostatic::Reset()
{
	J4 = 0.0;
}

void Electrostatic::ReadCheckpoint(std::ifstream& inFile)
{
	FieldSolver::ReadCheckpoint(inFile);

	phi.ReadCheckpoint(inFile);
	J4.ReadCheckpoint(inFile);
}

void Electrostatic::WriteCheckpoint(std::ofstream& outFile)
{
	FieldSolver::WriteCheckpoint(outFile);

	phi.WriteCheckpoint(outFile);
	J4.WriteCheckpoint(outFile);
}

void Electrostatic::Update()
{
	J4.DepositFromNeighbors();
	J4.ApplyFoldingCondition();
	conserved_current_to_dens<0,1,2,3>(J4,*owner);
	J4.ApplyBoundaryCondition();

	for (auto c : conductor)
		ellipticSolver->FixPotential(phi,c->theRgn,c->Voltage(owner->elapsedTime));

	CopyFieldData(source,Element(0),J4,Element(0));
	ellipticSolver->Solve(phi,source,-1.0);

	ComputeFinalFields();
}

void Electrostatic::ComputeFinalFields()
{
	assign_grad<0,0,1,2>(phi,Ef,*owner,-1.0);

	Ef.CopyFromNeighbors();
	Ef.ApplyBoundaryCondition();
}

void Electrostatic::Report(Diagnostic& diagnostic)
{
	FieldSolver::Report(diagnostic);
	ScalarField energyDensity;
	energyDensity.Initialize(*this,owner);
	for (auto cell : EntireCellRange(*this))
		energyDensity(cell) = 0.5*Norm(Ef(cell));
	diagnostic.VolumeIntegral("FieldEnergy",energyDensity,0);
	diagnostic.VolumeIntegral("TotalCharge",J4,0);
	diagnostic.Field("phi",phi,0);
	diagnostic.Field("Ex",Ef,0);
	diagnostic.Field("Ey",Ef,1);
	diagnostic.Field("Ez",Ef,2);
	diagnostic.Field("rho",J4,0);
	diagnostic.Field("Jx",J4,1);
	diagnostic.Field("Jy",J4,2);
	diagnostic.Field("Jz",J4,3);
}

void Electrostatic::SetupInitialPotential()
{
	// As an optimization, something can go here to try to make an initial guess other than zero.
}
