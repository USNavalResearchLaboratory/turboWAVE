#include "../simulation.h"
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
	if (native.native!=tw::units::plasma)
		throw tw::FatalError("Electrostatic module requires <native units = plasma>");

	phi.Initialize(*this,owner);
	source.Initialize(*this,owner);
	F.Initialize(6,*this,owner); // only use elements 0,1,2; set 3,4,5 to zero.
	J4.Initialize(4,*this,owner);

	F.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	F.SetBoundaryConditions(Element(0),tw::grid::x,fld::none,fld::normalFluxFixed);

	F.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	F.SetBoundaryConditions(Element(1),tw::grid::y,fld::none,fld::normalFluxFixed);

	F.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	F.SetBoundaryConditions(Element(2),tw::grid::z,fld::none,fld::normalFluxFixed);

	J4.SetBoundaryConditions(tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Element(1),tw::grid::x,fld::normalFluxFixed,fld::normalFluxFixed);

	J4.SetBoundaryConditions(tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Element(2),tw::grid::y,fld::normalFluxFixed,fld::normalFluxFixed);

	J4.SetBoundaryConditions(tw::grid::z,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Element(3),tw::grid::z,fld::normalFluxFixed,fld::normalFluxFixed);
}

void Electrostatic::Initialize()
{
	F = 0.0; // F(3:6) must remain 0
	FieldSolver::Initialize();
	ellipticSolver->SetFieldsBoundaryConditions(phi,Element(0));
	for (auto c : conductor)
		ellipticSolver->FixPotential(phi,c->theRgn,c->Voltage(owner->elapsedTime));
	SetupInitialPotential();

	Update();
}

void Electrostatic::ExchangeResources()
{
	PublishResource(&F,"electromagnetic:F"); // use electromagnetic format with B=0
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
	assign_grad<0,0,1,2>(phi,F,*owner,-1.0);

	F.CopyFromNeighbors();
	F.ApplyBoundaryCondition();
}

void Electrostatic::Report(Diagnostic& diagnostic)
{
	FieldSolver::Report(diagnostic);
	ScalarField energyDensity;
	energyDensity.Initialize(*this,owner);
	for (auto cell : EntireCellRange(*this))
		energyDensity(cell) = 0.5*Norm(F.Vec3(cell,0));
	diagnostic.VolumeIntegral("FieldEnergy",energyDensity,0);
	diagnostic.VolumeIntegral("TotalCharge",J4,0);
	diagnostic.Field("phi",phi,0,tw::dims::scalar_potential,"$\\phi$");
	diagnostic.Field("Ex",F,0,tw::dims::electric_field,"$E_x$");
	diagnostic.Field("Ey",F,1,tw::dims::electric_field,"$E_y$");
	diagnostic.Field("Ez",F,2,tw::dims::electric_field,"$E_z$");
	diagnostic.Field("rho",J4,0,tw::dims::charge_density,"$\\rho$");
	diagnostic.Field("Jx",J4,1,tw::dims::current_density,"$j_x$");
	diagnostic.Field("Jy",J4,2,tw::dims::current_density,"$j_y$");
	diagnostic.Field("Jz",J4,3,tw::dims::current_density,"$j_z$");
}

void Electrostatic::SetupInitialPotential()
{
	// As an optimization, something can go here to try to make an initial guess other than zero.
}
