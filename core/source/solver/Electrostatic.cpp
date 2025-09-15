module;

#include "tw_includes.h"

export module electrostatic;
import twmodule;
import fields;
import elliptic;
import field_solve;
import diagnostics;

using namespace tw::bc;

export struct Electrostatic:FieldSolver
{
	Field F,J4;
	ScalarField phi,source;

	Electrostatic(const std::string& name,Simulation* sim);
	virtual void Initialize();
	virtual void ExchangeResources();
	virtual void Reset();
	virtual void Update();
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void SetupInitialPotential();
	virtual void ComputeFinalFields();

	virtual void Report(Diagnostic&);
};

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

	F.SetBoundaryConditions(All(F),tw::grid::x,fld::neumannWall,fld::neumannWall);
	F.SetBoundaryConditions(Rng(0),tw::grid::x,fld::none,fld::normalFluxFixed);

	F.SetBoundaryConditions(All(F),tw::grid::y,fld::neumannWall,fld::neumannWall);
	F.SetBoundaryConditions(Rng(1),tw::grid::y,fld::none,fld::normalFluxFixed);

	F.SetBoundaryConditions(All(F),tw::grid::z,fld::neumannWall,fld::neumannWall);
	F.SetBoundaryConditions(Rng(2),tw::grid::z,fld::none,fld::normalFluxFixed);

	J4.SetBoundaryConditions(All(J4),tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Rng(1),tw::grid::x,fld::normalFluxFixed,fld::normalFluxFixed);

	J4.SetBoundaryConditions(All(J4),tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Rng(2),tw::grid::y,fld::normalFluxFixed,fld::normalFluxFixed);

	J4.SetBoundaryConditions(All(J4),tw::grid::z,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Rng(3),tw::grid::z,fld::normalFluxFixed,fld::normalFluxFixed);
}

void Electrostatic::Initialize()
{
	F = 0.0; // F(3:6) must remain 0
	FieldSolver::Initialize();
	ellipticSolver->SetFieldsBoundaryConditions(phi,Rng(0));
	for (auto c : conductors)
		ellipticSolver->FixPotential(phi,c->theRgn,c->Voltage(owner->WindowPos(0)));
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
	J4.DepositFromNeighbors(All(J4));
	J4.ApplyFoldingCondition(All(J4));
	conserved_current_to_dens<0,1,2,3>(1,J4,*owner);
	J4.ApplyBoundaryCondition(All(J4));

	for (auto c : conductors)
		ellipticSolver->FixPotential(phi,c->theRgn,c->Voltage(owner->WindowPos(0)));

	CopyFieldData(source,Rng(0),J4,Rng(0));
	ellipticSolver->Solve(phi,source,-1.0);

	ComputeFinalFields();
}

void Electrostatic::ComputeFinalFields()
{
	assign_grad<0,0,1,2>(1,phi,F,*owner,-1.0);

	F.CopyFromNeighbors(All(F));
	F.ApplyBoundaryCondition(All(F));
}

void Electrostatic::Report(Diagnostic& diagnostic)
{
	FieldSolver::Report(diagnostic);
	ScalarField energyDensity;
	energyDensity.Initialize(*this,owner);
	for (auto cell : EntireCellRange(*this,1))
		energyDensity(cell) = 0.5*Norm(F.Vec3(cell,0));
	diagnostic.VolumeIntegral("FieldEnergy",energyDensity,1,0);
	diagnostic.VolumeIntegral("TotalCharge",J4,1,0);
	diagnostic.ReportField("phi",phi,1,0,tw::dims::scalar_potential,"$\\phi$");
	diagnostic.ReportField("Ex",F,1,0,tw::dims::electric_field,"$E_x$");
	diagnostic.ReportField("Ey",F,1,1,tw::dims::electric_field,"$E_y$");
	diagnostic.ReportField("Ez",F,1,2,tw::dims::electric_field,"$E_z$");
	diagnostic.ReportField("rho",J4,1,0,tw::dims::charge_density,"$\\rho$");
	diagnostic.ReportField("Jx",J4,1,1,tw::dims::current_density,"$j_x$");
	diagnostic.ReportField("Jy",J4,1,2,tw::dims::current_density,"$j_y$");
	diagnostic.ReportField("Jz",J4,1,3,tw::dims::current_density,"$j_z$");
}

void Electrostatic::SetupInitialPotential()
{
	// As an optimization, something can go here to try to make an initial guess other than zero.
}
