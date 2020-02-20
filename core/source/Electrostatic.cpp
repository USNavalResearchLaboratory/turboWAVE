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

	Ef.SetBoundaryConditions(tw::dom::xAxis,fld::neumannWall,fld::neumannWall);
	Ef.SetBoundaryConditions(Element(0),tw::dom::xAxis,fld::none,fld::normalFluxFixed);

	Ef.SetBoundaryConditions(tw::dom::yAxis,fld::neumannWall,fld::neumannWall);
	Ef.SetBoundaryConditions(Element(1),tw::dom::yAxis,fld::none,fld::normalFluxFixed);

	Ef.SetBoundaryConditions(tw::dom::zAxis,fld::neumannWall,fld::neumannWall);
	Ef.SetBoundaryConditions(Element(2),tw::dom::zAxis,fld::none,fld::normalFluxFixed);

	J4.SetBoundaryConditions(tw::dom::xAxis,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Element(1),tw::dom::xAxis,fld::normalFluxFixed,fld::normalFluxFixed);

	J4.SetBoundaryConditions(tw::dom::yAxis,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Element(2),tw::dom::yAxis,fld::normalFluxFixed,fld::normalFluxFixed);

	J4.SetBoundaryConditions(tw::dom::zAxis,fld::dirichletCell,fld::dirichletCell);
	J4.SetBoundaryConditions(Element(3),tw::dom::zAxis,fld::normalFluxFixed,fld::normalFluxFixed);
}

void Electrostatic::Initialize()
{
	FieldSolver::Initialize();
	ellipticSolver->SetBoundaryConditions(phi);
	for (auto conductor : owner->conductor)
		ellipticSolver->FixPotential(phi,conductor->theRgn,conductor->Voltage(owner->elapsedTime));
	SetupInitialPotential();

	if (!owner->restarted)
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

void Electrostatic::ReadData(std::ifstream& inFile)
{
	FieldSolver::ReadData(inFile);

	phi.ReadData(inFile);
	J4.ReadData(inFile);
}

void Electrostatic::WriteData(std::ofstream& outFile)
{
	FieldSolver::WriteData(outFile);

	phi.WriteData(outFile);
	J4.WriteData(outFile);
}

void Electrostatic::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "FieldEnergy TotalCharge Current ";
}

void Electrostatic::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	tw::Int i,j,k;
	tw::Float fieldEnergy,totalCharge,current;
	tw::vec3 pos;

	tw::Int x0,x1,y0,y1,z0,z1;
	theRgn.GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	fieldEnergy = 0.0;
	totalCharge = 0.0;
	current = 0.0;
	for (k=z0;k<=z1;k++)
		for (j=y0;j<=y1;j++)
			for (i=x0;i<=x1;i++)
			{
				pos = owner->Pos(i,j,k);
				if (theRgn.Inside(pos,*owner))
				{
					fieldEnergy += 0.5 * owner->dS(i,j,k,0) * Norm(Ef(i,j,k));
					totalCharge += owner->dS(i,j,k,0) * J4(i,j,k,0);
					current += owner->dS(i,j,k,3) * 0.5 * (J4(i,j,k-1,3) + J4(i,j,k,3)) / tw::Float(z1 - z0 + 1);
				}
			}

	cols.push_back(fieldEnergy); avg.push_back(false);
	cols.push_back(totalCharge); avg.push_back(false);
	cols.push_back(current); avg.push_back(false);
}

void Electrostatic::Update()
{
	J4.DepositFromNeighbors();
	J4.ApplyFoldingCondition();
	conserved_current_to_dens<0,1,2,3>(J4,*owner);
	J4.ApplyBoundaryCondition();

	for (auto conductor : owner->conductor)
		ellipticSolver->FixPotential(phi,conductor->theRgn,conductor->Voltage(owner->elapsedTime));

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

void Electrostatic::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader("phi",box);
	owner->WriteBoxDataHeader("Ex",box);
	owner->WriteBoxDataHeader("Ey",box);
	owner->WriteBoxDataHeader("Ez",box);
	owner->WriteBoxDataHeader("rho",box);
	owner->WriteBoxDataHeader("Jx",box);
	owner->WriteBoxDataHeader("Jy",box);
	owner->WriteBoxDataHeader("Jz",box);
}

void Electrostatic::BoxDiagnose(GridDataDescriptor* box)
{
	owner->WriteBoxData("phi",box,&phi(0,0,0),phi.Stride());
	owner->WriteBoxData("Ex",box,&Ef(0,0,0).x,Ef.Stride());
	owner->WriteBoxData("Ey",box,&Ef(0,0,0).y,Ef.Stride());
	owner->WriteBoxData("Ez",box,&Ef(0,0,0).z,Ef.Stride());
	owner->WriteBoxData("rho",box,&J4(0,0,0,0),J4.Stride());
	owner->WriteBoxData("Jx",box,&J4(0,0,0,1),J4.Stride());
	owner->WriteBoxData("Jy",box,&J4(0,0,0,2),J4.Stride());
	owner->WriteBoxData("Jz",box,&J4(0,0,0,3),J4.Stride());
}

void Electrostatic::PointDiagnosticHeader(std::ofstream& outFile)
{
	outFile << "phi rho ";
}

void Electrostatic::PointDiagnose(std::ofstream& outFile,const weights_3D& w)
{
	tw::Float valNow;
	std::valarray<tw::Float> j4(4);

	phi.Interpolate(&valNow,w);
	outFile << valNow << " ";
	J4.Interpolate(j4,w);
	outFile << j4[0] << " ";
}

void Electrostatic::SetupInitialPotential()
{
	// As an optimization, something can go here to try to make an initial guess other than zero.
}
