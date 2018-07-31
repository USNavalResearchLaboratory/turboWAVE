#include "sim.h"
#include "fieldSolve.h"
#include "electrostatic.h"


///////////////////////////////
//                           //
// ELECTROSTATIC FIELD GROUP //
// (General Curvilinear)     //
//                           //
///////////////////////////////


Electrostatic::Electrostatic(Grid* theGrid):FieldSolver(theGrid)
{
	name = "Electrostatic";
	typeCode = electrostatic;
	phi.Initialize(*this,owner);
	source.Initialize(*this,owner);
	Ef.Initialize(*this,owner);
	sources.Initialize(4,*this,owner);
	lbc.resize(Num(1));
	rbc.resize(Num(1));

	Ef.SetBoundaryConditions(xAxis,neumannWall,neumannWall);
	Ef.SetBoundaryConditions(yAxis,neumannWall,neumannWall);
	Ef.SetBoundaryConditions(zAxis,neumannWall,neumannWall);

	sources.SetBoundaryConditions(xAxis,dirichletCell,dirichletCell);
	sources.SetBoundaryConditions(Element(1),xAxis,normalFluxFixed,normalFluxFixed);

	sources.SetBoundaryConditions(yAxis,dirichletCell,dirichletCell);
	sources.SetBoundaryConditions(Element(2),yAxis,normalFluxFixed,normalFluxFixed);

	sources.SetBoundaryConditions(zAxis,dirichletCell,dirichletCell);
	sources.SetBoundaryConditions(Element(3),zAxis,normalFluxFixed,normalFluxFixed);
	
	if (owner->gridGeometry!=cartesian)
	{
		Ef.SetBoundaryConditions(Element(0),xAxis,dirichletWall,neumannWall);
		sources.SetBoundaryConditions(Element(1),xAxis,dirichletWall,neumannWall);
	}

	electrodeRadius = 0.0; // indicates "no electrode"
	electrodePotential = 0.0;
	slewRate = 0.0;

	x0 = x1 = y0 = y1 = neumannWall;
	z0 = natural;
	z1 = natural;
	lbc = 0.0;
	rbc = 0.0;
}

void Electrostatic::Initialize()
{
	FieldSolver::Initialize();

	ellipticSolver->SetBoundaryConditions(phi,x0,x1,y0,y1,z0,z1);

	SetupElectrodePotential();
	SetupInitialPotential();

	if (!owner->restarted)
		Update();
}

void Electrostatic::ExchangeResources()
{
	PublishResource(&Ef,"electrostatic:E");
	PublishResource(&phi,"electrostatic:phi");
	PublishResource(&sources,"electromagnetic:sources"); // use electromagnetic format for sources
	PublishResource(this,"electrostatic:module");
}

void Electrostatic::Reset()
{
	sources = 0.0;
}

void Electrostatic::ReadInputFileTerm(std::stringstream& inputString,std::string& command)
{
	std::string word;
	tw::Float lbc0,rbc0;
	FieldSolver::ReadInputFileTerm(inputString,command);
	if (command=="external") // eg, external potential = ( 0.0 , 1.0 )
	{
		inputString >> word >> word >> lbc0 >> rbc0;
		lbc = lbc0;
		rbc = rbc0;
	}
	if (command=="electrode") // eg, electrode : radius = 1e5 , potential = 0.5 , slew rate = 0.01
	{
		inputString >> word >> word >> word >> electrodeRadius;
		inputString >> word >> word >> electrodePotential;
		inputString >> word >> word >> word >> slewRate;
	}
}

void Electrostatic::ReadData(std::ifstream& inFile)
{
	FieldSolver::ReadData(inFile);
	inFile.read((char *)&lbc[0],Num(1)*sizeof(tw::Float));
	inFile.read((char *)&rbc[0],Num(1)*sizeof(tw::Float));
	inFile.read((char *)&electrodeRadius,sizeof(tw::Float));
	inFile.read((char *)&electrodePotential,sizeof(tw::Float));
	inFile.read((char *)&slewRate,sizeof(tw::Float));

	phi.ReadData(inFile);
	sources.ReadData(inFile);
}

void Electrostatic::WriteData(std::ofstream& outFile)
{
	FieldSolver::WriteData(outFile);
	outFile.write((char *)&lbc[0],Num(1)*sizeof(tw::Float));
	outFile.write((char *)&rbc[0],Num(1)*sizeof(tw::Float));
	outFile.write((char *)&electrodeRadius,sizeof(tw::Float));
	outFile.write((char *)&electrodePotential,sizeof(tw::Float));
	outFile.write((char *)&slewRate,sizeof(tw::Float));

	phi.WriteData(outFile);
	sources.WriteData(outFile);
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
					totalCharge += owner->dS(i,j,k,0) * sources(i,j,k,0);
					current += owner->dS(i,j,k,3) * 0.5 * (sources(i,j,k-1,3) + sources(i,j,k,3)) / tw::Float(z1 - z0 + 1);
				}
			}

	cols.push_back(fieldEnergy); avg.push_back(false);
	cols.push_back(totalCharge); avg.push_back(false);
	cols.push_back(current); avg.push_back(false);
}

void Electrostatic::Update()
{
	tw::Int i;
	tw::Float t = owner->elapsedTime;

	sources.DepositFromNeighbors();
	sources.ApplyFoldingCondition();
	conserved_current_to_dens<0,1,2,3>(sources,*owner);
	sources.ApplyBoundaryCondition();

	if (suppressNextUpdate)
	{
		suppressNextUpdate = false;
		return;
	}

	SetupElectrodePotential();
	ellipticSolver->lbc = lbc;
	ellipticSolver->rbc = rbc;
	ellipticSolver->TransformBoundaryValues();

	for (i=0;i<owner->conductor.size();i++)
		ellipticSolver->FixPotential(phi,owner->conductor[i]->theRgn,owner->conductor[i]->Voltage(t));

	CopyFieldData(source,Element(0),sources,Element(0));
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
	owner->WriteBoxData("rho",box,&sources(0,0,0,0),sources.Stride());
	owner->WriteBoxData("Jx",box,&sources(0,0,0,1),sources.Stride());
	owner->WriteBoxData("Jy",box,&sources(0,0,0,2),sources.Stride());
	owner->WriteBoxData("Jz",box,&sources(0,0,0,3),sources.Stride());
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
	sources.Interpolate(j4,w);
	outFile << j4[0] << " ";
}

void Electrostatic::SetupElectrodePotential()
{
	tw::Int i;
	tw::Float rNow,r0,phiNow;
	r0 = electrodeRadius;
	phiNow = electrodePotential + slewRate*owner->elapsedTime;
	if (r0!=0.0)
		for (i=N0(1);i<=N1(1);i++)
		{
			rNow = owner->X(i,1);
			if (rNow/r0<0.5*pi)
				lbc[i] = phiNow*sqr(cos(rNow/r0));
			else
				lbc[i] = 0.0;
		}
}

void Electrostatic::SetupInitialPotential()
{
	tw::Int i,j,k;

	for (CellIterator cell(*this,true);cell<cell.end();++cell)
	{
		cell.Decode(&i,&j,&k);
		phi(i,j,k) = lbc[i] + (rbc[i]-lbc[i])*(owner->X(k,3)-GlobalCorner(*owner).z)/GlobalPhysicalSize(*owner).z;
	}
}
