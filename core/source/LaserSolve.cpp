#include "sim.h"
#include "fieldSolve.h"
#include "laserSolve.h"


//////////////////////////////
//                          //
//     LASER SOLVER BASE    //
//                          //
//////////////////////////////


LaserSolver::LaserSolver(Grid* theGrid):Module(theGrid)
{
	updateSequencePriority = 3;
	laserFreq = 10.0;
	polarizationType = linearPolarization;
	propagator = (LaserPropagator*)owner->AddPrivateTool(adiPropagator);

	a0.Initialize(*this,owner);
	a1.Initialize(*this,owner);
	chi.Initialize(*this,owner);

	// set boundary conditions in propagators and child classes
}

LaserSolver::~LaserSolver()
{
	owner->RemoveTool(propagator);
}

void LaserSolver::ExchangeResources()
{
	PublishResource(&a0,"laser:a0");
	PublishResource(&a1,"laser:a1");
	PublishResource(&chi,"laser:chi");
	PublishResource(&laserFreq,"laser:carrierFrequency");
	PublishResource(&polarizationType,"laser:polarizationType");
}

void LaserSolver::Initialize()
{
	tw::Int i,j,k,s;
	tw::vec3 pos;
	tw::Float polarizationFactor;

	Module::Initialize();
	propagator->SetData(laserFreq,dt,polarizationType,owner->movingWindow);
	propagator->SetBoundaryConditions(a0,a1,chi);

	if (owner->restarted)
		return;

	if (polarizationType==circularPolarization)
		polarizationFactor = 1.414;
	else
		polarizationFactor = 1.0;

	for (CellIterator cell(*this,true);cell<cell.end();++cell)
		for (s=0;s<owner->pulse.size();s++)
		{
			pos = owner->Pos(cell);
			pos.z = owner->ToLab(pos.z,-dth);
			a0(cell) += polarizationFactor*owner->pulse[s]->VectorPotentialEnvelope(-dth,pos);
			pos = owner->Pos(cell);
			pos.z = owner->ToLab(pos.z,dth);
			a1(cell) += polarizationFactor*owner->pulse[s]->VectorPotentialEnvelope(dth,pos);
		}
}

void LaserSolver::Update()
{
	propagator->Advance(a0,a1,chi);
}

void LaserSolver::Reset()
{
	chi = tw::Complex(0,0);
}

void LaserSolver::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;
	Module::ReadInputFileDirective(inputString,command);
	if (command=="propagator") // eg, propagator = eigenmode
	{
		inputString >> word >> word;
		if (word=="spectral" || word=="FFT" || word=="eigenmode")
		{
			*owner->tw_out << "Use eigenmode propagator" << std::endl;
			owner->RemoveTool(propagator);
			propagator = (LaserPropagator*)owner->AddPrivateTool(eigenmodePropagator);
		}
		if (word=="ADI" || word=="adi")
		{
			*owner->tw_out << "Use ADI propagator" << std::endl;
			owner->RemoveTool(propagator);
			propagator = (LaserPropagator*)owner->AddPrivateTool(adiPropagator);
		}
	}
	if (command=="carrier") // eg, carrier frequency = 30.0
	{
		inputString >> word >> word;
		inputString >> laserFreq;
		if (propagator!=NULL)
			propagator->w0 = laserFreq;
	}
	if (command=="polarization") // eg, polarization = linear
	{
		inputString >> word >> word;
		if (word=="linear")
			polarizationType = linearPolarization;
		if (word=="circular")
			polarizationType = circularPolarization;
		if (word=="radial")
			polarizationType = radialPolarization;
		if (word!="linear" && word!="circular" && word!="radial")
		{
			*owner->tw_out << "abort: invalid polarization" << std::endl;
			abort();
		}
	}
}

void LaserSolver::ReadData(std::ifstream& inFile)
{
	Module::ReadData(inFile);
	inFile.read((char *)&laserFreq,sizeof(tw::Float));
	inFile.read((char *)&polarizationType,sizeof(tw_polarization_type));
	a0.ReadData(inFile);
	a1.ReadData(inFile);

 	tw_tool propagatorType;
 	inFile.read((char*)&propagatorType,sizeof(tw_tool));
	owner->RemoveTool(propagator);
 	propagator = (LaserPropagator*)owner->AddPrivateTool(propagatorType);
}

void LaserSolver::WriteData(std::ofstream& outFile)
{
	Module::WriteData(outFile);
	outFile.write((char *)&laserFreq,sizeof(tw::Float));
	outFile.write((char *)&polarizationType,sizeof(tw_polarization_type));
	a0.WriteData(outFile);
	a1.WriteData(outFile);
	outFile.write((char*)&propagator->typeCode,sizeof(tw_tool));
}

void LaserSolver::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader("a2",box);
	owner->WriteBoxDataHeader("e_real",box);
	owner->WriteBoxDataHeader("e_imag",box);
	owner->WriteBoxDataHeader("chi_real",box);
	owner->WriteBoxDataHeader("chi_imag",box);
}

void LaserSolver::PointDiagnosticHeader(std::ofstream& outFile)
{
	outFile << "ar ai ";
}

void LaserSolver::PointDiagnose(std::ofstream& outFile,const weights_3D& w)
{
	tw::Complex aNow;
	a1.Interpolate(&aNow,w);
	outFile << real(aNow) << " " << imag(aNow) << " ";
}





////////////////////////
//                    //
//  QS Laser Solver  //
//                    //
////////////////////////


QSSolver::QSSolver(Grid* theGrid) : LaserSolver(theGrid)
{
	name = "Quasistatic-Laser";
	typeCode = qsLaser;
}

void QSSolver::EnergyHeadings(std::ofstream& outFile)
{
}

void QSSolver::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
}

void QSSolver::BoxDiagnose(GridDataDescriptor* box)
{
	ScalarField er,ei;
	er.Initialize(*this,owner);
	ei.Initialize(*this,owner);
	tw::Complex dadt,anow;
	for (CellIterator cell(*this,false);cell<cell.end();++cell)
	{
		dadt = -tw::Complex(a1(cell,0,3),a1(cell,1,3));
		anow = a1(cell);
		er(cell) = -real(dadt - ii*laserFreq*anow);
		ei(cell) = -imag(dadt - ii*laserFreq*anow);
	}

	owner->WriteBoxData("e_real",box,&er(0,0,0),er.Stride());
	owner->WriteBoxData("e_imag",box,&ei(0,0,0),ei.Stride());

	for (CellIterator cell(*this,false);cell<cell.end();++cell)
		er(cell) = norm(a1(cell));

	owner->WriteBoxData("a2",box,&er(0,0,0),er.Stride());
}


//////////////////////////////
//                          //
//     PGC LASER SOLVER     //
//                          //
//////////////////////////////


PGCSolver::PGCSolver(Grid* theGrid):LaserSolver(theGrid)
{
	name = "PGC-Laser";
	typeCode = pgcLaser;

	F.Initialize(8,*this,owner);
}

void PGCSolver::ExchangeResources()
{
	LaserSolver::ExchangeResources();
	PublishResource(&F,"laser:F");
}

void PGCSolver::Initialize()
{
	LaserSolver::Initialize();

	// Here we deal with boundary conditions particular to PGC
	// B.C.'s for laser amplitude and sources are dealt with in propagator objects

	// longitudinal boundary conditions
	if (owner->movingWindow)
	{
		F.SetBoundaryConditions(Element(0,5),zAxis,neumannWall,dirichletWall);
	}
	else
	{
		F.SetBoundaryConditions(Element(0,5),zAxis,neumannWall,neumannWall);
		F.SetBoundaryConditions(Element(2),zAxis,dirichletWall,dirichletWall);
		F.SetBoundaryConditions(Element(5),zAxis,dirichletWall,dirichletWall);
	}
	F.SetBoundaryConditions(Element(6,7),zAxis,neumannWall,neumannWall);

	// transverse boundary conditions
	if (owner->bc0[1]==axisymmetric || owner->bc0[1]==reflecting)
	{
		F.SetBoundaryConditions(Element(0,5),xAxis, neumannWall, neumannWall);
		F.SetBoundaryConditions(Element(0),xAxis,dirichletWall,neumannWall);
		F.SetBoundaryConditions(Element(3),xAxis,dirichletWall,neumannWall);
	}
	else
	{
		F.SetBoundaryConditions(Element(0,5),xAxis, neumannWall, neumannWall);
	}
	if (owner->bc0[2]==axisymmetric || owner->bc0[2]==reflecting)
	{
		F.SetBoundaryConditions(Element(0,5),yAxis, neumannWall, neumannWall);
		F.SetBoundaryConditions(Element(1),yAxis,dirichletWall,neumannWall);
		F.SetBoundaryConditions(Element(4),yAxis,dirichletWall,neumannWall);
	}
	else
	{
		F.SetBoundaryConditions(Element(0,5),yAxis, neumannWall, neumannWall);
	}
	F.SetBoundaryConditions(Element(6,7),xAxis,neumannWall,neumannWall);
	F.SetBoundaryConditions(Element(6,7),yAxis,neumannWall,neumannWall);

	ComputeFinalFields();
}

void PGCSolver::MoveWindow()
{
	for (StripIterator s(*this,3,strongbool::yes);s<s.end();++s)
		F.Shift(s,-1,0.0);
	F.DownwardCopy(zAxis,1);
}

void PGCSolver::AntiMoveWindow()
{
	for (StripIterator s(*this,3,strongbool::yes);s<s.end();++s)
	{
		tw::Float polarizationFactor = polarizationType==circularPolarization ? 1.414 : 1.0;
		tw::Complex incoming0(0,0);
		tw::Complex incoming1(0,0);
		for (tw::Int p=0;p<owner->pulse.size();p++)
		{
			tw::vec3 pos = owner->Pos(s,0);
			pos.z = owner->ToLab(pos.z,-dth);
			incoming0 += polarizationFactor*owner->pulse[p]->VectorPotentialEnvelope(owner->elapsedTime-dth,pos);
			pos = owner->Pos(s,0);
			pos.z = owner->ToLab(pos.z,dth);
			incoming1 += polarizationFactor*owner->pulse[p]->VectorPotentialEnvelope(owner->elapsedTime+dth,pos);
		}
		a0.Shift(s,1,incoming0);
		a1.Shift(s,1,incoming1);
	}
	a0.UpwardCopy(zAxis,1);
	a1.UpwardCopy(zAxis,1);
}

void PGCSolver::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "LaserEnergy WaveAction ";
}

void PGCSolver::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	tw::Int i,j,k;
	tw::Int x0,x1,y0,y1,z0,z1;
	tw::Float fieldEnergy,waveAction;
	tw::vec3 pos;
	tw::Complex aNow,dtau,dzeta;
	tw::Complex eNow,bNow;

	fieldEnergy = 0.0;
	waveAction = 0.0;
	theRgn.GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (i=x0;i<=x1;i++)
		for (j=y0;j<=y1;j++)
			for (k=z0;k<=z1;k++)
			{
				pos = owner->Pos(i,j,k);
				if (theRgn.Inside(pos,*owner))
				{
					aNow = half*(a0(i,j,k) + a1(i,j,k));
					dtau = dti*(a1(i,j,k) - a0(i,j,k));
					dzeta = half*(a0(i,j,k+1) + a1(i,j,k+1) - a0(i,j,k-1) - a1(i,j,k-1)) / (owner->X(k+1,3) - owner->X(k-1,3));

					eNow = ii*laserFreq*aNow - (dtau-dzeta);
					bNow = ii*laserFreq*aNow + dzeta;
					fieldEnergy += half*(norm(eNow) + norm(bNow))*owner->dS(i,j,k,0);
					waveAction += imag( conj(aNow)*bNow - aNow*conj(bNow) ) * owner->dS(i,j,k,0);
				}
			}

	cols.push_back(fieldEnergy); avg.push_back(false);
	cols.push_back(waveAction); avg.push_back(false);
}

void PGCSolver::Update()
{
	tw::Int i,j,k;
	//chi.ImpliedDeposit();
	chi.DepositFromNeighbors();
	chi.DivideCellVolume(*owner);
	chi.ApplyBoundaryCondition();
	if (owner->smoothing>0)
		chi.Smooth(*owner,owner->smoothing,owner->compensation);
	for (i=lb[1];i<=ub[1];i++)
		for (j=lb[2];j<=ub[2];j++)
			for (k=1;k<=dim[3];k++)
				chi(i,j,k) = owner->ValueOnLightGrid<ComplexField,tw::Complex>(chi,i,j,k,dth);
	chi.DownwardCopy(zAxis,1);
	chi.UpwardCopy(zAxis,1);
	chi.ApplyBoundaryCondition();

	propagator->Advance(a0,a1,chi);
	ComputeFinalFields();
}

void PGCSolver::ComputeFinalFields()
{
	tw::Int i,j,k;

	for (i=lb[1];i<=ub[1];i++)
		for (j=lb[2];j<=ub[2];j++)
			for (k=1;k<=dim[3];k++)
			{
				F(i,j,k,7) = norm(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,i,j,k,dth));
				F(i,j,k,6) = norm(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,i,j,k,-dth));
				F(i,j,k,6) = 0.5*(F(i,j,k,6) + F(i,j,k,7));
			}

	F.DownwardCopy(zAxis,Element(6,7),1);
	F.UpwardCopy(zAxis,Element(6,7),1);

	for (i=1;i<=dim[1];i++)
		for (j=1;j<=dim[2];j++)
			for (k=1;k<=dim[3];k++)
			{
				F(i,j,k,0) = (F(i+1,j,k,6) - F(i-1,j,k,6)) / owner->dL(i,j,k,1);
				F(i,j,k,1) = (F(i,j+1,k,6) - F(i,j-1,k,6)) / owner->dL(i,j,k,2);
				F(i,j,k,2) = (F(i,j,k+1,6) - F(i,j,k-1,6)) / owner->dL(i,j,k,3);

				F(i,j,k,3) = (F(i+1,j,k,7) - F(i-1,j,k,7)) / owner->dL(i,j,k,1);
				F(i,j,k,4) = (F(i,j+1,k,7) - F(i,j-1,k,7)) / owner->dL(i,j,k,2);
				F(i,j,k,5) = (F(i,j,k+1,7) - F(i,j,k-1,7)) / owner->dL(i,j,k,3);
			}

	F.CopyFromNeighbors(Element(0,5));
	F.ApplyBoundaryCondition();
}

void PGCSolver::BoxDiagnose(GridDataDescriptor* box)
{
	ScalarField er,ei;
	er.Initialize(*this,owner);
	ei.Initialize(*this,owner);
	tw::Complex dadt,anow;
	tw::Int i,j,k;
	for (i=1;i<=dim[1];i++)
		for (j=1;j<=dim[2];j++)
			for (k=1;k<=dim[3];k++)
			{
				dadt = dti*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,i,j,k,dth) - owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,i,j,k,-dth));
				anow = tw::Float(0.5)*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,i,j,k,-dth) + owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,i,j,k,dth));
				er(i,j,k) = -real(dadt - ii*laserFreq*anow);
				ei(i,j,k) = -imag(dadt - ii*laserFreq*anow);
			}

	owner->WriteBoxData("e_real",box,&er(0,0,0),er.Stride());
	owner->WriteBoxData("e_imag",box,&ei(0,0,0),ei.Stride());
	owner->WriteBoxData("a2",box,&F(0,0,0,7),F.Stride());
	owner->WriteBoxData("chi_real",box,&chi(0,0,0,0),chi.Stride());
	owner->WriteBoxData("chi_imag",box,&chi(0,0,0,1),chi.Stride());
}
