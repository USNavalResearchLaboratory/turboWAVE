#include "simulation.h"
#include "fieldSolve.h"
#include "laserSolve.h"
using namespace tw::bc;


//////////////////////////////
//                          //
//     LASER SOLVER BASE    //
//                          //
//////////////////////////////


LaserSolver::LaserSolver(const std::string& name,Simulation* sim):Module(name,sim)
{
	updateSequencePriority = 3;
	laserFreq = 10.0;
	polarizationType = linearPolarization;
	propagator = NULL;

	a0.Initialize(*this,owner);
	a1.Initialize(*this,owner);
	chi.Initialize(*this,owner);

	directives.Add("carrier frequency",new tw::input::Float(&laserFreq));
	std::map<std::string,tw_polarization_type> pol = {{"linear",linearPolarization},{"circular",circularPolarization},{"radial",radialPolarization}};
	directives.Add("polarization",new tw::input::Enums<tw_polarization_type>(pol,&polarizationType));
}

LaserSolver::~LaserSolver()
{
	if (propagator!=NULL)
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

	for (auto cell : EntireCellRange(*this))
		for (auto pulse : owner->pulse)
		{
			pos = owner->Pos(cell);
			pos.z = owner->ToLab(pos.z,-dth);
			a0(cell) += polarizationFactor*pulse->VectorPotentialEnvelope(-dth,pos,laserFreq);
			pos = owner->Pos(cell);
			pos.z = owner->ToLab(pos.z,dth);
			a1(cell) += polarizationFactor*pulse->VectorPotentialEnvelope(dth,pos,laserFreq);
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

void LaserSolver::VerifyInput()
{
	Module::VerifyInput();
	for (auto tool : moduleTool)
	{
		propagator = dynamic_cast<LaserPropagator*>(tool);
		if (propagator!=NULL)
			break;
	}
	if (propagator==NULL)
		propagator = (LaserPropagator*)owner->CreateTool("default_adi",tw::tool_type::adiPropagator);
}

void LaserSolver::ReadData(std::ifstream& inFile)
{
	Module::ReadData(inFile);
	inFile.read((char *)&laserFreq,sizeof(tw::Float));
	inFile.read((char *)&polarizationType,sizeof(tw_polarization_type));
	a0.ReadData(inFile);
	a1.ReadData(inFile);
}

void LaserSolver::WriteData(std::ofstream& outFile)
{
	Module::WriteData(outFile);
	outFile.write((char *)&laserFreq,sizeof(tw::Float));
	outFile.write((char *)&polarizationType,sizeof(tw_polarization_type));
	a0.WriteData(outFile);
	a1.WriteData(outFile);
}

void LaserSolver::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader("a2",box);
	owner->WriteBoxDataHeader("e_real",box);
	owner->WriteBoxDataHeader("e_imag",box);
	owner->WriteBoxDataHeader("j_real",box);
	owner->WriteBoxDataHeader("j_imag",box);
	owner->WriteBoxDataHeader("a_real_raw",box);
	owner->WriteBoxDataHeader("a_imag_raw",box);
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


//////////////////////////////
//                          //
// Quasistatic Laser Solver //
//                          //
//////////////////////////////


QSSolver::QSSolver(const std::string& name,Simulation* sim):LaserSolver(name,sim)
{
	typeCode = tw::module_type::qsLaser;
}


//////////////////////////////
//                          //
//     PGC LASER SOLVER     //
//                          //
//////////////////////////////


PGCSolver::PGCSolver(const std::string& name,Simulation* sim):LaserSolver(name,sim)
{
	typeCode = tw::module_type::pgcLaser;

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
		F.SetBoundaryConditions(Element(0,5),tw::dom::zAxis,fld::neumannWall,fld::dirichletWall);
	}
	else
	{
		F.SetBoundaryConditions(Element(0,5),tw::dom::zAxis,fld::neumannWall,fld::neumannWall);
		F.SetBoundaryConditions(Element(2),tw::dom::zAxis,fld::dirichletWall,fld::dirichletWall);
		F.SetBoundaryConditions(Element(5),tw::dom::zAxis,fld::dirichletWall,fld::dirichletWall);
	}
	F.SetBoundaryConditions(Element(6,7),tw::dom::zAxis,fld::neumannWall,fld::neumannWall);

	// transverse boundary conditions
	if (owner->bc0[1]==par::axisymmetric || owner->bc0[1]==par::reflecting)
	{
		F.SetBoundaryConditions(Element(0,5),tw::dom::xAxis, fld::neumannWall, fld::neumannWall);
		F.SetBoundaryConditions(Element(0),tw::dom::xAxis,fld::dirichletWall,fld::neumannWall);
		F.SetBoundaryConditions(Element(3),tw::dom::xAxis,fld::dirichletWall,fld::neumannWall);
	}
	else
	{
		F.SetBoundaryConditions(Element(0,5),tw::dom::xAxis, fld::neumannWall, fld::neumannWall);
	}
	if (owner->bc0[2]==par::axisymmetric || owner->bc0[2]==par::reflecting)
	{
		F.SetBoundaryConditions(Element(0,5),tw::dom::yAxis, fld::neumannWall, fld::neumannWall);
		F.SetBoundaryConditions(Element(1),tw::dom::yAxis,fld::dirichletWall,fld::neumannWall);
		F.SetBoundaryConditions(Element(4),tw::dom::yAxis,fld::dirichletWall,fld::neumannWall);
	}
	else
	{
		F.SetBoundaryConditions(Element(0,5),tw::dom::yAxis, fld::neumannWall, fld::neumannWall);
	}
	F.SetBoundaryConditions(Element(6,7),tw::dom::xAxis,fld::neumannWall,fld::neumannWall);
	F.SetBoundaryConditions(Element(6,7),tw::dom::yAxis,fld::neumannWall,fld::neumannWall);

	ComputeFinalFields();
}

void PGCSolver::MoveWindow()
{
	for (auto s : StripRange(*this,3,strongbool::yes))
		F.Shift(s,-1,0.0);
	F.DownwardCopy(tw::dom::zAxis,1);
}

void PGCSolver::AntiMoveWindow()
{
	for (auto s : StripRange(*this,3,strongbool::yes))
	{
		tw::Float polarizationFactor = polarizationType==circularPolarization ? 1.414 : 1.0;
		tw::Complex incoming0(0,0);
		tw::Complex incoming1(0,0);
		for (tw::Int p=0;p<owner->pulse.size();p++)
		{
			tw::vec3 pos = owner->Pos(s,0);
			pos.z = owner->ToLab(pos.z,-dth);
			incoming0 += polarizationFactor*owner->pulse[p]->VectorPotentialEnvelope(owner->elapsedTime-dth,pos,laserFreq);
			pos = owner->Pos(s,0);
			pos.z = owner->ToLab(pos.z,dth);
			incoming1 += polarizationFactor*owner->pulse[p]->VectorPotentialEnvelope(owner->elapsedTime+dth,pos,laserFreq);
		}
		a0.Shift(s,1,incoming0);
		a1.Shift(s,1,incoming1);
	}
	a0.UpwardCopy(tw::dom::zAxis,1);
	a1.UpwardCopy(tw::dom::zAxis,1);
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
					fieldEnergy += 0.25*(norm(eNow) + norm(bNow))*owner->dS(i,j,k,0);
					waveAction += imag( conj(aNow)*bNow - aNow*conj(bNow) ) * owner->dS(i,j,k,0);
				}
			}

	cols.push_back(fieldEnergy); avg.push_back(false);
	cols.push_back(waveAction); avg.push_back(false);
}

void PGCSolver::Update()
{
	chi.DepositFromNeighbors();
	chi.ApplyFoldingCondition();
	chi.DivideCellVolume(*owner);
	chi.ApplyBoundaryCondition();
	chi.Smooth(*owner,smoothing,compensation);
	#pragma omp parallel
	{
		for (auto s : StripRange(*this,3,strongbool::yes))
		{
			for (tw::Int k=1;k<=dim[3];k++)
				chi(s,k) = owner->ValueOnLightGrid<ComplexField,tw::Complex>(chi,s,k,dth);
		}
	}
	chi.DownwardCopy(tw::dom::zAxis,1);
	chi.UpwardCopy(tw::dom::zAxis,1);
	chi.ApplyBoundaryCondition();

	propagator->Advance(a0,a1,chi);
	ComputeFinalFields();
}

void PGCSolver::ComputeFinalFields()
{
	#pragma omp parallel
	{
		for (auto s : StripRange(*this,3,strongbool::yes))
		{
			for (tw::Int k=1;k<=dim[3];k++)
			{
				F(s,k,7) = norm(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth));
				F(s,k,6) = norm(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth));
				F(s,k,6) = 0.5*(F(s,k,6) + F(s,k,7));
			}
		}
	}

	F.DownwardCopy(tw::dom::zAxis,Element(6,7),1);
	F.UpwardCopy(tw::dom::zAxis,Element(6,7),1);

	#pragma omp parallel
	{
		for (auto cell : InteriorCellRange(*this))
		{
			F(cell,0) = F(cell,6,1);
			F(cell,1) = F(cell,6,2);
			F(cell,2) = F(cell,6,3);

			F(cell,3) = F(cell,7,1);
			F(cell,4) = F(cell,7,2);
			F(cell,5) = F(cell,7,3);
		}
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
	for (auto s : StripRange(*this,3,strongbool::no))
	{
		for (tw::Int k=1;k<=dim[3];k++)
		{
			dadt = dti*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth) - owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth));
			anow = tw::Float(0.5)*(owner->ValueOnLabGrid<ComplexField,tw::Complex>(a0,s,k,-dth) + owner->ValueOnLabGrid<ComplexField,tw::Complex>(a1,s,k,dth));
			er(s,k) = -real(dadt - ii*laserFreq*anow);
			ei(s,k) = -imag(dadt - ii*laserFreq*anow);
		}
	}

	owner->WriteBoxData("a_real_raw",box,&a1(0,0,0,0),a1.Stride());
	owner->WriteBoxData("a_imag_raw",box,&a1(0,0,0,1),a1.Stride());
	owner->WriteBoxData("e_real",box,&er(0,0,0),er.Stride());
	owner->WriteBoxData("e_imag",box,&ei(0,0,0),ei.Stride());
	owner->WriteBoxData("a2",box,&F(0,0,0,7),F.Stride());
	owner->WriteBoxData("j_real",box,&chi(0,0,0,0),chi.Stride());
	owner->WriteBoxData("j_imag",box,&chi(0,0,0,1),chi.Stride());
}
