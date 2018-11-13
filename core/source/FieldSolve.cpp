#include "sim.h"
#include "fieldSolve.h"

/////////////////////////////
// FIELD SOLVER BASE CLASS //
/////////////////////////////


FieldSolver::FieldSolver(Grid* theGrid):Module(theGrid)
{
	updateSequencePriority = 3;
	ellipticSolver = NULL;
}

FieldSolver::~FieldSolver()
{
	if (ellipticSolver!=NULL)
		owner->RemoveTool(ellipticSolver);
}

void FieldSolver::Initialize()
{
	Module::Initialize();
	if (ellipticSolver==NULL)
		ellipticSolver = (EllipticSolver*)owner->CreateTool("facr_poisson_solver",facrPoissonSolver);
}

void FieldSolver::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;

	Module::ReadInputFileDirective(inputString,command);

	ellipticSolver = (EllipticSolver*)owner->ToolFromDirective(inputString,command);

	if (ellipticSolver!=NULL)
		ellipticSolver->ReadInputFileDirective(inputString,command);
}

void FieldSolver::ReadData(std::ifstream& inFile)
{
	Module::ReadData(inFile);
	ellipticSolver = (EllipticSolver*)owner->LoadRestartedTool(inFile);
}

void FieldSolver::WriteData(std::ofstream& outFile)
{
	Module::WriteData(outFile);
	ellipticSolver->SaveToolReference(outFile);
}


////////////////////////////////
// ELECTROMAGNETIC BASE CLASS //
////////////////////////////////


Electromagnetic::Electromagnetic(Grid* theGrid):FieldSolver(theGrid)
{
	if (1.0/theGrid->dt <= sqrt(
		(theGrid->globalCells[1]==1 ? 0.0 : 1.0)/sqr(dx(*theGrid)) +
		(theGrid->globalCells[2]==1 ? 0.0 : 1.0)/sqr(dy(*theGrid)) +
		(theGrid->globalCells[3]==1 ? 0.0 : 1.0)/sqr(dz(*theGrid))))
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

	sources.SetBoundaryConditions(xAxis,dirichletCell,dirichletCell);
	sources.SetBoundaryConditions(Element(1),xAxis,normalFluxFixed,normalFluxFixed);

	sources.SetBoundaryConditions(yAxis,dirichletCell,dirichletCell);
	sources.SetBoundaryConditions(Element(2),yAxis,normalFluxFixed,normalFluxFixed);

	sources.SetBoundaryConditions(zAxis,dirichletCell,dirichletCell);
	sources.SetBoundaryConditions(Element(3),zAxis,normalFluxFixed,normalFluxFixed);

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
	if (owner->smoothing>0)
		sources.Smooth(*owner,owner->smoothing,owner->compensation);
}
#endif

void Electromagnetic::MoveWindow()
{
	for (StripIterator s(*this,3,strongbool::yes);s<s.end();++s)
		F.Shift(s,-1,0.0);
	F.DownwardCopy(zAxis,1);
}

void Electromagnetic::InitializeConductors()
{
	// Conductor resides in cells shifted back by 1/2
	// These have E known along upper edges and B normal to upper walls
	// The conductor fills the whole cell or none of the cell
	tw::vec3 shiftedCenter;
	//#pragma omp parallel for private(i,j,k,s,shiftedCenter) collapse(3) schedule(static)
	for (CellIterator cell(*this,true);cell<cell.end();++cell)
	{
		conductorMask(cell) = 1.0;
		shiftedCenter = owner->Pos(cell) - 0.5*owner->dPos(cell);
		for (tw::Int s=0;s<owner->conductor.size();s++)
			if (owner->conductor[s]->affectsA && owner->conductor[s]->theRgn->Inside(shiftedCenter,*owner))
				conductorMask(cell) = 0.0;
	}
}

void Electromagnetic::SetExteriorBoundaryConditionsE(Field& A,const Element& ex,const Element& ey,const Element& ez)
{
	// in the following, it doesn't hurt to zero out low-side cell walls in low-side ghost cells
	A.SetBoundaryConditions(ex,xAxis,dirichletCell,none);
	A.SetBoundaryConditions(ex,yAxis,dirichletCell,dirichletCell);
	A.SetBoundaryConditions(ex,zAxis,dirichletCell,dirichletCell);

	A.SetBoundaryConditions(ey,xAxis,dirichletCell,dirichletCell);
	A.SetBoundaryConditions(ey,yAxis,dirichletCell,none);
	A.SetBoundaryConditions(ey,zAxis,dirichletCell,dirichletCell);

	A.SetBoundaryConditions(ez,xAxis,dirichletCell,dirichletCell);
	A.SetBoundaryConditions(ez,yAxis,dirichletCell,dirichletCell);
	A.SetBoundaryConditions(ez,zAxis,dirichletCell,none);
}

void Electromagnetic::SetExteriorBoundaryConditionsB(Field& A,const Element& bx,const Element& by,const Element& bz)
{
	// in the following, it doesn't hurt to zero out low-side cell walls in low-side ghost cells
	A.SetBoundaryConditions(bx,xAxis,dirichletCell,dirichletCell);
	A.SetBoundaryConditions(bx,yAxis,dirichletCell,none);
	A.SetBoundaryConditions(bx,zAxis,dirichletCell,none);

	A.SetBoundaryConditions(by,xAxis,dirichletCell,none);
	A.SetBoundaryConditions(by,yAxis,dirichletCell,dirichletCell);
	A.SetBoundaryConditions(by,zAxis,dirichletCell,none);

	A.SetBoundaryConditions(bz,xAxis,dirichletCell,none);
	A.SetBoundaryConditions(bz,yAxis,dirichletCell,none);
	A.SetBoundaryConditions(bz,zAxis,dirichletCell,dirichletCell);
}

void Electromagnetic::ForceQuasistaticVectorPotential(Field& A4,ScalarField& DtPhi)
{
	// only works in cartesian, assumes coulomb gauge
	tw::Int i,j,k,ax;
	tw::Float beta = sqrt(1.0 - 1.0/sqr(gammaBeam));
	tw::Float w = beta*dt*freq.z;
	VectorizingIterator<3> v(*this,true);
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
				v.SetStrip(i,j,0);
				for (k=1;k<=dim[3];k++)
					rhs(v,k,0) = DtPhi(v,k,0,ax) - sources(v,k,ax);
			}
		rhs.CopyFromNeighbors();

		ellipticSolver->Solve(ans,rhs,1.0);
		for (CellIterator cell(*this,true);cell<cell.end();++cell)
			A4(cell,ax+4) = ans(cell);

		// Find Solution at last time

		for (i=1;i<=dim[1];i++)
			for (j=1;j<=dim[2];j++)
			{
				v.SetStrip(i,j,0);
				for (k=1;k<=dim[3];k++)
					rhs(v,k,0) = DtPhi(v,k,0,ax) - sources(v,k,ax);
			}
		rhs.CopyFromNeighbors();

		for (k=lb[3];k<=dim[3];k++)
			for (j=lb[2];j<=ub[2];j++)
				for (i=lb[1];i<=ub[1];i++)
					rhs(i,j,k) = (1.0 - w)*rhs(i,j,k) + w*rhs(i,j,k+1);
		rhs.DownwardCopy(zAxis,1);

		ellipticSolver->Solve(ans,rhs,1.0);
		for (CellIterator cell(*this,true);cell<cell.end();++cell)
			A4(cell,ax) = ans(cell);
	}

	ellipticSolver->gammaBeam = 1.0;
}

void Electromagnetic::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;

	FieldSolver::ReadInputFileDirective(inputString,command);

	if (command=="dipole") // eg, dipole center = 0 0 0
	{
		inputString >> word >> word >> dipoleCenter.x >> dipoleCenter.y >> dipoleCenter.z;
	}
	if (command=="gamma") // eg, gamma beam = 10.0
	{
		inputString >> word >> word >> gammaBeam;
	}
	if (command=="far-field") // eg, new far-field diagnostic { ... }
	{
		inputString >> word	>> word;
		farFieldDetector.push_back(new FarFieldDetectorDescriptor(owner));
		farFieldDetector.back()->ReadInputFile(inputString);
	}
}

void Electromagnetic::ReadData(std::ifstream& inFile)
{
	tw::Int i,num;

	FieldSolver::ReadData(inFile);
	inFile.read((char*)&dipoleCenter,sizeof(tw::vec3));
	inFile.read((char*)&gammaBeam,sizeof(tw::Float));
	F.ReadData(inFile);

	inFile.read((char *)&num,sizeof(tw::Int));
	(*owner->tw_out) << "Read " << num << " far-field diagnostics" << std::endl;
	for (i=0;i<num;i++)
	{
		farFieldDetector.push_back(new FarFieldDetectorDescriptor(owner));
		farFieldDetector.back()->ReadData(inFile);
	}
}

void Electromagnetic::WriteData(std::ofstream& outFile)
{
	tw::Int i;

	FieldSolver::WriteData(outFile);
	outFile.write((char*)&dipoleCenter,sizeof(tw::vec3));
	outFile.write((char*)&gammaBeam,sizeof(tw::Float));
	F.WriteData(outFile);

	i = farFieldDetector.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<farFieldDetector.size();i++)
		farFieldDetector[i]->WriteData(outFile);
}

void Electromagnetic::StartDiagnostics()
{
	#ifdef USE_OPENCL
	F.ReceiveFromComputeBuffer();
	sources.ReceiveFromComputeBuffer();
	#endif
}

void Electromagnetic::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "FieldEnergy J.E ExB.dS TotalCharge Current Px Py Pz ";
}

void Electromagnetic::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	tw::Int i,j,k;
	tw::Float fieldEnergy,JE,fluxOut,totalCharge,current,dV;
	tw::vec3 dipoleMoment,pos,PoyntingVector,r1,r2;
	tw::vec3 Ef,Bf,J;

	fieldEnergy = 0.0;
	JE = 0.0;
	fluxOut = 0.0;
	totalCharge = 0.0;
	current = 0.0;

	tw::Int xl,xh,yl,yh,zl,zh;
	tw::Int x0,x1,y0,y1,z0,z1;

	theRgn.GetRawCellBounds(&xl,&xh,&yl,&yh,&zl,&zh);
	theRgn.GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (k=z0;k<=z1;k++)
		for (j=y0;j<=y1;j++)
			for (i=x0;i<=x1;i++)
			{
				pos = owner->Pos(i,j,k);
				if (theRgn.Inside(pos,*owner))
				{
					dV = owner->dS(i,j,k,0);
					Ef = F.Vec3(i,j,k,0);
					Bf = F.Vec3(i,j,k,3);
					J = sources.Vec3(i,j,k,1);
					owner->CurvilinearToCartesian(&(r1=dipoleCenter));
					owner->CurvilinearToCartesian(&(r2=pos));
					dipoleMoment += (r2 - r1) * sources(i,j,k,0) * dV;
					fieldEnergy += 0.5 * (Norm(Ef) + Norm(Bf)) * dV;
					JE += J ^ Ef * dV;
					totalCharge += sources(i,j,k,0) * dV;
					current += owner->dS(i,j,k,3) * 0.5 * (sources(i,j,k-1,3) + sources(i,j,k,3)) / tw::Float(z1 - z0 + 1);
				}
			}

	if (x0<=x1 && y0<=y1 && z0<=z1) // region intersects this domain
	{
		for (k=z0;k<=z1;k++)
			for (j=y0;j<=y1;j++)
			{
				pos = owner->Pos(x0,j,k);
				if (theRgn.Inside(pos,*owner) && xl>=1)
				{
					PoyntingVector = F.Vec3(x0,j,k,0) | F.Vec3(x0,j,k,3);
					fluxOut -= owner->dS(x0,j,k,1) * PoyntingVector.x;
				}
				pos = owner->Pos(x1,j,k);
				if (theRgn.Inside(pos,*owner) && xh<=dim[1])
				{
					PoyntingVector = F.Vec3(x1,j,k,0) | F.Vec3(x1,j,k,3);
					fluxOut += owner->dS(x1,j,k,1) * PoyntingVector.x;
				}
			}
		for (k=z0;k<=z1;k++)
			for (i=x0;i<=x1;i++)
			{
				pos = owner->Pos(i,y0,k);
				if (theRgn.Inside(pos,*owner) && yl>=1)
				{
					PoyntingVector = F.Vec3(i,y0,k,0) | F.Vec3(i,y0,k,3);
					fluxOut -= owner->dS(i,y0,k,2) * PoyntingVector.y;
				}
			   pos = owner->Pos(i,y1,k);
				if (theRgn.Inside(pos,*owner) && yh<=dim[2])
				{
					PoyntingVector = F.Vec3(i,y1,k,0) | F.Vec3(i,y1,k,3);
					fluxOut += owner->dS(i,y1,k,2) * PoyntingVector.y;
				}
			}
		for (j=y0;j<=y1;j++)
			for (i=x0;i<=x1;i++)
			{
				pos = owner->Pos(i,j,z0);
				if (theRgn.Inside(pos,*owner) && zl>=1)
				{
					PoyntingVector = F.Vec3(i,j,z0,0) | F.Vec3(i,j,z0,3);
					fluxOut -= owner->dS(i,j,z0,3) * PoyntingVector.z;
				}
				pos = owner->Pos(i,j,z1);
				if (theRgn.Inside(pos,*owner) && zh<=dim[3])
				{
					PoyntingVector = F.Vec3(i,j,z1,0) | F.Vec3(i,j,z1,3);
					fluxOut += owner->dS(i,j,z1,3) * PoyntingVector.z;
				}
			}
	}

	cols.push_back(fieldEnergy); avg.push_back(false);
	cols.push_back(JE); avg.push_back(false);
	cols.push_back(fluxOut); avg.push_back(false);
	cols.push_back(totalCharge); avg.push_back(false);
	cols.push_back(current); avg.push_back(false);
	cols.push_back(dipoleMoment.x); avg.push_back(false);
	cols.push_back(dipoleMoment.y); avg.push_back(false);
	cols.push_back(dipoleMoment.z); avg.push_back(false);
}

void Electromagnetic::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader("Ex",box);
	owner->WriteBoxDataHeader("Ey",box);
	owner->WriteBoxDataHeader("Ez",box);
	owner->WriteBoxDataHeader("Bx",box);
	owner->WriteBoxDataHeader("By",box);
	owner->WriteBoxDataHeader("Bz",box);
	owner->WriteBoxDataHeader("rho",box);
	owner->WriteBoxDataHeader("Jx",box);
	owner->WriteBoxDataHeader("Jy",box);
	owner->WriteBoxDataHeader("Jz",box);
}

void Electromagnetic::BoxDiagnose(GridDataDescriptor* box)
{
	owner->WriteBoxData("Ex",box,&F(0,0,0,0),F.Stride());
	owner->WriteBoxData("Ey",box,&F(0,0,0,1),F.Stride());
	owner->WriteBoxData("Ez",box,&F(0,0,0,2),F.Stride());
	owner->WriteBoxData("Bx",box,&F(0,0,0,3),F.Stride());
	owner->WriteBoxData("By",box,&F(0,0,0,4),F.Stride());
	owner->WriteBoxData("Bz",box,&F(0,0,0,5),F.Stride());
	owner->WriteBoxData("rho",box,&sources(0,0,0,0),sources.Stride());
	owner->WriteBoxData("Jx",box,&sources(0,0,0,1),sources.Stride());
	owner->WriteBoxData("Jy",box,&sources(0,0,0,2),sources.Stride());
	owner->WriteBoxData("Jz",box,&sources(0,0,0,3),sources.Stride());
}

void Electromagnetic::PointDiagnosticHeader(std::ofstream& outFile)
{
	outFile << "Ex Ey Ez Bx By Bz ";
}

void Electromagnetic::PointDiagnose(std::ofstream& outFile,const weights_3D& w)
{
	std::valarray<tw::Float> EM(6);
	F.Interpolate(EM,w);
	for (tw::Int i=0;i<6;i++)
		outFile << EM[i] << " ";
}

void Electromagnetic::CustomDiagnose()
{
	tw::Int i,j,k,s;

	std::stringstream fileName;
	tw::Int master = 0;

	if (owner->IsFirstStep() && !(owner->appendMode&&owner->restarted) && owner->strip[0].Get_rank()==master)
	{
		for (s=0;s<farFieldDetector.size();s++)
		{
			fileName.str("");
			fileName << farFieldDetector[s]->filename << "-Atheta.dvdat";
			farFieldDetector[s]->AthetaFile.open(fileName.str().c_str(),std::ios::binary);
			WriteDVHeader(farFieldDetector[s]->AthetaFile,2,
				farFieldDetector[s]->timePts,farFieldDetector[s]->thetaPts,farFieldDetector[s]->phiPts,
				farFieldDetector[s]->t0,farFieldDetector[s]->t1,
				farFieldDetector[s]->theta0,farFieldDetector[s]->theta1,
				farFieldDetector[s]->phi0,farFieldDetector[s]->phi1);
			farFieldDetector[s]->AthetaFile.close();

			fileName.str("");
			fileName << farFieldDetector[s]->filename << "-Aphi.dvdat";
			farFieldDetector[s]->AphiFile.open(fileName.str().c_str(),std::ios::binary);
			WriteDVHeader(farFieldDetector[s]->AphiFile,2,
				farFieldDetector[s]->timePts,farFieldDetector[s]->thetaPts,farFieldDetector[s]->phiPts,
				farFieldDetector[s]->t0,farFieldDetector[s]->t1,
				farFieldDetector[s]->theta0,farFieldDetector[s]->theta1,
				farFieldDetector[s]->phi0,farFieldDetector[s]->phi1);
			farFieldDetector[s]->AphiFile.close();
		}
	}

	for (s=0;s<farFieldDetector.size();s++)
	{
		if (farFieldDetector[s]->WriteThisStep(owner->elapsedTime,owner->dt,owner->stepNow))
			farFieldDetector[s]->AccumulateField(sources);
	}

	if (owner->IsLastStep())
	{
		Vec3Field accum;
		tw::vec3 ACG,rVec,thetaVec,phiVec;
		tw::Float theta,phi,dtheta,dphi;
		FarFieldDetectorDescriptor *det;
		float data;

		for (s=0;s<farFieldDetector.size();s++)
		{
			det = farFieldDetector[s];
			dtheta = (det->theta1-det->theta0)/tw::Float(det->thetaPts);
			dphi = (det->phi1-det->phi0)/tw::Float(det->phiPts);
			DiscreteSpace ff_layout(det->timePts,det->thetaPts,det->phiPts,tw::vec3(0.0,0.0,0.0),tw::vec3(det->t1-det->t0,det->theta1-det->theta0,det->phi1-det->phi0),1);
			accum.Initialize(ff_layout,owner);
			owner->strip[0].Sum(&det->A(0,0,0),&accum(0,0,0),sizeof(tw::vec3)*accum.TotalCells(),master);

			if (owner->strip[0].Get_rank()==master)
			{
				// put vector potential in coulomb gauge and spherical coordinates
				for (k=1;k<=det->phiPts;k++)
					for (j=1;j<=det->thetaPts;j++)
						for (i=1;i<=det->timePts;i++)
						{
							theta = det->theta0 + (tw::Float(j)-0.5)*dtheta;
							phi = det->phi0 + (tw::Float(k)-0.5)*dphi;
							rVec = tw::vec3( sin(theta)*cos(phi) , sin(theta)*sin(phi) , cos(theta) );
							thetaVec = tw::vec3( cos(theta)*cos(phi) , cos(theta)*sin(phi) , -sin(theta) );
							phiVec = tw::vec3( -sin(phi) , cos(phi) , 0.0 );
							ACG = rVec | (accum(i,j,k) | rVec); // form coulomb gauge A
							accum(i,j,k) = tw::vec3( ACG^rVec , ACG^thetaVec , ACG^phiVec ); // put in spherical coordinates
						}

				fileName.str("");
				fileName << det->filename << "-Atheta.dvdat";
				det->AthetaFile.open(fileName.str().c_str(),std::ios::binary | std::ios::app);
				for (k=1;k<=det->phiPts;k++)
					for (j=1;j<=det->thetaPts;j++)
						for (i=1;i<=det->timePts;i++)
						{
							data = accum(i,j,k).y;
							WriteBigEndian((char *)&data,sizeof(float),0,det->AthetaFile);
						}
				det->AthetaFile.close();

				fileName.str("");
				fileName << det->filename << "-Aphi.dvdat";
				det->AphiFile.open(fileName.str().c_str(),std::ios::binary | std::ios::app);
				for (k=1;k<=det->phiPts;k++)
					for (j=1;j<=det->thetaPts;j++)
						for (i=1;i<=det->timePts;i++)
						{
							data = accum(i,j,k).z;
							WriteBigEndian((char *)&data,sizeof(float),0,det->AphiFile);
						}
				det->AphiFile.close();
			}
		}
	}
}



////////////////////////////////////////////////
// COULOMB GAUGE FIELD SOLVER (original WAVE) //
////////////////////////////////////////////////


CoulombSolver::CoulombSolver(Grid* theGrid):Electromagnetic(theGrid)
{
	name = "Coulomb-EM";
	typeCode = tw::module_type::coulombSolver;
	A4.Initialize(8,*this,owner);

	L1.Initialize(theGrid,theGrid,&theGrid->wave,zAxis,lowSide,1,5);
	L2.Initialize(theGrid,theGrid,&theGrid->wave,zAxis,lowSide,2,6);
	R1.Initialize(theGrid,theGrid,&theGrid->wave,zAxis,highSide,1,5);
	R2.Initialize(theGrid,theGrid,&theGrid->wave,zAxis,highSide,2,6);

	z0 = natural;
	z1 = natural;
}

void CoulombSolver::ExchangeResources()
{
	Electromagnetic::ExchangeResources();
}

void CoulombSolver::Initialize()
{
	Electromagnetic::Initialize();

	// Overrides the boundary conditions set by the tool
	ellipticSolver->SetBoundaryConditions(periodic,periodic,periodic,periodic,z0,z1);
	ellipticSolver->SetBoundaryConditions(scratch2);

	SetExteriorBoundaryConditionsE(A4,Element(1),Element(2),Element(3));
	SetExteriorBoundaryConditionsE(A4,Element(5),Element(6),Element(7));

	if (owner->restarted)
		return;

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

void CoulombSolver::ReadData(std::ifstream& inFile)
{
	Electromagnetic::ReadData(inFile);

	A4.ReadData(inFile);
	L1.ReadData(inFile);
	L2.ReadData(inFile);
	R1.ReadData(inFile);
	R2.ReadData(inFile);
}

void CoulombSolver::WriteData(std::ofstream& outFile)
{
	Electromagnetic::WriteData(outFile);

	A4.WriteData(outFile);
	L1.WriteData(outFile);
	L2.WriteData(outFile);
	R1.WriteData(outFile);
	R2.WriteData(outFile);
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
		for (VectorizingIterator<3> v(*this,false);v<v.end();++v)
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
		for (StripIterator strip(A4,3,strongbool::no);strip<strip.end();++strip)
			A4(strip,0,3) = A4(strip,1,3) + spacing.z*(A4.dfwd(strip,0,1,1) + A4.dfwd(strip,0,2,2));
	}

	if (owner->n1[3]==MPI_PROC_NULL)
	{
		R1.Set(A4,owner->elapsedTime+dth,dt);
		R2.Set(A4,owner->elapsedTime+dth,dt);
		for (StripIterator strip(A4,3,strongbool::no);strip<strip.end();++strip)
			A4(strip,dim[3]+1,3) = A4(strip,dim[3],3) - spacing.z*(A4.dfwd(strip,dim[3],1,1) + A4.dfwd(strip,dim[3],2,2));
	}

	A4.DownwardCopy(xAxis,Element(1,3),1);
	A4.UpwardCopy(xAxis,Element(1,3),1);
	A4.DownwardCopy(yAxis,Element(1,3),1);
	A4.UpwardCopy(yAxis,Element(1,3),1);

	// Swap A0 and A1 so A1 has most recent data

	A4.Swap(Element(1,3),Element(5,7));

	// COMPUTE E AND B FIELDS FOR PUSHER

	ComputeFinalFields();
}

void CoulombSolver::ComputeFinalFields()
{
	#pragma omp parallel
	{
		for (VectorizingIterator<3> v(*this,false);v<v.end();++v)
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

void CoulombSolver::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	Electromagnetic::BoxDiagnosticHeader(box);

	owner->WriteBoxDataHeader("phi",box);
	owner->WriteBoxDataHeader("Ax",box);
	owner->WriteBoxDataHeader("Ay",box);
	owner->WriteBoxDataHeader("Az",box);

	owner->WriteBoxDataHeader("AxDot",box);
	owner->WriteBoxDataHeader("AyDot",box);
	owner->WriteBoxDataHeader("AzDot",box);
}

void CoulombSolver::BoxDiagnose(GridDataDescriptor* box)
{
	// Upon entry : A0(-1/2) , A1(1/2) , phi0(-1) , phi1(0)
	Electromagnetic::BoxDiagnose(box);

	owner->WriteBoxData("phi",box,&A4(0,0,0,4),A4.Stride());

	scratch1 = 0.0;
	AddMulFieldData(scratch1,Element(0),A4,Element(1),0.5);
	AddMulFieldData(scratch1,Element(0),A4,Element(5),0.5);
	owner->WriteBoxData("Ax",box,&scratch1(0,0,0),scratch1.Stride());

	scratch1 = 0.0;
	AddMulFieldData(scratch1,Element(0),A4,Element(2),0.5);
	AddMulFieldData(scratch1,Element(0),A4,Element(6),0.5);
	owner->WriteBoxData("Ay",box,&scratch1(0,0,0),scratch1.Stride());

	scratch1 = 0.0;
	AddMulFieldData(scratch1,Element(0),A4,Element(3),0.5);
	AddMulFieldData(scratch1,Element(0),A4,Element(7),0.5);
	owner->WriteBoxData("Az",box,&scratch1(0,0,0),scratch1.Stride());

	scratch1 = 0.0;
	AddMulFieldData(scratch1,Element(0),A4,Element(1),-1.0/dt);
	AddMulFieldData(scratch1,Element(0),A4,Element(5),1.0/dt);
	owner->WriteBoxData("AxDot",box,&scratch1(0,0,0),scratch1.Stride());

	scratch1 = 0.0;
	AddMulFieldData(scratch1,Element(0),A4,Element(2),-1.0/dt);
	AddMulFieldData(scratch1,Element(0),A4,Element(6),1.0/dt);
	owner->WriteBoxData("AyDot",box,&scratch1(0,0,0),scratch1.Stride());

	scratch1 = 0.0;
	AddMulFieldData(scratch1,Element(0),A4,Element(3),-1.0/dt);
	AddMulFieldData(scratch1,Element(0),A4,Element(7),1.0/dt);
	owner->WriteBoxData("AzDot",box,&scratch1(0,0,0),scratch1.Stride());
}

void CoulombSolver::PointDiagnosticHeader(std::ofstream& outFile)
{
	Electromagnetic::PointDiagnosticHeader(outFile);
	outFile << "phi ";
}

void CoulombSolver::PointDiagnose(std::ofstream& outFile,const weights_3D& w)
{
	Electromagnetic::PointDiagnose(outFile,w);

	std::valarray<tw::Float> A4Now(8);

	A4.Interpolate(A4Now,w);
	outFile << A4Now[4] << " ";
}



///////////////////////
//  DIRECT SOLVER    //
//  ( with PML's )   //
///////////////////////


DirectSolver::DirectSolver(Grid* theGrid):Electromagnetic(theGrid)
{
	name = "Direct-EM";
	typeCode = tw::module_type::directSolver;
	yeeTool = (YeePropagatorPML*)owner->AddPrivateTool(yeePropagatorPML);
	enforceChargeConservation = true;
	layerThicknessX0 = 0;
	layerThicknessX1 = 0;
	layerThicknessY0 = 0;
	layerThicknessY1 = 0;
	layerThicknessZ0 = 0;
	layerThicknessZ1 = 0;
	reflectionCoefficient = 0.01;
	A.Initialize(12,*this,owner);
	DiscreteSpace pml_layout;
	pml_layout.Resize(dim[1],1,1,corner,size);
	PMLx.Initialize(6,pml_layout,owner); // sx,tx,jx,sxstar,txstar,jxstar
	pml_layout.Resize(dim[2],1,1,corner,size);
	PMLy.Initialize(6,pml_layout,owner);
	pml_layout.Resize(dim[3],1,1,corner,size);
	PMLz.Initialize(6,pml_layout,owner);

	#ifdef USE_OPENCL
	A.InitializeComputeBuffer();
	PMLx.InitializeComputeBuffer();
	PMLy.InitializeComputeBuffer();
	PMLz.InitializeComputeBuffer();
	yeeTool->SetupComputeKernels(F,A,PMLx,PMLy,PMLz,sources);
	#endif
}

DirectSolver::~DirectSolver()
{
	owner->RemoveTool(yeeTool);
}

void DirectSolver::SetupPML(Field& pml,tw::Int g0,tw::Int gN,tw::Int L0,tw::Int L1,tw::Float R,tw::Float ds)
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

	for (i=pml.N0(1);i<=pml.N1(1);i++)
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
		maxConductivity = -1.5*log(reflectionCoefficient)/delta;
		p0 = delta - 0.5*ds;
		for (i=pml.N0(1);i<=pml.N1(1);i++)
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
		maxConductivity = -1.5*log(reflectionCoefficient)/delta;
		p0 = ds*N1 - delta - 0.5*ds;
		for (i=pml.N0(1);i<=pml.N1(1);i++)
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

	if (owner->restarted)
		return;

	// Setup PML media

	SetupPML(PMLx,owner->GlobalCellIndex(0,1),owner->globalCells[1],layerThicknessX0,layerThicknessX1,reflectionCoefficient,spacing.x);
	SetupPML(PMLy,owner->GlobalCellIndex(0,2),owner->globalCells[2],layerThicknessY0,layerThicknessY1,reflectionCoefficient,spacing.y);
	SetupPML(PMLz,owner->GlobalCellIndex(0,3),owner->globalCells[3],layerThicknessZ0,layerThicknessZ1,reflectionCoefficient,spacing.z);

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
		for (VectorizingIterator<3> v(*this,true);v<v.end();++v)
		{
			#pragma omp simd
			for (tw::Int k=lb[3];k<=ub[3];k++)
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
		for (VectorizingIterator<3> v(*this,true);v<v.end();++v)
		{
			#pragma omp simd
			for (tw::Int k=lb[3];k<=ub[3];k++)
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
	for (StripIterator s(*this,3,strongbool::yes);s<s.end();++s)
		A.Shift(s,-1,0.0);
	A.DownwardCopy(zAxis,1);
}

void DirectSolver::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;

	Electromagnetic::ReadInputFileDirective(inputString,command);
	if (command=="enforce") // eg, enforce charge conservation = yes
	{
		inputString >> word >> word >> word >> word;
		enforceChargeConservation = (word=="yes" || word=="true" || word=="on");
	}
	if (command=="layer") // eg, layer thickness = 8
	{
		inputString >> word >> word >> layerThicknessX0;
		layerThicknessX1 = layerThicknessY0 = layerThicknessY1 = layerThicknessZ0 = layerThicknessZ1 = layerThicknessX0;
		if (owner->movingWindow)
			layerThicknessZ0 = layerThicknessZ1 = 0;
	}
	if (command=="layers") // eg, layers = ( 4 , 4 , 4 , 4 , 8 , 8 )
	{
		inputString >> word >> layerThicknessX0 >> layerThicknessX1;
		inputString >> layerThicknessY0 >> layerThicknessY1;
		inputString >> layerThicknessZ0 >> layerThicknessZ1;
	}
	if (command=="reflection") // eg, reflection coefficient = 0.01
	{
		inputString >> word >> word >> reflectionCoefficient;
	}
}

void DirectSolver::ReadData(std::ifstream& inFile)
{
	Electromagnetic::ReadData(inFile);
	A.ReadData(inFile);
	PMLx.ReadData(inFile);
	PMLy.ReadData(inFile);
	PMLz.ReadData(inFile);
	inFile.read((char *)&enforceChargeConservation,sizeof(bool));
	inFile.read((char *)&layerThicknessX0,sizeof(tw::Int));
	inFile.read((char *)&layerThicknessX1,sizeof(tw::Int));
	inFile.read((char *)&layerThicknessY0,sizeof(tw::Int));
	inFile.read((char *)&layerThicknessY1,sizeof(tw::Int));
	inFile.read((char *)&layerThicknessZ0,sizeof(tw::Int));
	inFile.read((char *)&layerThicknessZ1,sizeof(tw::Int));
	inFile.read((char *)&reflectionCoefficient,sizeof(tw::Float));
}

void DirectSolver::WriteData(std::ofstream& outFile)
{
	Electromagnetic::WriteData(outFile);
	A.WriteData(outFile);
	PMLx.WriteData(outFile);
	PMLy.WriteData(outFile);
	PMLz.WriteData(outFile);
	outFile.write((char *)&enforceChargeConservation,sizeof(bool));
	outFile.write((char *)&layerThicknessX0,sizeof(tw::Int));
	outFile.write((char *)&layerThicknessX1,sizeof(tw::Int));
	outFile.write((char *)&layerThicknessY0,sizeof(tw::Int));
	outFile.write((char *)&layerThicknessY1,sizeof(tw::Int));
	outFile.write((char *)&layerThicknessZ0,sizeof(tw::Int));
	outFile.write((char *)&layerThicknessZ1,sizeof(tw::Int));
	outFile.write((char *)&reflectionCoefficient,sizeof(tw::Float));
}

void DirectSolver::Update()
{
	tw::Int r;

	// Add electric antenna currents

	if (owner->conductor.size())
	{
		#ifdef USE_OPENCL
		sources.ReceiveFromComputeBuffer();
		#endif
		for (r=0;r<owner->conductor.size();r++)
			if (owner->conductor[r]->electricCurrent)
				owner->conductor[r]->DepositSources(sources, *owner, owner->elapsedTime, dt);
		#ifdef USE_OPENCL
		sources.SendToComputeBuffer();
		#endif
	}

	Electromagnetic::Update();

	// Advance the electric field

	yeeTool->AdvanceE(A,PMLx,PMLy,PMLz,sources);
	if (owner->conductor.size())
		yeeTool->UpdateInteriorBoundaryE(A,conductorMask);

	// Save the old magnetic field so final fields can be centered

	yeeTool->PrepCenteredFields(F,A);

	// Advance the magnetic field

	yeeTool->AdvanceB(A,PMLx,PMLy,PMLz);
	if (owner->conductor.size())
		yeeTool->UpdateInteriorBoundaryB(A,conductorMask);

	// Setup the final field for the particle push

	yeeTool->CenteredFields(F,A);
}



///////////////////////////////////
//   CURVILINEAR DIRECT SOLVER   //
///////////////////////////////////


CurvilinearDirectSolver::CurvilinearDirectSolver(Grid* theGrid):DirectSolver(theGrid)
{
	name = "Curvilinear-EM";
	typeCode = tw::module_type::curvilinearDirectSolver;
	A.Initialize(6,*this,owner);
}

void CurvilinearDirectSolver::Initialize()
{
	Field A0;

	Electromagnetic::Initialize();
	A0.Initialize(3,*this,owner);

	SetExteriorBoundaryConditionsE(A,Element(0),Element(1),Element(2));
	SetExteriorBoundaryConditionsE(A,Element(3),Element(4),Element(5));

	if (owner->gridGeometry==cylindrical)
	{
		A.SetBoundaryConditions(Element(0),xAxis,none,none);
		A.SetBoundaryConditions(Element(1),xAxis,dirichletWall,dirichletCell);
		A.SetBoundaryConditions(Element(2),xAxis,neumannWall,dirichletCell);
		A.SetBoundaryConditions(Element(3),xAxis,dirichletWall,dirichletCell);
		A.SetBoundaryConditions(Element(4,5),xAxis,none,none);
	}

	if (owner->restarted)
		return;

	// Initialize radiation fields

	LoadVectorPotential<0,1,2>(A0,-dth);
	#pragma omp parallel
	{
		for (VectorizingIterator<3> v(*this,true);v<v.end();++v)
		{
			#pragma omp simd
			for (tw::Int k=lb[3];k<=ub[3];k++)
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
		for (VectorizingIterator<3> v(*this,true);v<v.end();++v)
		{
			#pragma omp simd
			for (tw::Int k=lb[3];k<=ub[3];k++)
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
		for (VectorizingIterator<3> v(*this,true);v<v.end();++v)
		{
			#pragma omp simd
			for (tw::Int k=lb[3];k<=ub[3];k++)
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
	if (owner->gridGeometry==cylindrical && owner->X(0,1)<0.0)
	{
		for (k=lb[3];k<=ub[3];k++)
			for (j=lb[2];j<=ub[2];j++)
				A(1,j,k,0) = 0.0;
	}
	if (owner->gridGeometry==spherical && owner->X(0,1)<0.0)
	{
		for (k=lb[3];k<=ub[3];k++)
			for (j=lb[2];j<=ub[2];j++)
				A(1,j,k,0) = 0.0;
	}
	if (owner->gridGeometry==spherical && owner->X(0,2)<0.0)
	{
		for (k=lb[3];k<=ub[3];k++)
			for (i=lb[1];i<=ub[1];i++)
				A(i,1,k,1) = 0.0;
	}
	if (owner->gridGeometry==spherical && owner->X(ub[2],2)>pi)
	{
		for (k=lb[3];k<=ub[3];k++)
			for (i=lb[1];i<=ub[1];i++)
				A(i,ub[2],k,1) = 0.0;
	}
}

void CurvilinearDirectSolver::SetSingularPointsB()
{
	// Assume that F holds the old staggered B-field at this point
	// Then use curlE = -dB/dt to update B
	tw::Int i,j,k;
	if (owner->gridGeometry==cylindrical && owner->X(0,1)<0.0)
	{
		for (k=lb[3];k<=ub[3];k++)
			for (j=lb[2];j<=ub[2];j++)
			{
				A(1,j,k,4) = 0.0;
				A(1,j,k,5) = F(1,j,k,5) - 2.0*dt*A(1,j,k,1)/owner->X(1,1);
			}
	}
	if (owner->gridGeometry==spherical && owner->X(0,2)<0.0)
	{
		for (k=lb[3];k<=ub[3];k++)
			for (i=lb[1];i<=ub[1];i++)
			{
				A(i,1,k,5) = 0.0;
				A(i,1,k,3) = F(i,1,k,3) - 2.0*dt*A(i,1,k,2)/(owner->X(i,1)*sin(owner->X(1,2)));
			}
	}
	if (owner->gridGeometry==spherical && owner->X(ub[2],2)>pi)
	{
		for (k=lb[3];k<=ub[3];k++)
			for (i=lb[1];i<=ub[1];i++)
			{
				A(i,ub[2],k,5) = 0.0;
				A(i,ub[2],k,3) = F(i,ub[2],k,3) + 2.0*dt*A(i,dim[2],k,2)/(owner->X(i,1)*sin(owner->X(dim[2],2)));
			}
	}
}

void CurvilinearDirectSolver::Update()
{
	tw::Int i,j,k,r,kg;
	tw::vec3 S;

	// Add electric antenna currents

	if (owner->conductor.size())
	{
		#ifdef USE_OPENCL
		sources.ReceiveFromComputeBuffer();
		#endif
		for (r=0;r<owner->conductor.size();r++)
			if (owner->conductor[r]->electricCurrent)
				owner->conductor[r]->DepositSources(sources, *owner, owner->elapsedTime, dt);
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
	A.UpwardCopy(xAxis,Element(3,5),1);
	A.UpwardCopy(yAxis,Element(3,5),1);
	A.UpwardCopy(zAxis,Element(3,5),1);

	// Setup the final field for the particle push

	#pragma omp parallel
	{
		for (VectorizingIterator<3> v(*this,true);v<v.end();++v)
		{
			#pragma omp simd
			for (tw::Int k=lb[3];k<=ub[3];k++)
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
