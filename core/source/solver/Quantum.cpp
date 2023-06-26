#include "simulation.h"
#include "quantum.h"
using namespace tw::bc;


////////////////////////////////
//                            //
//   ATOMIC PHYSICS MODULE    //
//                            //
////////////////////////////////


AtomicPhysics::AtomicPhysics(const std::string& name,Simulation* sim):Module(name,sim)
{
	keepA2Term = true;
	dipoleApproximation = true;
	alpha = 0.0072973525664;
	timeRelaxingToGround = 0.0;

	H.form = qo::schroedinger;
	H.qorb = -1.0;
	H.morb = 1.0;
	H.qnuc = 1.0;
	H.rnuc = 0.01;
	H.B0 = 0.0;

	A4.Initialize(4,*this,owner,tw::grid::x);
	Ao4.Initialize(4,*this,owner,tw::grid::x);
	J4.Initialize(4,*this,owner,tw::grid::x);

	photonPropagator = NULL;

	#ifdef USE_OPENCL
	InitializeCLProgram("quantum.cl");
	A4.InitializeComputeBuffer();
	Ao4.InitializeComputeBuffer();
	J4.InitializeComputeBuffer();
	#endif

	directives.Add("orbiting charge",new tw::input::Float(&H.qorb),false);
	directives.Add("orbiting mass",new tw::input::Float(&H.morb),false);
	directives.Add("B0",new tw::input::Vec3(&H.B0),false);
	directives.Add("keep a2 term",new tw::input::Bool(&keepA2Term),false);
	directives.Add("dipole approximation",new tw::input::Bool(&dipoleApproximation),false);
	directives.Add("relaxation time",new tw::input::Float(&timeRelaxingToGround),false);
	directives.Add("soft core potential charge",new tw::input::Custom,false);
	directives.Add("bachelet potential",new tw::input::Custom,false);
}

AtomicPhysics::~AtomicPhysics()
{
	if (photonPropagator!=NULL)
		owner->RemoveTool(photonPropagator);
	// release any base class OpenCL kernels here
}

tw::Float AtomicPhysics::GetSphericalPotential(tw::Float r) const
{
	return H.qnuc/sqrt(sqr(H.rnuc) + sqr(r));
	//return (H.qnuc / r) * (H.c1*erf(sqrt(H.a1)*r) + H.c2*erf(sqrt(H.a2)*r));
}

void AtomicPhysics::Initialize()
{
	Module::Initialize();
	#ifdef USE_OPENCL
	photonPropagator->SetupComputeKernels(A4,Ao4,J4);
	#endif

	// Boundary conditions should preserve hermiticity
	// One way is to have A = 0 and grad(psi)=0 for components normal to boundary

	tw::bc::fld psiDefaultBC,A4DefaultBC;
	psiDefaultBC = fld::neumannWall;
	A4DefaultBC = fld::dirichletWall;
	psi_r.SetBoundaryConditions(tw::grid::x,psiDefaultBC,psiDefaultBC);
	psi_r.SetBoundaryConditions(tw::grid::y,psiDefaultBC,psiDefaultBC);
	psi_r.SetBoundaryConditions(tw::grid::z,psiDefaultBC,psiDefaultBC);
	psi_i.SetBoundaryConditions(tw::grid::x,psiDefaultBC,psiDefaultBC);
	psi_i.SetBoundaryConditions(tw::grid::y,psiDefaultBC,psiDefaultBC);
	psi_i.SetBoundaryConditions(tw::grid::z,psiDefaultBC,psiDefaultBC);
	A4.SetBoundaryConditions(tw::grid::x,A4DefaultBC,A4DefaultBC);
	A4.SetBoundaryConditions(tw::grid::y,A4DefaultBC,A4DefaultBC);
	A4.SetBoundaryConditions(tw::grid::z,A4DefaultBC,A4DefaultBC);
	J4.SetBoundaryConditions(tw::grid::x,A4DefaultBC,A4DefaultBC);
	J4.SetBoundaryConditions(tw::grid::y,A4DefaultBC,A4DefaultBC);
	J4.SetBoundaryConditions(tw::grid::z,A4DefaultBC,A4DefaultBC);

	switch (owner->gridGeometry)
	{
		case tw::grid::cartesian:
			break;
		case tw::grid::cylindrical:
			A4.SetBoundaryConditions(Element(0),tw::grid::x,fld::neumannWall,fld::dirichletWall);
			A4.SetBoundaryConditions(Element(3),tw::grid::x,fld::neumannWall,fld::dirichletWall);
			J4.SetBoundaryConditions(Element(0),tw::grid::x,fld::neumannWall,fld::dirichletWall);
			J4.SetBoundaryConditions(Element(3),tw::grid::x,fld::neumannWall,fld::dirichletWall);
			break;
		case tw::grid::spherical:
			A4.SetBoundaryConditions(Element(0),tw::grid::y,fld::neumannWall,fld::neumannWall);
			A4.SetBoundaryConditions(Element(1),tw::grid::y,fld::neumannWall,fld::neumannWall);
			J4.SetBoundaryConditions(Element(0),tw::grid::y,fld::neumannWall,fld::neumannWall);
			J4.SetBoundaryConditions(Element(1),tw::grid::y,fld::neumannWall,fld::neumannWall);
			break;
	}

	// Error check quantum numbers , print state information

	for (auto w : waveFunction)
	{
		*owner->tw_out << w->name << ": energy = " << w->Energy(H) << " , Cn = " << w->NormalizationConstant(H) << std::endl;
		if (!w->GoodQuantumNumbers(H))
			throw tw::FatalError("Bad quantum numbers detected in module <"+name+">");
	}
}

void AtomicPhysics::ExchangeResources()
{
	Module::ExchangeResources();
	PublishResource(&J4,"qo:j4");
	PublishResource(&H,"qo:H");
}

void AtomicPhysics::FormPotentials(tw::Float t)
{
	// residual charge has to be in desired units (it is Z for atomic units, Z*sqrt(alpha) for natural)
	// vector potential in wave block must also be in desired units, unlike earlier versions which used mc^2/e in every case

	#pragma omp parallel firstprivate(t)
	{
		tw::Float phiNow,r;
		tw::vec3 A0,A1,r_curv,r_cart;
		for (auto cell : EntireCellRange(*this))
		{
			r = owner->SphericalRadius(owner->Pos(cell));
			phiNow = GetSphericalPotential(r);
			r_cart = dipoleApproximation ? tw::vec3(0,0,0) : owner->Pos(cell);
			owner->CurvilinearToCartesian(&r_cart);
			A0 = A1 = tw::vec3(-0.5*r_cart.y*H.B0.z,0.5*r_cart.x*H.B0.z,0.0);
			for (tw::Int s=0;s<wave.size();s++)
			{
				A0 += wave[s]->VectorPotential(t-dt,r_cart);
				A1 += wave[s]->VectorPotential(t,r_cart);
			}
			r_curv = owner->Pos(cell);
			owner->TangentVectorToCurvilinear(&A0,r_curv);
			owner->TangentVectorToCurvilinear(&A1,r_curv);

			Ao4(cell,0) = phiNow;
			Ao4(cell,1) = A0.x;
			Ao4(cell,2) = A0.y;
			Ao4(cell,3) = A0.z;

			A4(cell,0) = phiNow;
			A4(cell,1) = A1.x;
			A4(cell,2) = A1.y;
			A4(cell,3) = A1.z;
		}
	}

	Ao4.CopyFromNeighbors();
	//Ao4.ApplyBoundaryCondition();
	A4.CopyFromNeighbors();
	//A4.ApplyBoundaryCondition();
}

void AtomicPhysics::FormGhostCellPotentials(tw::Float t)
{
	for (tw::Int ax=1;ax<=3;ax++)
		if (A4.Dim(ax)>1)
			#pragma omp parallel firstprivate(t,ax)
			{
				for (auto s : StripRange(*this,ax,strongbool::no))
					for (tw::Int ghostCell=0;ghostCell<=Dim(s.Axis())+1;ghostCell+=Dim(s.Axis())+1)
					{
						tw::vec3 pos(owner->Pos(s,ghostCell));
						tw::vec3 A3(-0.5*pos.y*H.B0.z,0.5*pos.x*H.B0.z,0.0);
						for (tw::Int wv=0;wv<wave.size();wv++)
							A3 += wave[wv]->VectorPotential(t,pos);
						if ((ghostCell==0 && owner->n0[ax]==MPI_PROC_NULL) || (ghostCell!=0 && owner->n1[ax]==MPI_PROC_NULL))
						{
							A4(s,ghostCell,1) = A3.x;
							A4(s,ghostCell,2) = A3.y;
							A4(s,ghostCell,3) = A3.z;
						}
					}
			}
}

tw::vec4 AtomicPhysics::GetA4AtOrigin()
{
	tw::Int s;
	tw::vec3 aNow,r_cart;
	tw::vec4 A;

	aNow = 0.0;
	r_cart = tw::vec3(0,0,0);
	A[0] = GetSphericalPotential(0.0);
	for (s=0;s<wave.size();s++)
		aNow += wave[s]->VectorPotential(owner->elapsedTime,r_cart);
	A[1] = aNow.x;
	A[2] = aNow.y;
	A[3] = aNow.z;
	return A;
}

void AtomicPhysics::VerifyInput()
{
	Module::VerifyInput();
	if (owner->gridGeometry!=tw::grid::cartesian)
		if (Norm(H.B0)!=0.0)
			throw tw::FatalError("Static B field assumes Cartesian geometry.");
	for (auto tool : moduleTool)
	{
		QState *state = dynamic_cast<QState*>(tool);
		if (state) waveFunction.push_back(state);
	}
	for (auto tool : moduleTool)
	{
		photonPropagator = dynamic_cast<LorentzPropagator*>(tool);
		if (photonPropagator!=NULL)
			break;
	}
	if (photonPropagator==NULL)
		photonPropagator = (LorentzPropagator*)owner->CreateTool("default_photons",tw::tool_type::lorentzPropagator);
}

void AtomicPhysics::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	tw::dnum q,r;
	std::string word;
	Module::ReadInputFileDirective(inputString,command);
	// note: examples of charge are geared toward atomic units
	// if using natural units, unit of charge is sqrt(alpha) ~ 0.085
	if (command=="soft core potential charge") // eg, soft core potential , charge = 1.0 , radius = 0.01
	{
		inputString >> word >> q >> word >> word >> r;
		H.qnuc = q >> native;
		H.rnuc = r >> native;
	}
	if (command=="bachelet potential") // eg, bachelet potential = 1.0 1.0 1.0 0.1 0.5
	{
		inputString >> word;
		inputString >> q;
		inputString >> H.c1 >> H.c2 >> H.a1 >> H.a2;
		H.qnuc = q >> native;
	}
}

void AtomicPhysics::ReadCheckpoint(std::ifstream& inFile)
{
	Module::ReadCheckpoint(inFile);
	psi_r.ReadCheckpoint(inFile);
	psi_i.ReadCheckpoint(inFile);
	J4.ReadCheckpoint(inFile);
	Ao4.ReadCheckpoint(inFile);
	A4.ReadCheckpoint(inFile);
	#ifdef USE_OPENCL
	psi_r.SendToComputeBuffer();
	psi_i.SendToComputeBuffer();
	J4.SendToComputeBuffer();
	Ao4.SendToComputeBuffer();
	A4.SendToComputeBuffer();
	#endif
}

void AtomicPhysics::WriteCheckpoint(std::ofstream& outFile)
{
	#ifdef USE_OPENCL
	psi_r.ReceiveFromComputeBuffer();
	psi_i.ReceiveFromComputeBuffer();
	J4.ReceiveFromComputeBuffer();
	Ao4.ReceiveFromComputeBuffer();
	A4.ReceiveFromComputeBuffer();
	#endif
	Module::WriteCheckpoint(outFile);
	psi_r.WriteCheckpoint(outFile);
	psi_i.WriteCheckpoint(outFile);
	J4.WriteCheckpoint(outFile);
	Ao4.WriteCheckpoint(outFile);
	A4.WriteCheckpoint(outFile);
}

void AtomicPhysics::Report(Diagnostic& diagnostic)
{
	tw::vec3 ENow,ANow;
	for (auto w : wave)
	{
		ANow += w->VectorPotential(owner->elapsedTime,tw::vec3(0,0,0));
		ENow -= dti*w->VectorPotential(owner->elapsedTime+0.5*dt,tw::vec3(0,0,0));
		ENow += dti*w->VectorPotential(owner->elapsedTime-0.5*dt,tw::vec3(0,0,0));
	}
	diagnostic.Float("Ex",ENow.x,true);
	diagnostic.Float("Ey",ENow.y,true);
	diagnostic.Float("Ez",ENow.z,true);
	diagnostic.Float("Ax",ANow.x,true);
	diagnostic.Float("Ay",ANow.y,true);
	diagnostic.Float("Az",ANow.z,true);

	diagnostic.Field("rho",J4,0,tw::dims::charge_density,"$\\rho$");
	diagnostic.Field("Jx",J4,1,tw::dims::current_density,"$j_x$");
	diagnostic.Field("Jy",J4,2,tw::dims::current_density,"$j_y$");
	diagnostic.Field("Jz",J4,3,tw::dims::current_density,"$j_z$");

	diagnostic.Field("phi",A4,0,tw::dims::scalar_potential,"$\\phi$");
	diagnostic.Field("Ax",A4,1,tw::dims::vector_potential,"$A_x$");
	diagnostic.Field("Ay",A4,2,tw::dims::vector_potential,"$A_y$");
	diagnostic.Field("Az",A4,3,tw::dims::vector_potential,"$A_z$");

	// ScalarField temp;
	// temp.Initialize(*this,owner);
	// for (tw::Int k=1;k<=dim[3];k++)
	// 	for (tw::Int j=1;j<=dim[2];j++)
	// 		for (tw::Int i=1;i<=dim[1];i++)
	// 			temp(i,j,k) = div<1,2,3>(J4,i,j,k,*owner);
	// diagnostic.Field("divJ",temp,0,tw::dims::none,"$\\nabla\\cdot {\\bf j}$");
}



///////////////////////////////////////
//                                   //
// NON-RELATIVISTIC PARTICLE (TDSE)  //
//                                   //
///////////////////////////////////////


Schroedinger::Schroedinger(const std::string& name,Simulation* sim):AtomicPhysics(name,sim)
{
	if (native.native!=tw::units::atomic)
		throw tw::FatalError("Schroedinger module requires <native units = atomic>");

	// Should move OpenCL stuff into the propagator tool
	H.form = qo::schroedinger;
	#ifndef USE_OPENCL
	propagator = NULL;
	#endif
	psi0.Initialize(*this,owner);
	psi1.Initialize(*this,owner);
	v.Initialize(*this,owner);
	w.Initialize(*this,owner);
	scratch.Initialize(*this,owner);

	#ifdef USE_OPENCL

	cl_int err;
	applyNumerator = clCreateKernel(program,"ApplyNumerator",&err);
	applyDenominator = clCreateKernel(program,"ApplyDenominator",&err);
	chargeKernel = clCreateKernel(program,"DepositCharge",&err);
	currentKernel = clCreateKernel(program,"DepositCurrent",&err);
	stitchKernel = clCreateKernel(program,"Stitch",&err);

	scratch.InitializeComputeBuffer();
	v.InitializeComputeBuffer();
	w.InitializeComputeBuffer();
	psi0.InitializeComputeBuffer();
	psi1.InitializeComputeBuffer();
	J4.InitializeComputeBuffer();

	clSetKernelArg(applyNumerator,0,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(applyNumerator,1,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(applyNumerator,2,sizeof(cl_mem),&owner->metricsBuffer);

	clSetKernelArg(applyDenominator,0,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(applyDenominator,1,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(applyDenominator,2,sizeof(cl_mem),&owner->metricsBuffer);
	clSetKernelArg(applyDenominator,6,sizeof(cl_mem),&scratch.computeBuffer);
	clSetKernelArg(applyDenominator,7,sizeof(cl_mem),&v.computeBuffer);
	clSetKernelArg(applyDenominator,8,sizeof(cl_mem),&w.computeBuffer);

	clSetKernelArg(chargeKernel,0,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(chargeKernel,1,sizeof(cl_mem),&J4.computeBuffer);

	clSetKernelArg(currentKernel,0,sizeof(cl_mem),&psi0.computeBuffer);
	clSetKernelArg(currentKernel,1,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(currentKernel,2,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(currentKernel,3,sizeof(cl_mem),&J4.computeBuffer);
	clSetKernelArg(currentKernel,4,sizeof(cl_mem),&owner->metricsBuffer);

	clSetKernelArg(stitchKernel,0,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(stitchKernel,2,sizeof(cl_mem),&v.computeBuffer);
	clSetKernelArg(stitchKernel,3,sizeof(cl_mem),&w.computeBuffer);

	#endif
}

Schroedinger::~Schroedinger()
{
	#ifndef USE_OPENCL
	if (propagator!=NULL)
		owner->RemoveTool(propagator);
	#endif
	#ifdef USE_OPENCL
	clReleaseKernel(applyNumerator);
	clReleaseKernel(applyDenominator);
	clReleaseKernel(chargeKernel);
	clReleaseKernel(currentKernel);
	clReleaseKernel(stitchKernel);
	#endif
}

void Schroedinger::ExchangeResources()
{
	AtomicPhysics::ExchangeResources();
	PublishResource(&psi1,"qo:psi");
}

void Schroedinger::VerifyInput()
{
	AtomicPhysics::VerifyInput();
	#ifndef USE_OPENCL
	if (propagator==NULL)
		propagator = (SchroedingerPropagator*)owner->CreateTool("TDSE",tw::tool_type::schroedingerPropagator);
	#endif
}

void Schroedinger::Initialize()
{
	AtomicPhysics::Initialize();

	const tw::bc::fld psiDefaultBC = fld::neumannWall;
	psi0.SetBoundaryConditions(tw::grid::x,psiDefaultBC,psiDefaultBC);
	psi0.SetBoundaryConditions(tw::grid::y,psiDefaultBC,psiDefaultBC);
	psi0.SetBoundaryConditions(tw::grid::z,psiDefaultBC,psiDefaultBC);
	psi1.SetBoundaryConditions(tw::grid::x,psiDefaultBC,psiDefaultBC);
	psi1.SetBoundaryConditions(tw::grid::y,psiDefaultBC,psiDefaultBC);
	psi1.SetBoundaryConditions(tw::grid::z,psiDefaultBC,psiDefaultBC);

	// Solve for the lowest energy s-state on a spherical grid.
	// This is used only to print the numerical ground state energy level.
	const tw::Float maxR = owner->SphericalRadius(GlobalCorner(*owner)+GlobalPhysicalSize(*owner));
	const tw::Float dr = dx(*owner) * owner->ScaleFactor(1,tw::vec3(tw::small_pos,0.0,0.0));
	const tw::Float r = maxR>30.0 ? 30.0 : maxR;
	const tw::Int dim = MyCeil(r/dr);
	std::valarray<tw::Float> eigenvector(dim),phi_r(dim);
	for (tw::Int i=0;i<dim;i++)
		phi_r[i] = GetSphericalPotential((tw::Float(i)+0.5)*dr);
	tw::Float groundStateEnergy = GetSphericalGroundState(eigenvector,phi_r,dr);
	(*owner->tw_out) << "Numerical ground state energy = " << groundStateEnergy << std::endl;

	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*this))
			for (auto w : waveFunction)
				psi1(cell) += w->Amplitude(H,owner->Pos(cell),0.0,0);
	}
	Normalize();
	psi0 = psi1;
	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*this))
			J4(cell,0) = norm(psi1(cell));
	}

	FormPotentials(owner->elapsedTime);

	#ifdef USE_OPENCL
	psi1.SendToComputeBuffer();
	A4.SendToComputeBuffer();
	#endif
}

void Schroedinger::ReadCheckpoint(std::ifstream& inFile)
{
	AtomicPhysics::ReadCheckpoint(inFile);
	psi0.ReadCheckpoint(inFile);
	psi1.ReadCheckpoint(inFile);
	#ifdef USE_OPENCL
	psi0.SendToComputeBuffer();
	psi1.SendToComputeBuffer();
	#endif
}

void Schroedinger::WriteCheckpoint(std::ofstream& outFile)
{
	AtomicPhysics::WriteCheckpoint(outFile);
	#ifdef USE_OPENCL
	psi0.ReceiveFromComputeBuffer();
	psi1.ReceiveFromComputeBuffer();
	#endif
	psi0.WriteCheckpoint(outFile);
	psi1.WriteCheckpoint(outFile);
}

#ifdef USE_OPENCL
void Schroedinger::Update()
{
	Field mpi_packet;
	tw::Float partitionFactor,relax,totalProbability;
	partitionFactor = 1.0/tw::Float(owner->Dimensionality());
	relax = owner->elapsedTime < timeRelaxingToGround ? 1.0 : 0.0;

	// Clear source vector
	J4.MADDComputeBuffer(0.0,0.0);

	// Update vector potential
	tw::vec4 A40 = GetA4AtOrigin();
	A4.FillComputeBufferVec4(Element(1,3),A40);

	// Set up compute arguments
	// argument 3 is the strip argument
	clSetKernelArg(applyNumerator,4,sizeof(tw::Float),&partitionFactor);
	clSetKernelArg(applyNumerator,5,sizeof(tw::Float),&relax);
	clSetKernelArg(applyDenominator,4,sizeof(tw::Float),&partitionFactor); // 3 is the strip
	clSetKernelArg(applyDenominator,5,sizeof(tw::Float),&relax);

	// Start Charge Deposition
	owner->CellUpdateProtocol(chargeKernel);

	// Loop over Dimensions Updating Wavefunction and Currents
	for (tw::Int i=1;i<=3;i++)
	{
		if (owner->localCells[i]>1)
		{
			CopyComputeBuffer(psi0,psi1);
			tw::Int j = i > 2 ? i-2 : i+1; // 2 , 3 , 1
			tw::Int k = i > 1 ? i-1 : i+2; // 3 , 1 , 2
			DiscreteSpace mpi_layout;
			tw::Int dim[4] = { 1, owner->localCells[j], owner->localCells[k], 1 };
			tw::Int gdim[4] = { 1, owner->globalCells[j], owner->globalCells[k], 1 };
			tw::Int dom[4] = { 1, owner->domains[j], owner->domains[k], 1 };
			mpi_layout.Resize(dim,gdim,dom,globalCorner,globalSize);
			mpi_packet.Initialize(16,mpi_layout,owner);
			mpi_packet.InitializeComputeBuffer();
			clSetKernelArg(applyDenominator,9,sizeof(cl_mem),&mpi_packet.computeBuffer);
			owner->StripUpdateProtocol(applyNumerator,i,3);
			owner->StripUpdateProtocol(applyDenominator,i,3);

			// Global Integration
			mpi_packet.ReceiveFromComputeBuffer();
			ComputeAlphasAndBetas<tw::Complex>(&owner->strip[i],owner->localCells2[j]*owner->localCells2[k],(tw::Complex*)&mpi_packet(0,0,0,0));
			mpi_packet.SendToComputeBuffer();
			clSetKernelArg(stitchKernel,1,sizeof(cl_mem),&owner->stripBuffer[i]);
			clSetKernelArg(stitchKernel,4,sizeof(cl_mem),&mpi_packet.computeBuffer);
			owner->LocalUpdateProtocol(stitchKernel);
			psi1.UpdateGhostCellsInComputeBuffer();

			if (owner->elapsedTime > timeRelaxingToGround)
				owner->StripUpdateProtocol(currentKernel,i,5);
		}
	}

	// Finish Charge Deposition
	owner->CellUpdateProtocol(chargeKernel);

	// Normalize After Relaxation Step
	if (owner->elapsedTime < timeRelaxingToGround)
	{
		CopyComputeBuffer(scratch,psi1);
		scratch.DestructiveComplexMod2ComputeBuffer();
		scratch.ZeroGhostCellsInComputeBuffer();
		scratch.WeightComputeBufferByVolume(*owner,0.0);
		totalProbability = scratch.DestructiveSumComputeBuffer();
		owner->strip[0].AllSum(&totalProbability,&totalProbability,sizeof(tw::Float),0);
		psi1.MADDComputeBuffer(1.0/sqrt(totalProbability),0.0);
	}

	// Get J4 for Bohmian Trajectories
	J4.UpdateGhostCellsInComputeBuffer();
	J4.ReceiveFromComputeBuffer();
}
#else
void Schroedinger::Update()
{
	tw::Complex dtc = owner->elapsedTime < timeRelaxingToGround ? -ii*dt : dt;
	FormPotentials(owner->elapsedTime);
	J4 = 0.0;

	propagator->DepositCurrent(tw::grid::t,psi0,psi1,A4,J4,dtc);

	psi0 = psi1;
	propagator->ApplyNumerator(tw::grid::x,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::grid::x,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::grid::x,psi0,psi1,A4,J4,dtc);

	psi0 = psi1;
	propagator->ApplyNumerator(tw::grid::y,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::grid::y,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::grid::y,psi0,psi1,A4,J4,dtc);

	psi0 = psi1;
	propagator->ApplyNumerator(tw::grid::z,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::grid::z,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::grid::z,psi0,psi1,A4,J4,dtc);

	propagator->DepositCurrent(tw::grid::t,psi0,psi1,A4,J4,dtc);

	J4.CopyFromNeighbors();
	J4.ApplyBoundaryCondition();
	if (owner->elapsedTime < timeRelaxingToGround)
		Normalize();
}
// void Schroedinger::Update()
// {
// 	// This advances the wavefunction conservatively, but uses the naive
// 	// centered differenced current deposition
// 	tw::Complex dtc = owner->elapsedTime < timeRelaxingToGround ? -ii*dt : dt;
// 	FormPotentials(owner->elapsedTime);
//
// 	psi0 = psi1;
// 	propagator->ApplyNumerator(tw::grid::x,psi1,A4,keepA2Term,dtc);
// 	propagator->ApplyDenominator(tw::grid::x,psi1,A4,keepA2Term,dtc);
// 	propagator->ApplyNumerator(tw::grid::y,psi1,A4,keepA2Term,dtc);
// 	propagator->ApplyDenominator(tw::grid::y,psi1,A4,keepA2Term,dtc);
// 	propagator->ApplyNumerator(tw::grid::z,psi1,A4,keepA2Term,dtc);
// 	propagator->ApplyDenominator(tw::grid::z,psi1,A4,keepA2Term,dtc);
// 	if (owner->elapsedTime < timeRelaxingToGround)
// 		Normalize();
// 	else
// 		UpdateJ4();
// }
#endif

void Schroedinger::UpdateJ4()
{
	// Here is the naive centered differencing to determine the non-relativistic probability current.
	// The results are usually pathological.  The conservative scheme should be used instead.
	tw::Complex psiNow,psi_x,psi_y,psi_z;
	if (owner->elapsedTime < timeRelaxingToGround)
		return;

	for (tw::Int k=1;k<=dim[3];k++)
		for (tw::Int j=1;j<=dim[2];j++)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				psiNow = psi1(i,j,k);
				psi_x = (psi1(i+1,j,k) - psi1(i-1,j,k))/owner->dL(i,j,k,1);
				psi_y = (psi1(i,j+1,k) - psi1(i,j-1,k))/owner->dL(i,j,k,2);
				psi_z = (psi1(i,j,k+1) - psi1(i,j,k-1))/owner->dL(i,j,k,3);

				J4(i,j,k,0) = norm(psiNow);
				J4(i,j,k,1) = -real((half*ii/H.morb) * (conj(psiNow)*psi_x - conj(psi_x)*psiNow)) - H.qorb*A4(i,j,k,1)*norm(psiNow)/H.morb;
				J4(i,j,k,2) = -real((half*ii/H.morb) * (conj(psiNow)*psi_y - conj(psi_y)*psiNow)) - H.qorb*A4(i,j,k,2)*norm(psiNow)/H.morb;
				J4(i,j,k,3) = -real((half*ii/H.morb) * (conj(psiNow)*psi_z - conj(psi_z)*psiNow)) - H.qorb*A4(i,j,k,3)*norm(psiNow)/H.morb;
			}
	J4.CopyFromNeighbors();
	J4.ApplyBoundaryCondition();
}

void Schroedinger::Normalize()
{
	tw::Float totalProbability = 0.0;
	for (auto cell : InteriorCellRange(*this))
		totalProbability += norm(psi1(cell)) * owner->dS(cell,0);
	owner->strip[0].AllSum(&totalProbability,&totalProbability,sizeof(tw::Float),0);
	psi1 *= 1.0/sqrt(totalProbability);
	psi1.CopyFromNeighbors();
	psi1.ApplyBoundaryCondition();
}

void Schroedinger::Report(Diagnostic& diagnostic)
{
	AtomicPhysics::Report(diagnostic);

	const tw::vec3 r0 = 0.0;
	ScalarField temp;
	temp.Initialize(*this,owner);

	for (auto cell : InteriorCellRange(*this))
		temp(cell) = norm(psi1(cell));
	diagnostic.VolumeIntegral("TotalProb",temp,0);
	diagnostic.FirstMoment("Dx",temp,0,r0,tw::grid::x);
	diagnostic.FirstMoment("Dy",temp,0,r0,tw::grid::y);
	diagnostic.FirstMoment("Dz",temp,0,r0,tw::grid::z);

	diagnostic.Field("psi_r",psi1,0,tw::dims::none,"$\\Re\\psi$");
	diagnostic.Field("psi_i",psi1,1,tw::dims::none,"$\\Im\\psi$");
}

void Schroedinger::StartDiagnostics()
{
	#ifdef USE_OPENCL
	psi1.ReceiveFromComputeBuffer();
	A4.ReceiveFromComputeBuffer();
	#endif
}


/////////////////////////////////////////
//                                     //
// Dual TDSE Coupled by Pauli Matrices //
//                                     //
/////////////////////////////////////////

Pauli::Pauli(const std::string& name,Simulation* sim):AtomicPhysics(name,sim)
{
	throw tw::FatalError("Pauli module is not supported in this version of TW.");
	if (native.native!=tw::units::plasma)
		throw tw::FatalError("Pauli module requires <native units = atomic>");
	H.form = qo::pauli;
	#ifndef USE_OPENCL
	propagator = NULL;
	#endif
	psi0.Initialize(*this,owner);
	psi1.Initialize(*this,owner);
	chi0.Initialize(*this,owner);
	chi1.Initialize(*this,owner);
	v.Initialize(*this,owner);
	w.Initialize(*this,owner);
	scratch.Initialize(*this,owner);
	#ifdef USE_OPENCL
	scratch.InitializeComputeBuffer();
	v.InitializeComputeBuffer();
	w.InitializeComputeBuffer();
	psi0.InitializeComputeBuffer();
	psi1.InitializeComputeBuffer();
	chi0.InitializeComputeBuffer();
	chi1.InitializeComputeBuffer();
	J4.InitializeComputeBuffer();
	#endif
}

Pauli::~Pauli()
{
	#ifndef USE_OPENCL
	if (propagator!=NULL)
		owner->RemoveTool(propagator);
	#endif
}

void Pauli::VerifyInput()
{
	AtomicPhysics::VerifyInput();
	#ifndef USE_OPENCL
	if (propagator==NULL)
		propagator = (SchroedingerPropagator*)owner->CreateTool("TDSE",tw::tool_type::schroedingerPropagator);
	#endif
}

void Pauli::Initialize()
{
	AtomicPhysics::Initialize();
}

void Pauli::ReadCheckpoint(std::ifstream& inFile)
{
	AtomicPhysics::ReadCheckpoint(inFile);
	psi0.ReadCheckpoint(inFile);
	psi1.ReadCheckpoint(inFile);
	chi0.ReadCheckpoint(inFile);
	chi1.ReadCheckpoint(inFile);
	#ifdef USE_OPENCL
	psi0.SendToComputeBuffer();
	psi1.SendToComputeBuffer();
	chi0.SendToComputeBuffer();
	chi1.SendToComputeBuffer();
	#endif
}

void Pauli::WriteCheckpoint(std::ofstream& outFile)
{
	AtomicPhysics::WriteCheckpoint(outFile);
	#ifdef USE_OPENCL
	psi0.ReceiveFromComputeBuffer();
	psi1.ReceiveFromComputeBuffer();
	chi0.ReceiveFromComputeBuffer();
	chi1.ReceiveFromComputeBuffer();
	#endif
	psi0.WriteCheckpoint(outFile);
	psi1.WriteCheckpoint(outFile);
	chi0.WriteCheckpoint(outFile);
	chi1.WriteCheckpoint(outFile);
}

#ifdef USE_OPENCL
void Pauli::Update()
{
	cl_int err;
	cl_kernel spinKernel = clCreateKernel(program,"UpdateSpin",&err);

	// Spin
	clSetKernelArg(spinKernel,0,sizeof(cl_mem),&psi1.computeBuffer);
	clSetKernelArg(spinKernel,1,sizeof(cl_mem),&chi1.computeBuffer);
	clSetKernelArg(spinKernel,2,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(spinKernel,3,sizeof(cl_mem),&owner->metricsBuffer);
	owner->LocalUpdateProtocol(spinKernel);

	// ADD: something to handle normalization after a relaxation step

	clReleaseKernel(spinKernel);
}
#else
void Pauli::Update()
{
	// KNOWN PROBLEM: we have to add spin term to current density

	tw::Complex dtc = owner->elapsedTime < timeRelaxingToGround ? -ii*dt : dt;
	FormPotentials(owner->elapsedTime);
	J4 = 0.0;

	propagator->DepositCurrent(tw::grid::t,psi0,psi1,A4,J4,dtc);
	propagator->DepositCurrent(tw::grid::t,chi0,chi1,A4,J4,dtc);

	psi0 = psi1;
	chi0 = chi1;
	propagator->ApplyNumerator(tw::grid::x,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::grid::x,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::grid::x,psi0,psi1,A4,J4,dtc);
	propagator->ApplyNumerator(tw::grid::x,chi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::grid::x,chi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::grid::x,chi0,chi1,A4,J4,dtc);

	psi0 = psi1;
	chi0 = chi1;
	propagator->ApplyNumerator(tw::grid::y,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::grid::y,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::grid::y,psi0,psi1,A4,J4,dtc);
	propagator->ApplyNumerator(tw::grid::y,chi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::grid::y,chi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::grid::y,chi0,chi1,A4,J4,dtc);

	psi0 = psi1;
	chi0 = chi1;
	propagator->ApplyNumerator(tw::grid::z,psi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::grid::z,psi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::grid::z,psi0,psi1,A4,J4,dtc);
	propagator->ApplyNumerator(tw::grid::z,chi1,A4,keepA2Term,dtc);
	propagator->ApplyDenominator(tw::grid::z,chi1,A4,keepA2Term,dtc);
	propagator->DepositCurrent(tw::grid::z,chi0,chi1,A4,J4,dtc);

	propagator->DepositCurrent(tw::grid::t,psi0,psi1,A4,J4,dtc);
	propagator->DepositCurrent(tw::grid::t,chi0,chi1,A4,J4,dtc);

   	propagator->UpdateSpin(psi1,chi1,A4,alpha*dt);

	if (owner->elapsedTime < timeRelaxingToGround)
		Normalize();
}
#endif

void Pauli::Normalize()
{
	tw::Float totalProbability = 0.0;
	for (auto cell : InteriorCellRange(*this))
		totalProbability += (norm(psi1(cell))+norm(chi1(cell))) * owner->dS(cell,0);
	owner->strip[0].AllSum(&totalProbability,&totalProbability,sizeof(tw::Float),0);
	psi1 *= 1.0/sqrt(totalProbability);
	chi1 *= 1.0/sqrt(totalProbability);
	psi1.CopyFromNeighbors();
	psi1.ApplyBoundaryCondition();
	chi1.CopyFromNeighbors();
	chi1.ApplyBoundaryCondition();
}

void Pauli::Report(Diagnostic& diagnostic)
{
	AtomicPhysics::Report(diagnostic);

	const tw::vec3 r0 = 0.0;
	ScalarField temp;
	temp.Initialize(*this,owner);

	for (auto cell : InteriorCellRange(*this))
		temp(cell) = norm(psi1(cell)) + norm(chi1(cell));
	diagnostic.VolumeIntegral("TotalProb",temp,0);
	diagnostic.FirstMoment("Dx",temp,0,r0,tw::grid::x);
	diagnostic.FirstMoment("Dy",temp,0,r0,tw::grid::y);
	diagnostic.FirstMoment("Dz",temp,0,r0,tw::grid::z);

	diagnostic.Field("psi_r",psi1,0,tw::dims::none,"$\\Re\\psi$");
	diagnostic.Field("psi_i",psi1,1,tw::dims::none,"$\\Im\\psi$");
	diagnostic.Field("chi_r",chi1,0,tw::dims::none,"$\\Re\\chi$");
	diagnostic.Field("chi_i",chi1,1,tw::dims::none,"$\\Im\\chi$");

	for (auto cell : InteriorCellRange(*this))
		temp(cell) = norm(psi1(cell)) - norm(chi1(cell));
	diagnostic.Field("Sz",temp,0);
	diagnostic.VolumeIntegral("Sz",temp,0);
}

void Pauli::StartDiagnostics()
{
	#ifdef USE_OPENCL
	psi0.ReceiveFromComputeBuffer();
	psi1.ReceiveFromComputeBuffer();
	chi0.ReceiveFromComputeBuffer();
	chi1.ReceiveFromComputeBuffer();
	A4.ReceiveFromComputeBuffer();
	J4.ReceiveFromComputeBuffer();
	#endif
}


///////////////////////////////////////
//                                   //
//  RELATIVISTIC PARTICLE (KG EQN)   //
//                                   //
///////////////////////////////////////


KleinGordon::KleinGordon(const std::string& name,Simulation* sim) : AtomicPhysics(name,sim)
{
	// Wavefunction is in Hamiltonian 2-component representation.
	// This is the preliminary Hamiltonian form from Feshbach-Villars, Eq. 2.12.
	// Unlike the symmetric form, this allows a leap-frog scheme to be employed.
	// In this representation, component 0 is the usual scalar wavefunction.
	// Component 1 is psi1 = (id/dt - q*phi)*psi0/m
	// The symmetric form is obtained from (psi0+psi1)/sqrt(2) , (psi0-psi1)/sqrt(2)

	if (native.native!=tw::units::natural)
		throw tw::FatalError("KleinGordon module requires <native units = natural>");

	H.form = qo::klein_gordon;
	H.morb = 1.0;
	H.qorb = -sqrt(alpha);
	H.qnuc = sqrt(alpha);
	dipoleApproximation = false;
	psi_r.Initialize(2,*this,owner,tw::grid::x);
	psi_i.Initialize(2,*this,owner,tw::grid::x);

	#ifdef USE_OPENCL
	cl_int err;
	updatePsi = clCreateKernel(program,"KGAdvance1",&err);
	updateChi = clCreateKernel(program,"KGAdvance2",&err);

	psi_r.InitializeComputeBuffer();
	psi_i.InitializeComputeBuffer();

	clSetKernelArg(updatePsi,0,sizeof(cl_mem),&psi_r.computeBuffer);
	clSetKernelArg(updatePsi,1,sizeof(cl_mem),&psi_i.computeBuffer);
	clSetKernelArg(updatePsi,2,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(updatePsi,3,sizeof(cl_mem),&owner->metricsBuffer);

	clSetKernelArg(updateChi,0,sizeof(cl_mem),&psi_r.computeBuffer);
	clSetKernelArg(updateChi,1,sizeof(cl_mem),&psi_i.computeBuffer);
	clSetKernelArg(updateChi,2,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(updateChi,3,sizeof(cl_mem),&owner->metricsBuffer);
	#endif
}

KleinGordon::~KleinGordon()
{
	#ifdef USE_OPENCL
	clReleaseKernel(updatePsi);
	clReleaseKernel(updateChi);
	#endif
}

void KleinGordon::Initialize()
{
	AtomicPhysics::Initialize();

	for (auto cell : InteriorCellRange(*this))
	{
		const tw::vec3 pos = owner->Pos(cell);
		for (auto w : waveFunction)
		{
			psi_r(cell,0) += real(w->Amplitude(H,pos,0.0,0));
			psi_i(cell,0) += imag(w->Amplitude(H,pos,0.0,0));
			psi_r(cell,1) += real(w->Amplitude(H,pos,dth,1));
			psi_i(cell,1) += imag(w->Amplitude(H,pos,dth,1));
		}
	}

	psi_r.CopyFromNeighbors();
	psi_i.CopyFromNeighbors();
	FormPotentials(owner->elapsedTime);
	Normalize();
	UpdateJ4();

	#ifdef USE_OPENCL

	psi_r.SendToComputeBuffer();
	psi_i.SendToComputeBuffer();
	Ao4.SendToComputeBuffer();
	A4.SendToComputeBuffer();
	clSetKernelArg(updatePsi,4,sizeof(H.morb),&H.morb);
	clSetKernelArg(updatePsi,5,sizeof(H.qorb),&H.qorb);
	clSetKernelArg(updateChi,4,sizeof(H.morb),&H.morb);
	clSetKernelArg(updateChi,5,sizeof(H.qorb),&H.qorb);

	#endif
}

tw::Float KleinGordon::ComputeRho(const tw::cell& cell)
{
	return H.qorb*(norm(FV(cell,1.0)) - norm(FV(cell,-1.0)));
}

void KleinGordon::UpdateJ4()
{
	#pragma omp parallel
	{
		for (auto v : VectorStripRange<1>(*this,false))
		{
			for (tw::Int i=1;i<=dim[1];i++)
			{
				J4(v,i,0) = H.qorb*(norm(FV(v,i,1.0)) - norm(FV(v,i,-1.0)));
				J4(v,i,1) = 0.0;
				J4(v,i,2) = 0.0;
				J4(v,i,3) = 0.0;
			}
		}
	}
	J4.CopyFromNeighbors();
	J4.ApplyBoundaryCondition();
}

void KleinGordon::Normalize()
{
	tw::Float totalCharge = 0.0;
	for (auto cell : InteriorCellRange(*this))
		totalCharge += ComputeRho(cell) * owner->dS(cell,0);
	owner->strip[0].AllSum(&totalCharge,&totalCharge,sizeof(tw::Float),0);
	psi_r *= sqrt(fabs(H.qorb/totalCharge));
	psi_i *= sqrt(fabs(H.qorb/totalCharge));
	psi_r.CopyFromNeighbors();
	psi_r.ApplyBoundaryCondition();
	psi_i.CopyFromNeighbors();
	psi_i.ApplyBoundaryCondition();
}

#ifdef USE_OPENCL
void KleinGordon::Update()
{
	// Update wavefunction
	owner->LocalUpdateProtocol(updatePsi);
	// Boundary conditions and communications
	psi_r.UpdateGhostCellsInComputeBuffer(Element(0));
	psi_i.UpdateGhostCellsInComputeBuffer(Element(0));

	photonPropagator->Advance(A4,Ao4,J4,0.0,dt);
	FormGhostCellPotentials(owner->elapsedTime+dt);
	A4.SendGhostCellsToComputeBuffer();
	photonPropagator->MidstepEstimate(A4,Ao4);

	// Update auxiliary wavefunction
	owner->LocalUpdateProtocol(updateChi);
	// Boundary conditions and communications
	psi_r.UpdateGhostCellsInComputeBuffer(Element(1));
	psi_i.UpdateGhostCellsInComputeBuffer(Element(1));

	photonPropagator->UndoMidstepEstimate(A4,Ao4);
}
#else
void KleinGordon::Update()
{
	// Solve in the Hamiltonian leap-frog representation
	// dpsi/dt = -i*q*phi*psi - i*m*chi
	// dchi/dt = -i*q*phi*chi - i*m*psi + i*(Dk)^2(psi)/m
	// Leapfrog Scheme:
	// At start we know psi(n-1/2),chi(n)
	// Advance psi to n+1/2 using chi(n),A(n)
	// Advance chi to n+1 using psi(n+1/2),A(n+1/2)

	static const tw::Float q0 = H.qorb;
	static const tw::Float m0 = H.morb;
	static const tw::Int AB = tw::vec_align_bytes;
	#pragma omp parallel
	{
		alignas(AB) tw::Float Ur[dim[1]],Ui[dim[1]],Dr[dim[1]],Di[dim[1]];
		// Update psi
		for (auto v : VectorStripRange<1>(*this,false))
		{
			#pragma omp simd aligned(Ur,Ui:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				// unitary operator of time translation for diagonal part of Hamiltonian
				const tw::Float dq = dt*q0*A4(v,i,0);
				Ur[i-1] = cos(dq);
				Ui[i-1] = -sin(dq);
			}
			#pragma omp simd aligned(Dr,Di:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				Dr[i-1] = m0*psi_i(v,i,1);
				Di[i-1] = -m0*psi_r(v,i,1);
			}
			#pragma omp simd aligned(Ur,Ui,Dr,Di:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				psi_r(v,i,0) += dth*Dr[i-1];
				psi_i(v,i,0) += dth*Di[i-1];
				complex_multiply_assign(psi_r(v,i,0),psi_i(v,i,0),Ur[i-1],Ui[i-1]);
				psi_r(v,i,0) += dth*Dr[i-1];
				psi_i(v,i,0) += dth*Di[i-1];
			}
		}
	}

	psi_r.CopyFromNeighbors(Element(0));
	psi_r.ApplyBoundaryCondition(Element(0));
	psi_i.CopyFromNeighbors(Element(0));
	psi_i.ApplyBoundaryCondition(Element(0));

	photonPropagator->Advance(A4,Ao4,J4,0.0,dt);
	FormGhostCellPotentials(owner->elapsedTime+dt);
	photonPropagator->MidstepEstimate(A4,Ao4);

	#pragma omp parallel
	{
		alignas(AB) tw::Float Ur[dim[1]],Ui[dim[1]],Dr[dim[1]],Di[dim[1]];
		// Update chi
		for (auto v : VectorStripRange<1>(*this,false))
		{
			#pragma omp simd aligned(Ur,Ui:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				// unitary operator of time translation for diagonal part of Hamiltonian
				const tw::Float dq = dt*q0*A4(v,i,0);
				Ur[i-1] = cos(dq);
				Ui[i-1] = -sin(dq);
			}
			#pragma omp simd aligned(Dr,Di:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				// Evaluate i(Dk)^2 assuming div(A)=0
				// This adds i*del^2(psi)
				Dr[i-1] = -(psi_i.d2(v,i,0,1) + psi_i.d2(v,i,0,2) + psi_i.d2(v,i,0,3));
				Di[i-1] = (psi_r.d2(v,i,0,1) + psi_r.d2(v,i,0,2) + psi_r.d2(v,i,0,3));
				// This adds 2q*div(A*psi)
				Dr[i-1] += 2*q0*freq.x*(psi_r.sfwd(v,i,0,1)*A4.sfwd(v,i,1,1) - psi_r.sbak(v,i,0,1)*A4.sbak(v,i,1,1));
				Dr[i-1] += 2*q0*freq.y*(psi_r.sfwd(v,i,0,2)*A4.sfwd(v,i,2,2) - psi_r.sbak(v,i,0,2)*A4.sbak(v,i,2,2));
				Dr[i-1] += 2*q0*freq.z*(psi_r.sfwd(v,i,0,3)*A4.sfwd(v,i,3,3) - psi_r.sbak(v,i,0,3)*A4.sbak(v,i,3,3));
				Di[i-1] += 2*q0*freq.x*(psi_i.sfwd(v,i,0,1)*A4.sfwd(v,i,1,1) - psi_i.sbak(v,i,0,1)*A4.sbak(v,i,1,1));
				Di[i-1] += 2*q0*freq.y*(psi_i.sfwd(v,i,0,2)*A4.sfwd(v,i,2,2) - psi_i.sbak(v,i,0,2)*A4.sbak(v,i,2,2));
				Di[i-1] += 2*q0*freq.z*(psi_i.sfwd(v,i,0,3)*A4.sfwd(v,i,3,3) - psi_i.sbak(v,i,0,3)*A4.sbak(v,i,3,3));
				// This adds -i*q^2*A^2*psi
				Dr[i-1] += q0*q0*(sqr(A4(v,i,1))+sqr(A4(v,i,2))+sqr(A4(v,i,3)))*psi_i(v,i,0);
				Di[i-1] -= q0*q0*(sqr(A4(v,i,1))+sqr(A4(v,i,2))+sqr(A4(v,i,3)))*psi_r(v,i,0);
				// Finish by incorporating mass
				Dr[i-1] = Dr[i-1]/m0 + m0*psi_i(v,i,0);
				Di[i-1] = Di[i-1]/m0 - m0*psi_r(v,i,0);
			}
			#pragma omp simd aligned(Ur,Ui,Dr,Di:AB)
			for (tw::Int i=1;i<=dim[1];i++)
			{
				// Perform the integration
				psi_r(v,i,1) += dth*Dr[i-1];
				psi_i(v,i,1) += dth*Di[i-1];
				complex_multiply_assign(psi_r(v,i,1),psi_i(v,i,1),Ur[i-1],Ui[i-1]);
				psi_r(v,i,1) += dth*Dr[i-1];
				psi_i(v,i,1) += dth*Di[i-1];
			}
		}
	}

	psi_r.CopyFromNeighbors(Element(1));
	psi_r.ApplyBoundaryCondition(Element(1));
	psi_i.CopyFromNeighbors(Element(1));
	psi_i.ApplyBoundaryCondition(Element(1));

	photonPropagator->UndoMidstepEstimate(A4,Ao4);
}
#endif

void KleinGordon::Report(Diagnostic& diagnostic)
{
	AtomicPhysics::Report(diagnostic);

	const tw::vec3 r0 = 0.0;
	diagnostic.VolumeIntegral("TotalCharge",J4,0);
	diagnostic.FirstMoment("Dx",J4,0,r0,tw::grid::x);
	diagnostic.FirstMoment("Dy",J4,0,r0,tw::grid::y);
	diagnostic.FirstMoment("Dz",J4,0,r0,tw::grid::z);

	diagnostic.Field("psi0_r",psi_r,0,tw::dims::none,"$\\Re\\psi_0$");
	diagnostic.Field("psi1_r",psi_r,1,tw::dims::none,"$\\Re\\psi_1$");
}

void KleinGordon::StartDiagnostics()
{
	#ifdef USE_OPENCL
	psi_r.ReceiveFromComputeBuffer();
	psi_i.ReceiveFromComputeBuffer();
	A4.ReceiveFromComputeBuffer();
	#endif
	UpdateJ4();
}


///////////////////////////////////////
//                                   //
//   RELATIVISTIC PARTICLE (DIRAC)   //
//                                   //
///////////////////////////////////////


Dirac::Dirac(const std::string& name,Simulation* sim) : AtomicPhysics(name,sim)
{
	// Wavefunction is in standard representation
	// i.e., gamma^0 = diag(1,1,-1,-1)

	if (native.native!=tw::units::natural)
		throw tw::FatalError("Dirac module requires <native units = natural>");

	H.form = qo::dirac;
	H.morb = 1.0;
	H.qorb = -sqrt(alpha);
	H.qnuc = sqrt(alpha);
	dipoleApproximation = false;
	psi_r.Initialize(4,*this,owner,tw::grid::x);
	psi_i.Initialize(4,*this,owner,tw::grid::x);

	#ifdef USE_OPENCL
	cl_int err;
	leapFrog = clCreateKernel(program,"DiracAdvance",&err);

	psi_r.InitializeComputeBuffer();
	psi_i.InitializeComputeBuffer();

	#endif
}

Dirac::~Dirac()
{
	#ifdef USE_OPENCL
	clReleaseKernel(leapFrog);
	#endif
}

void Dirac::Initialize()
{
	AtomicPhysics::Initialize();
	psi_r.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::none);
	psi_i.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::none);

	#pragma omp parallel
	{
		for (auto cell : InteriorCellRange(*this))
		{
			tw::vec3 pos = owner->Pos(cell);
			for (auto w : waveFunction)
			{
				psi_r(cell,0) += real(w->Amplitude(H,pos,0.0,0));
				psi_r(cell,1) += real(w->Amplitude(H,pos,0.0,1));
				psi_r(cell,2) += real(w->Amplitude(H,pos,dth,2));
				psi_r(cell,3) += real(w->Amplitude(H,pos,dth,3));
				psi_i(cell,0) += imag(w->Amplitude(H,pos,0.0,0));
				psi_i(cell,1) += imag(w->Amplitude(H,pos,0.0,1));
				psi_i(cell,2) += imag(w->Amplitude(H,pos,dth,2));
				psi_i(cell,3) += imag(w->Amplitude(H,pos,dth,3));
			}
		}
	}

	psi_r.CopyFromNeighbors();
	psi_i.CopyFromNeighbors();
	FormPotentials(owner->elapsedTime);
	Normalize();
	UpdateJ4();

	#ifdef USE_OPENCL

	psi_r.SendToComputeBuffer();
	psi_i.SendToComputeBuffer();
	Ao4.SendToComputeBuffer();
	A4.SendToComputeBuffer();
	clSetKernelArg(leapFrog,0,sizeof(cl_mem),&psi_r.computeBuffer);
	clSetKernelArg(leapFrog,1,sizeof(cl_mem),&psi_i.computeBuffer);
	clSetKernelArg(leapFrog,2,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(leapFrog,3,sizeof(cl_mem),&owner->metricsBuffer);
	clSetKernelArg(leapFrog,4,sizeof(H.morb),&H.morb);
	clSetKernelArg(leapFrog,5,sizeof(H.qorb),&H.qorb);

	#endif
}

tw::Float Dirac::ComputeRho(const tw::cell& cell)
{
	tw::Float ans = 0.0;
	for (tw::Int c=0;c<4;c++)
		ans += sqr(psi_r(cell,c)) + sqr(psi_i(cell,c));
	return H.qorb*ans;
}

void Dirac::UpdateJ4()
{
	#pragma omp parallel
	{
		tw::Complex z0,z1,z2,z3;
		for (auto cell : InteriorCellRange(*this))
		{
			z0 = tw::Complex(psi_r(cell,0),psi_i(cell,0));
			z1 = tw::Complex(psi_r(cell,1),psi_i(cell,1));
			z2 = tw::Complex(psi_r(cell,2),psi_i(cell,2));
			z3 = tw::Complex(psi_r(cell,3),psi_i(cell,3));
			J4(cell,0) = H.qorb*(norm(z0)+norm(z1)+norm(z2)+norm(z3));
			J4(cell,1) = two*H.qorb*real(conj(z0)*z3 + z1*conj(z2));
			J4(cell,2) = two*H.qorb*imag(conj(z0)*z3 + z1*conj(z2));
			J4(cell,3) = two*H.qorb*real(z0*conj(z2) - z1*conj(z3));
		}
	}
	J4.CopyFromNeighbors();
	J4.ApplyBoundaryCondition();
}

void Dirac::Normalize()
{
	tw::Float totalCharge = 0.0;
	for (auto cell : InteriorCellRange(*this))
		totalCharge += ComputeRho(cell) * owner->dS(cell,0);
	owner->strip[0].AllSum(&totalCharge,&totalCharge,sizeof(tw::Float),0);
	psi_r *= sqrt(fabs(H.qorb/totalCharge));
	psi_i *= sqrt(fabs(H.qorb/totalCharge));
	psi_r.CopyFromNeighbors();
	psi_r.ApplyBoundaryCondition();
	psi_i.CopyFromNeighbors();
	psi_i.ApplyBoundaryCondition();
}

#ifdef USE_OPENCL
void Dirac::Update()
{
	// Solve in the standard representation
	// [gamma^mu(p_mu - qA_mu) - m]psi = 0
	// gamma^0 = | 1  0 |     gamma^i = |  0  t |
	//           | 0 -1 |               | -t  0 |
	// where t is the i^th Pauli matrix
	// Leapfrog Scheme:
	// At start we know psi0(n-1/2),psi1(n-1/2),psi2(n),psi3(n),A(n-1),A(n)
	// Advance psi0,psi1 to n+1/2 using psi2(n),psi3(n),A(n)
	// Advance psi2,psi3 to n+1 using psi0(n+1/2),psi1(n+1/2),A(n+1/2)

	cl_int components[4] = { 0,1,2,3 };
	tw::Float mass[2] = { H.morb , -H.morb };

	// Update upper pair in the bispinor
	clSetKernelArg(leapFrog,4,sizeof(tw::Float),&mass[0]);
	clSetKernelArg(leapFrog,6,sizeof(cl_int),&components[0]);
	clSetKernelArg(leapFrog,7,sizeof(cl_int),&components[1]);
	clSetKernelArg(leapFrog,8,sizeof(cl_int),&components[2]);
	clSetKernelArg(leapFrog,9,sizeof(cl_int),&components[3]);
	owner->LocalUpdateProtocol(leapFrog);
	// Boundary conditions and communications
	psi_r.UpdateGhostCellsInComputeBuffer(Element(0,1));
	psi_i.UpdateGhostCellsInComputeBuffer(Element(0,1));

	photonPropagator->Advance(A4,Ao4,J4,0.0,dt);
	FormGhostCellPotentials(owner->elapsedTime+dt);
	A4.SendGhostCellsToComputeBuffer();
	photonPropagator->MidstepEstimate(A4,Ao4);

	// Update lower pair in the bispinor
	clSetKernelArg(leapFrog,4,sizeof(tw::Float),&mass[1]);
	clSetKernelArg(leapFrog,6,sizeof(cl_int),&components[2]);
	clSetKernelArg(leapFrog,7,sizeof(cl_int),&components[3]);
	clSetKernelArg(leapFrog,8,sizeof(cl_int),&components[0]);
	clSetKernelArg(leapFrog,9,sizeof(cl_int),&components[1]);
	owner->LocalUpdateProtocol(leapFrog);
	// Boundary conditions and communications
	psi_r.UpdateGhostCellsInComputeBuffer(Element(2,3));
	psi_i.UpdateGhostCellsInComputeBuffer(Element(2,3));

	photonPropagator->UndoMidstepEstimate(A4,Ao4);
}
#else
void Dirac::Update()
{
	// Solve in the standard representation
	// [gamma^mu(p_mu - qA_mu) - m]psi = 0
	// gamma^0 = | 1  0 |     gamma^i = |  0  t |
	//           | 0 -1 |               | -t  0 |
	// where t is the i^th Pauli matrix
	// Leapfrog Scheme:
	// At start we know psi0(n-1/2),psi1(n-1/2),psi2(n),psi3(n),A(n-1),A(n)
	// Advance psi0,psi1 to n+1/2 using psi2(n),psi3(n),A(n)
	// Advance psi2,psi3 to n+1 using psi0(n+1/2),psi1(n+1/2),A(n+1/2)

	LeapFrog<0,1,2,3>(1.0);
	psi_r.CopyFromNeighbors(Element(0,1));
	psi_r.ApplyBoundaryCondition(Element(0,1));
	psi_i.CopyFromNeighbors(Element(0,1));
	psi_i.ApplyBoundaryCondition(Element(0,1));

	photonPropagator->Advance(A4,Ao4,J4,0.0,dt);
	FormGhostCellPotentials(owner->elapsedTime+dt);
	photonPropagator->MidstepEstimate(A4,Ao4);

	LeapFrog<2,3,0,1>(-1.0);
	psi_r.CopyFromNeighbors(Element(2,3));
	psi_r.ApplyBoundaryCondition(Element(2,3));
	psi_i.CopyFromNeighbors(Element(2,3));
	psi_i.ApplyBoundaryCondition(Element(2,3));

	photonPropagator->UndoMidstepEstimate(A4,Ao4);
}
#endif

void Dirac::Report(Diagnostic& diagnostic)
{
	AtomicPhysics::Report(diagnostic);

	const tw::vec3 r0 = 0.0;
	diagnostic.VolumeIntegral("TotalCharge",J4,0);
	diagnostic.FirstMoment("Dx",J4,0,r0,tw::grid::x);
	diagnostic.FirstMoment("Dy",J4,0,r0,tw::grid::y);
	diagnostic.FirstMoment("Dz",J4,0,r0,tw::grid::z);

	diagnostic.Field("psi0_r",psi_r,0,tw::dims::none,"$\\Re\\psi_0$");
	diagnostic.Field("psi1_r",psi_r,1,tw::dims::none,"$\\Re\\psi_1$");
	diagnostic.Field("psi2_r",psi_r,2,tw::dims::none,"$\\Re\\psi_2$");
	diagnostic.Field("psi3_r",psi_r,3,tw::dims::none,"$\\Re\\psi_3$");
	diagnostic.Field("psi0_i",psi_i,0,tw::dims::none,"$\\Im\\psi_0$");
	diagnostic.Field("psi1_i",psi_i,1,tw::dims::none,"$\\Im\\psi_1$");
	diagnostic.Field("psi2_i",psi_i,2,tw::dims::none,"$\\Im\\psi_2$");
	diagnostic.Field("psi3_i",psi_i,3,tw::dims::none,"$\\Im\\psi_3$");
}

void Dirac::StartDiagnostics()
{
	#ifdef USE_OPENCL
	psi_r.ReceiveFromComputeBuffer();
	psi_i.ReceiveFromComputeBuffer();
	A4.ReceiveFromComputeBuffer();
	#endif
	UpdateJ4();
}

PopulationDiagnostic::PopulationDiagnostic(const std::string& name,Simulation *sim) : Module(name,sim)
{
	updateSequencePriority = tw::priority::diagnostic;
	H = NULL;
	psi = NULL;
}

bool PopulationDiagnostic::InspectResource(void* resource,const std::string& description)
{
	// The Hamiltonian parameters and wavefunction we want to analyze come from some other module.
	// We get this data via the publisher-consumer resource exchange mechanism.
	if (description=="qo:H")
	{
		H = (HamiltonianParameters*)resource;
		return true;
	}
	if (description=="qo:psi")
	{
		psi = (ComplexField*)resource;
		return true;
	}
	return false;
}

void PopulationDiagnostic::VerifyInput()
{
	// The reference states are encapsulated as QState tools.
	// These are attached explicitly by the user in the input file.
	for (auto tool : moduleTool)
	{
		QState *state = dynamic_cast<QState*>(tool);
		if (state) refState.push_back(state);
	}
	// EM waves are also attached in the input file.
	// Module base class already populates the wave list.
	// However the user has to know to attach waves to all modules that need them.
}

void PopulationDiagnostic::Initialize()
{
	for (auto ref : refState)
		if (!ref->GoodQuantumNumbers(*H))
			throw tw::FatalError("Bad quantum numbers detected in module <"+name+">");
}

void PopulationDiagnostic::Report(Diagnostic& diagnostic)
{
	Module::Report(diagnostic);

	ScalarField temp;
	temp.Initialize(*this,owner);

	tw::vec3 ANow(0.0);
	for (auto w : wave)
		ANow += w->VectorPotential(owner->elapsedTime,tw::vec3(0,0,0));

	for (auto ref : refState)
	{
		// Real part
		for (auto cell : InteriorCellRange(*this))
		{
			const tw::vec3 pos = owner->Pos(cell);
			// the following is a gauge transformation to "remove" the uniform A
			const tw::Complex psiNow = (*psi)(cell)*std::exp(-ii*H->qorb*(ANow^pos));
			temp(cell) = real(conj(ref->Amplitude(*H,pos,owner->elapsedTime,0))*psiNow);
		}
		diagnostic.VolumeIntegral("real<"+ref->name+"|psi>",temp,0);

		// Imaginary part
		for (auto cell : InteriorCellRange(*this))
		{
			const tw::vec3 pos = owner->Pos(cell);
			// the following is a gauge transformation to "remove" the uniform A
			const tw::Complex psiNow = (*psi)(cell)*std::exp(-ii*H->qorb*(ANow^pos));
			temp(cell) = imag(conj(ref->Amplitude(*H,pos,owner->elapsedTime,0))*psiNow);
		}
		diagnostic.VolumeIntegral("imag<"+ref->name+"|psi>",temp,0);
	}

	// Relativistic scalar overlaps - enable this later.

	// for (tw::Int s=0;s<refState.size();s++)
	// {
	// 	for (auto cell : InteriorCellRange(*this))
	// 	{
	// 		const tw::vec3 pos = owner->Pos(cell);
	// 		// Feshbach-Villars decomposition, including gauge transformation
	// 		const tw::Complex estate = FV(cell,1.0)*std::exp(-ii*H.qorb*(ANow^pos));
	// 		const tw::Complex pstate = FV(cell,-1.0)*std::exp(-ii*H.qorb*(ANow^pos));
	// 		const tw::Complex psi_ref = refState[s].Amplitude(pos,0.0,0)/root2;
	// 		const tw::Complex chi_ref = refState[s].Amplitude(pos,0.0,1)/root2;
	// 		temp(cell) = real(conj(psi_ref+chi_ref)*estate - conj(psi_ref-chi_ref)*pstate);
	// 	}
	// 	diagnostic.VolumeIntegral("real<ref"+std::to_string(s)+"|psi>",temp,0);
	//
	// 	for (auto cell : InteriorCellRange(*this))
	// 	{
	// 		const tw::vec3 pos = owner->Pos(cell);
	// 		// Feshbach-Villars decomposition, including gauge transformation
	// 		const tw::Complex estate = FV(cell,1.0)*std::exp(-ii*H.qorb*(ANow^pos));
	// 		const tw::Complex pstate = FV(cell,-1.0)*std::exp(-ii*H.qorb*(ANow^pos));
	// 		const tw::Complex psi_ref = refState[s].Amplitude(pos,0.0,0)/root2;
	// 		const tw::Complex chi_ref = refState[s].Amplitude(pos,0.0,1)/root2;
	// 		temp(cell) = imag(conj(psi_ref+chi_ref)*estate - conj(psi_ref-chi_ref)*pstate);
	// 	}
	// 	diagnostic.VolumeIntegral("imag<ref"+std::to_string(s)+"|psi>",temp,0);
	// }
}
