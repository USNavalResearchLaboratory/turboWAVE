#include "simulation.h"
#include "solidState.h"

////////////////////////////
//                        //
// BOUND ELECTRONS MODULE //
//                        //
////////////////////////////


BoundElectrons::BoundElectrons(const std::string& name,Simulation* sim) : Module(name,sim)
{
	updateSequencePriority = tw::priority::source;
	subSequencePriority = 1;
	q0 = -1.0;
	m0 = 1.0;
	resFreq = 10.0;
	dampFreq = 0.0;
	oscStrength = 1.0;
	a1.resize(7); a1 = 0.0;
	a2.resize(7); a2 = 0.0;
	a3.resize(7); a3 = 0.0;
	b = d = 0.0;
	crystalBasis.u = tw::vec3(1,0,0);
	crystalBasis.v = tw::vec3(0,1,0);
	crystalBasis.w = tw::vec3(0,0,1);
	theta = 0.0;
	phi = 0.0;

	dens.Initialize(*this,owner);
	R0.Initialize(3,*this,owner);
	R1.Initialize(3,*this,owner);
	packet.resize(38);

	#ifdef USE_OPENCL
	cl_int err;
	InitializeCLProgram("particles.cl");
	R0.InitializeComputeBuffer();
	R1.InitializeComputeBuffer();
	dens.InitializeComputeBuffer();
	packet_buffer = clCreateBuffer(owner->context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,sizeof(tw::Float)*packet.size(),&packet[0],&err);
	k_update = clCreateKernel(program,"PushBound",&err);
	// setting kernel arguments has to wait until all modules create compute buffers
	// hence, do it in Initialize() below
	#endif

	directives.Add("charge",new tw::input::Float(&q0),false);
	directives.Add("mass",new tw::input::Float(&m0),false);
	directives.Add("resonance",new tw::input::Vec3(&resFreq));
	directives.Add("damping",new tw::input::Vec3(&dampFreq),false);
	directives.Add("strength",new tw::input::Vec3(&oscStrength),false);
	directives.Add("a1",new tw::input::Numbers<tw::Float>(&a1[1],6),false);
	directives.Add("a2",new tw::input::Numbers<tw::Float>(&a2[1],6),false);
	directives.Add("a3",new tw::input::Numbers<tw::Float>(&a3[1],6),false);
	directives.Add("b",new tw::input::Float(&b),false);
	directives.Add("d",new tw::input::Float(&d),false);
	directives.Add("theta",new tw::input::Float(&theta),false);
	directives.Add("phi",new tw::input::Float(&phi),false);
	directives.Add("basis",new tw::input::Custom,false);
}

BoundElectrons::~BoundElectrons()
{
	#ifdef USE_OPENCL
	clReleaseMemObject(packet_buffer);
	clReleaseKernel(k_update);
	#endif
}

void BoundElectrons::Initialize()
{
	tw::Int i,s;
	tw::vec3 pos;

	Module::Initialize();
	Normalize(crystalBasis);
	crystalBasis.u.RotateZ(phi);
	crystalBasis.v.RotateZ(phi);
	crystalBasis.w.RotateZ(phi);
	crystalBasis.u.RotateY(theta);
	crystalBasis.v.RotateY(theta);
	crystalBasis.w.RotateY(theta);

	for (auto cell : InteriorCellRange(*this))
	{
		for (s=0;s<3;s++)
		{
			R0(cell,s) = 0.0;
			R1(cell,s) = 0.0;
		}
		pos = owner->Pos(cell);
		dens(cell) = 0.0;
		for (s=0;s<profile.size();s++)
			dens(cell) += profile[s]->GetValue(pos,*owner);
	}

	dens.CopyFromNeighbors();

	// load packet containing crystal parameters
	for (i=0;i<3;i++)
	{
		packet[i] = resFreq[i];
		packet[i+3] = dampFreq[i];
		packet[i+6] = oscStrength[i];
		packet[i+29] = crystalBasis.u[i];
		packet[i+32] = crystalBasis.v[i];
		packet[i+35] = crystalBasis.w[i];
	}
	packet[27] = b;
	packet[28] = d;
	for (i=0;i<6;i++)
	{
		packet[i+9] = a1[i+1];
		packet[i+15] = a2[i+1];
		packet[i+21] = a3[i+1];
	}

	#ifdef USE_OPENCL

	dens.SendToComputeBuffer();
	clEnqueueWriteBuffer(owner->commandQueue,packet_buffer,CL_TRUE,0,sizeof(tw::Float)*packet.size(),&packet[0],  0,NULL,NULL);
	clFinish(owner->commandQueue);

	clSetKernelArg(k_update,0,sizeof(cl_mem),&dens.computeBuffer);
	clSetKernelArg(k_update,1,sizeof(cl_mem),&R0.computeBuffer);
	clSetKernelArg(k_update,2,sizeof(cl_mem),&R1.computeBuffer);
	clSetKernelArg(k_update,3,sizeof(cl_mem),&EM->computeBuffer);
	clSetKernelArg(k_update,4,sizeof(cl_mem),&sources->computeBuffer);
	clSetKernelArg(k_update,5,sizeof(cl_mem),&packet_buffer);
	clSetKernelArg(k_update,6,sizeof(cl_mem),&owner->metricsBuffer);
	clSetKernelArg(k_update,7,sizeof(q0),&q0);
	clSetKernelArg(k_update,8,sizeof(m0),&m0);
	clSetKernelArg(k_update,9,sizeof(dt),&dt);

	#endif
}

bool BoundElectrons::InspectResource(void* resource,const std::string& description)
{
	if (description=="electrostatic:E")
	{
		ESField = (Vec3Field*)resource;
		return true;
	}

	if (description=="electromagnetic:F")
	{
		EM = (Field*)resource;
		return true;
	}

	if (description=="electromagnetic:sources")
	{
		sources = (Field*)resource;
		return true;
	}

	if (description=="kinetics:fixed")
	{
		fixed = (ScalarField*)resource;
		return true;
	}

	if (description=="laser:F")
	{
		laser = (Field*)resource;
		return true;
	}

	if (description=="laser:chi")
	{
		chi = (ComplexField*)resource;
		return true;
	}

	if (description=="laser:carrierFrequency")
	{
		carrierFrequency = (tw::Float*)resource;
		return true;
	}

	if (description=="laser:circular")
	{
		circular = (bool*)resource;
		return true;
	}

	return false;
}

void BoundElectrons::MoveWindow()
{
	// must prepare before shift
	dens.DownwardCopy(tw::grid::z,1);
	R0.DownwardCopy(tw::grid::z,1);
	R1.DownwardCopy(tw::grid::z,1);

	// carry out shift
	for (auto s : StripRange(*this,3,strongbool::yes))
	{
		tw::vec3 pos = owner->Pos(s,Dim(s.Axis())+1);
		tw::Float incomingMaterial = 0.0;
		for (tw::Int p=0;p<profile.size();p++)
			incomingMaterial += profile[p]->GetValue(pos,*owner);
		dens.Shift(s,-1,incomingMaterial);
		R0.Shift(s,-1,0.0);
		R1.Shift(s,-1,0.0);
	}
}

#ifdef USE_OPENCL
void BoundElectrons::Update()
{
	owner->LocalUpdateProtocol(k_update);
	// J4 boundary conditions and message passing treated in field solver
}
#else
void BoundElectrons::Update()
{
	// Solve Dt^2(R) + 2*nu*Dt(R) + omega^2*R = FNL - E
	// FNL = -a*R*R + b*(R.R)*R (really F/m)
	// At start, E is known at t = 0, R0 at t = -1, R1 at t = 0

	const tw::Int zLast = owner->n1[3]==MPI_PROC_NULL ? ufg[3]-6 : dim[3];

	#pragma omp parallel
	{
		tw::vec3 qEm,r0,r1,FNL;
		tw::vec4 j4;
		tw::Float RR[7];

		for (auto v : VectorStripRange<3>(*this,false))
			for (tw::Int k=1;k<=zLast;k++)
				if (dens(v,k)!=0.0)
				{
					qEm = (q0/m0)*tw::vec3(EM->sfwd(v,k,0,1),EM->sfwd(v,k,1,2),EM->sfwd(v,k,2,3));
					crystalBasis.ExpressInBasis(&qEm);

					r0 = tw::vec3(R0(v,k,0),R0(v,k,1),R0(v,k,2));
					r1 = tw::vec3(R1(v,k,0),R1(v,k,1),R1(v,k,2));

					RR[1] = sqr(r1.x);
					RR[2] = sqr(r1.y);
					RR[3] = sqr(r1.z);
					RR[4] = 2.0 * r1.y * r1.z;
					RR[5] = 2.0 * r1.x * r1.z;
					RR[6] = 2.0 * r1.x * r1.y;
					FNL.x = -(a1[1]*RR[1] + a1[2]*RR[2] + a1[3]*RR[3] + a1[4]*RR[4] + a1[5]*RR[5] + a1[6]*RR[6]);
					FNL.y = -(a2[1]*RR[1] + a2[2]*RR[2] + a2[3]*RR[3] + a2[4]*RR[4] + a2[5]*RR[5] + a2[6]*RR[6]);
					FNL.z = -(a3[1]*RR[1] + a3[2]*RR[2] + a3[3]*RR[3] + a3[4]*RR[4] + a3[5]*RR[5] + a3[6]*RR[6]);
					FNL += b*(r1 ^ r1)*r1 + d*(r1 ^ r1)*(r1 ^ r1)*r1;

					for (tw::Int c=0;c<3;c++)
					{
						R0(v,k,c) = r1[c];
						R1(v,k,c) = r1[c]*(2.0 + 2.0*dampFreq[c]*dt - resFreq[c]*resFreq[c]*dt*dt) - r0[c] + qEm[c]*dt*dt + FNL[c]*dt*dt;
						R1(v,k,c) /= 1.0 + 2.0*dampFreq[c]*dt;
					}
					r0 = oscStrength*tw::vec3(R0(v,k,0),R0(v,k,1),R0(v,k,2));
					r1 = oscStrength*tw::vec3(R1(v,k,0),R1(v,k,1),R1(v,k,2));
					crystalBasis.ExpressInStdBasis(&r0);
					crystalBasis.ExpressInStdBasis(&r1);

					// charge in high-side adjacent cell last step
					const tw::vec3 Q0 = 0.5*q0*dens(v,k)*owner->dS(v,k,0)*r0*freq;
					// charge in high-side adjacent cell this step
					const tw::vec3 Q1 = 0.5*q0*dens(v,k)*owner->dS(v,k,0)*r1*freq;
					// current in either cell wall
					const tw::vec3 I3 = (Q1 - Q0)/dt;

					tw::Int i,j,kout;
					v.Decode(k,&i,&j,&kout);

					// Deposit cell charge using small displacement model
					#pragma omp atomic update
					(*sources)(i-1,j,k,0) -= Q1.x;
					#pragma omp atomic update
					(*sources)(i,j-1,k,0) -= Q1.y;
					(*sources)(i,j,k-1,0) -= Q1.z;
					(*sources)(i,j,k+1,0) += Q1.z;
					#pragma omp atomic update
					(*sources)(i,j+1,k,0) += Q1.y;
					#pragma omp atomic update
					(*sources)(i+1,j,k,0) += Q1.x;

					// Deposit wall current using small displacement model
					(*sources)(i,j,k,1) += I3.x;
					#pragma omp atomic update
					(*sources)(i+1,j,k,1) += I3.x;
					(*sources)(i,j,k,2) += I3.y;
					#pragma omp atomic update
					(*sources)(i,j+1,k,2) += I3.y;
					(*sources)(i,j,k,3) += I3.z;
					(*sources)(i,j,k+1,3) += I3.z;
				}
	}
}
#endif

void BoundElectrons::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;

	Module::ReadInputFileDirective(inputString,command);

	if (command=="basis")
	{
		inputString >> word;
		inputString >> crystalBasis.u.x >> crystalBasis.u.y >> crystalBasis.u.z;
		inputString >> crystalBasis.v.x >> crystalBasis.v.y >> crystalBasis.v.z;
		inputString >> crystalBasis.w.x >> crystalBasis.w.y >> crystalBasis.w.z;
	}
}

void BoundElectrons::ReadCheckpoint(std::ifstream& inFile)
{
	Module::ReadCheckpoint(inFile);
	R0.ReadCheckpoint(inFile);
	R1.ReadCheckpoint(inFile);
	dens.ReadCheckpoint(inFile);
}

void BoundElectrons::WriteCheckpoint(std::ofstream& outFile)
{
	Module::WriteCheckpoint(outFile);
	R0.WriteCheckpoint(outFile);
	R1.WriteCheckpoint(outFile);
	dens.WriteCheckpoint(outFile);
}

void BoundElectrons::StartDiagnostics()
{
	#ifdef USE_OPENCL
	R0.ReceiveFromComputeBuffer();
	R1.ReceiveFromComputeBuffer();
	#endif
}

void BoundElectrons::Report(Diagnostic& diagnostic)
{
	Module::Report(diagnostic);

	ScalarField temp;
	temp.Initialize(*this,owner);

	for (auto cell : InteriorCellRange(*this))
	{
		const tw::vec3 vel = oscStrength*(R1.Vec3(cell,0) - R0.Vec3(cell,0))/dt;
		temp(cell) = 0.5*m0*dens(cell)*Norm(vel);
	}
	diagnostic.VolumeIntegral("bound-KE",temp,0);

	for (auto cell : InteriorCellRange(*this))
		temp(cell) = 0.5*m0*dens(cell)*Norm(resFreq*R1.Vec3(cell,0));
	diagnostic.VolumeIntegral("bound-PE",temp,0);

	diagnostic.Field(name,dens,0);
	diagnostic.Field(name+"_x",R1,0);
	diagnostic.Field(name+"_y",R1,1);
	diagnostic.Field(name+"_z",R1,2);
}
