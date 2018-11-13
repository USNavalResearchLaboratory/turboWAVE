#include "sim.h"
#include "solidState.h"


////////////////////////////
//                        //
// BOUND ELECTRONS MODULE //
//                        //
////////////////////////////


BoundElectrons::BoundElectrons(Grid* theGrid) : Module(theGrid)
{
	name = "bound";
	typeCode = boundElectrons;
	updateSequencePriority = 2;
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

	if (owner->restarted)
		return;

	for (CellIterator cell(*this,false);cell<cell.end();++cell)
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
	dens.DownwardCopy(zAxis,1);
	R0.DownwardCopy(zAxis,1);
	R1.DownwardCopy(zAxis,1);

	// carry out shift
	for (StripIterator s(*this,3,strongbool::yes);s<s.end();++s)
	{
		tw::vec3 pos = owner->Pos(s,s.Dim()+1);
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

	const tw::Int zLast = owner->n1[3]==MPI_PROC_NULL ? ub[3]-6 : dim[3];

	#pragma omp parallel
	{
		tw::vec3 qEm,r0,r1,FNL;
		tw::vec4 j4;
		tw::Float RR[7];

		for (VectorizingIterator<3> v(*this,false);v<v.end();++v)
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

	if (command=="charge")
	{
		inputString >> word >> q0;
	}
	if (command=="mass")
	{
		inputString >> word >> m0;
	}
	if (command=="basis")
	{
		inputString >> word;
		inputString >> crystalBasis.u.x >> crystalBasis.u.y >> crystalBasis.u.z;
		inputString >> crystalBasis.v.x >> crystalBasis.v.y >> crystalBasis.v.z;
		inputString >> crystalBasis.w.x >> crystalBasis.w.y >> crystalBasis.w.z;
	}
	if (command=="resonance")
	{
		inputString >> word >> resFreq.x >> resFreq.y >> resFreq.z;
	}
	if (command=="damping")
	{
		inputString >> word >> dampFreq.x >> dampFreq.y >> dampFreq.z;
	}
	if (command=="strength")
	{
		inputString >> word >> oscStrength.x >> oscStrength.y >> oscStrength.z;
	}
	if (command=="a1") // eg, a1 = ( 0 , 0 , 0 , 1 , 0 , 0 )
	{
		inputString >> word >> a1[1] >> a1[2] >> a1[3] >> a1[4] >> a1[5] >> a1[6];
	}
	if (command=="a2") // eg, a2 = ( 0 , 0 , 0 , 0 , 1 , 0 )
	{
		inputString >> word >> a2[1] >> a2[2] >> a2[3] >> a2[4] >> a2[5] >> a2[6];
	}
	if (command=="a3") // eg, a3 = ( 0 , 0 , 0 , 0 , 0 , 1 )
	{
		inputString >> word >> a3[1] >> a3[2] >> a3[3] >> a3[4] >> a3[5] >> a3[6];
	}
	if (command=="b") // eg, b = 1
	{
		inputString >> word >> b;
	}
	if (command=="d") // eg, d = 1
	{
		inputString >> word >> d;
	}
	if (command=="theta")
	{
		inputString >> word >> theta;
		theta *= pi/180.0;
	}
	if (command=="phi")
	{
		inputString >> word >> phi;
		phi *= pi/180.0;
	}
}

void BoundElectrons::ReadData(std::ifstream& inFile)
{
	Module::ReadData(inFile);
	inFile.read((char*)&q0,sizeof(tw::Float));
	inFile.read((char*)&m0,sizeof(tw::Float));
	inFile.read((char*)&theta,sizeof(tw::Float));
	inFile.read((char*)&phi,sizeof(tw::Float));
	inFile.read((char*)&crystalBasis,sizeof(crystalBasis));
	inFile.read((char*)&resFreq,sizeof(tw::vec3));
	inFile.read((char*)&dampFreq,sizeof(tw::vec3));
	inFile.read((char*)&oscStrength,sizeof(tw::vec3));
	inFile.read((char *)&a1[1],sizeof(tw::Float)*6);
	inFile.read((char *)&a2[1],sizeof(tw::Float)*6);
	inFile.read((char *)&a3[1],sizeof(tw::Float)*6);
	inFile.read((char *)&b,sizeof(tw::Float));
	inFile.read((char *)&d,sizeof(tw::Float));

	R0.ReadData(inFile);
	R1.ReadData(inFile);
	dens.ReadData(inFile);
}

void BoundElectrons::WriteData(std::ofstream& outFile)
{
	Module::WriteData(outFile);
	outFile.write((char*)&q0,sizeof(tw::Float));
	outFile.write((char*)&m0,sizeof(tw::Float));
	outFile.write((char*)&theta,sizeof(tw::Float));
	outFile.write((char*)&phi,sizeof(tw::Float));
	outFile.write((char*)&crystalBasis,sizeof(crystalBasis));
	outFile.write((char*)&resFreq,sizeof(tw::vec3));
	outFile.write((char*)&dampFreq,sizeof(tw::vec3));
	outFile.write((char*)&oscStrength,sizeof(tw::vec3));
	outFile.write((char *)&a1[1],sizeof(tw::Float)*6);
	outFile.write((char *)&a2[1],sizeof(tw::Float)*6);
	outFile.write((char *)&a3[1],sizeof(tw::Float)*6);
	outFile.write((char *)&b,sizeof(tw::Float));
	outFile.write((char *)&d,sizeof(tw::Float));

	R0.WriteData(outFile);
	R1.WriteData(outFile);
	dens.WriteData(outFile);
}

void BoundElectrons::StartDiagnostics()
{
	#ifdef USE_OPENCL
	R0.ReceiveFromComputeBuffer();
	R1.ReceiveFromComputeBuffer();
	#endif
}

void BoundElectrons::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "bound-KE bound-PE ";
}

void BoundElectrons::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	tw::Int i,j,k;
	tw::Float kinetic,potential,number;
	tw::vec3 pos,r0,r1,vel;

	kinetic = 0.0;
	potential = 0.0;

	tw::Int x0,x1,y0,y1,z0,z1;
	theRgn.GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
	for (k=z0;k<=z1;k++)
		for (j=y0;j<=y1;j++)
			for (i=x0;i<=x1;i++)
			{
				pos = owner->Pos(i,j,k);
				if (theRgn.Inside(pos,*owner))
				{
					number = dens(i,j,k)*owner->dS(i,j,k,0);
					r0 = oscStrength*tw::vec3(R0(i,j,k,0),R0(i,j,k,1),R0(i,j,k,2));
 					r1 = oscStrength*tw::vec3(R1(i,j,k,0),R1(i,j,k,1),R1(i,j,k,2));
                    vel = (r1-r0)/dt; // in crystal basis
					kinetic += 0.5*m0*Norm(vel)*number;
					potential += 0.5*m0*Norm(resFreq*r1)*number;
					//diss += 2.0*m0*((dampFreq*vel)^vel)*number;
				}
			}

	cols.push_back(kinetic); avg.push_back(false);
	cols.push_back(potential); avg.push_back(false);
}

void BoundElectrons::BoxDiagnosticHeader(GridDataDescriptor *box)
{
	owner->WriteBoxDataHeader(name,box);
	owner->WriteBoxDataHeader(name+"_x",box);
	owner->WriteBoxDataHeader(name+"_y",box);
	owner->WriteBoxDataHeader(name+"_z",box);
}

void BoundElectrons::BoxDiagnose(GridDataDescriptor *box)
{
	owner->WriteBoxData(name,box,&dens(0,0,0),dens.Stride());
	owner->WriteBoxData(name+"_x",box,&R1(0,0,0,0),R1.Stride());
	owner->WriteBoxData(name+"_y",box,&R1(0,0,0,1),R1.Stride());
	owner->WriteBoxData(name+"_z",box,&R1(0,0,0,2),R1.Stride());
}

void BoundElectrons::PointDiagnosticHeader(std::ofstream& outFile)
{
}

void BoundElectrons::PointDiagnose(std::ofstream& outFile,const weights_3D& w)
{
}
