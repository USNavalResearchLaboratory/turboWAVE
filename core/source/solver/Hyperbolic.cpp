module;

#include "tw_includes.h"

export module hyperbolic;
import compute_tool;
import fields;

export struct YeePropagatorPML:ComputeTool
{
	YeePropagatorPML(const std::string& name,MetricSpace *m,Task *tsk);

	#ifdef USE_OPENCL
	virtual ~YeePropagatorPML();
	cl_kernel k_advanceE,k_prepCenteredFields,k_advanceB,k_centeredFields;
	void SetupComputeKernels(Field& F,Field& A,Field& PMLx,Field& PMLy,Field& PMLz,Field& j4);
	#endif

	void AdvanceE(Field& A,Field& PMLx,Field& PMLy,Field& PMLz,Field& j4);
	void AdvanceB(Field& A,Field& PMLx,Field& PMLy,Field& PMLz);
	void PrepCenteredFields(Field& F,Field& A,tw::vec3& E0,tw::vec3& B0);
	void CenteredFields(Field& F,Field& A);
	void UpdateInteriorBoundaryE(Field& A,const ScalarField& conductor);
	void UpdateInteriorBoundaryB(Field& A,const ScalarField& conductor);
	void UpdateExteriorBoundary(Field& A,Field& PMLx,Field& PMLy,Field& PMLz);
};

export struct LorentzPropagator:ComputeTool
{
	LorentzPropagator(const std::string& name,MetricSpace *m,Task *tsk);

	#ifdef USE_OPENCL
	virtual ~LorentzPropagator();
	cl_kernel k_advance,k_swap,k_midstep,k_undoMidstep;
	void SetupComputeKernels(Field& A4,Field& Ao4,Field& j4);
	#endif

	void Advance(Field& A4,Field& Ao4,Field& j4,const tw::Float mult,const tw::Float dt);
	void MidstepEstimate(Field& A4,Field& Ao4);
	void UndoMidstepEstimate(Field& A4,Field& Ao4);
};

YeePropagatorPML::YeePropagatorPML(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	InitializeCLProgram("hyperbolic.cl");
}

#ifdef USE_OPENCL

YeePropagatorPML::~YeePropagatorPML()
{
	clReleaseKernel(k_advanceE);
	clReleaseKernel(k_prepCenteredFields);
	clReleaseKernel(k_advanceB);
	clReleaseKernel(k_centeredFields);
}

void YeePropagatorPML::SetupComputeKernels(Field& F,Field& A,Field& PMLx,Field& PMLy,Field& PMLz,Field& j4)
{
	cl_int err;
	k_advanceE = clCreateKernel(program,"YeeUpdateE",&err);
	k_prepCenteredFields = clCreateKernel(program,"PrepCenteredFields",&err);
	k_advanceB = clCreateKernel(program,"YeeUpdateB",&err);
	k_centeredFields = clCreateKernel(program,"CenteredFields",&err);

	clSetKernelArg(k_advanceE,0,sizeof(cl_mem),&A.computeBuffer);
	clSetKernelArg(k_advanceE,1,sizeof(cl_mem),&j4.computeBuffer);
	clSetKernelArg(k_advanceE,2,sizeof(cl_mem),&PMLx.computeBuffer);
	clSetKernelArg(k_advanceE,3,sizeof(cl_mem),&PMLy.computeBuffer);
	clSetKernelArg(k_advanceE,4,sizeof(cl_mem),&PMLz.computeBuffer);
	clSetKernelArg(k_advanceE,5,sizeof(cl_mem),&space->metricsBuffer);

	clSetKernelArg(k_prepCenteredFields,0,sizeof(cl_mem),&F.computeBuffer);
	clSetKernelArg(k_prepCenteredFields,1,sizeof(cl_mem),&A.computeBuffer);

	clSetKernelArg(k_advanceB,0,sizeof(cl_mem),&A.computeBuffer);
	clSetKernelArg(k_advanceB,1,sizeof(cl_mem),&PMLx.computeBuffer);
	clSetKernelArg(k_advanceB,2,sizeof(cl_mem),&PMLy.computeBuffer);
	clSetKernelArg(k_advanceB,3,sizeof(cl_mem),&PMLz.computeBuffer);
	clSetKernelArg(k_advanceB,4,sizeof(cl_mem),&space->metricsBuffer);

	clSetKernelArg(k_centeredFields,0,sizeof(cl_mem),&F.computeBuffer);
	clSetKernelArg(k_centeredFields,1,sizeof(cl_mem),&A.computeBuffer);
}

void YeePropagatorPML::AdvanceE(Field& A,Field& PMLx,Field& PMLy,Field& PMLz,Field& j4)
{
	space->LocalUpdateProtocol(k_advanceE,task->commandQueue);

	A.ReceiveBoundaryCellsFromComputeBuffer();
	UpdateExteriorBoundary(A,PMLx,PMLy,PMLz);
	A.CopyFromNeighbors(Element(0,5));
	A.SendBoundaryCellsToComputeBuffer();
}

void YeePropagatorPML::AdvanceB(Field& A,Field& PMLx,Field& PMLy,Field& PMLz)
{
 	size_t global_offset[3] = { (size_t)A.Layers(1) , (size_t)A.Layers(2) , (size_t)A.Layers(3) };
	size_t global_work_range[3] = { (size_t)A.UNG(1) , (size_t)A.UNG(2) , (size_t)A.UNG(3) }; // dim+1, or 1
	clEnqueueNDRangeKernel(task->commandQueue,k_advanceB,3,global_offset,global_work_range,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);

	A.ReceiveBoundaryCellsFromComputeBuffer();
	A.UpwardCopy(tw::grid::x,Element(6,11),1);
	A.UpwardCopy(tw::grid::y,Element(6,11),1);
	A.UpwardCopy(tw::grid::z,Element(6,11),1);
	A.SendBoundaryCellsToComputeBuffer();
}

void YeePropagatorPML::PrepCenteredFields(Field& F,Field& A)
{
	// set up indexing for [0...dim] or [0]
	size_t global_work_range[3] = { (size_t)A.UNG(1) , (size_t)A.UNG(2) , (size_t)A.UNG(3) }; // dim+1, or 1
	clEnqueueNDRangeKernel(task->commandQueue,k_prepCenteredFields,3,NULL,global_work_range,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void YeePropagatorPML::CenteredFields(Field& F,Field& A)
{
	// set up indexing for [0...dim] or [0]
	size_t global_work_range[3] = { (size_t)A.UNG(1) , (size_t)A.UNG(2) , (size_t)A.UNG(3) }; // dim+1, or 1
	clEnqueueNDRangeKernel(task->commandQueue,k_centeredFields,3,NULL,global_work_range,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);

	F.ReceiveBoundaryCellsFromComputeBuffer();
	F.DownwardCopy(tw::grid::x,1);
	F.DownwardCopy(tw::grid::y,1);
	F.DownwardCopy(tw::grid::z,1);
	F.ApplyBoundaryCondition();
	F.SendBoundaryCellsToComputeBuffer();
}

#else

void YeePropagatorPML::AdvanceE(Field& A,Field& PMLx,Field& PMLy,Field& PMLz,Field& j4)
{
	//const tw::Int xDim = A.Dim(1);
	//const tw::Int yDim = A.Dim(2);
	const tw::Int zDim = A.Dim(3);

	tw::Float sx,tx,sy,ty;
	std::valarray<tw::Float> sz(space->Num(3)),tz(space->Num(3));
	PMLz.GetStrip(sz,tw::strip(1,PMLz,0,0,0),0);
	PMLz.GetStrip(tz,tw::strip(1,PMLz,0,0,0),1);

	#pragma omp parallel private(sx,sy,tx,ty)
	{
		for (auto v : VectorStripRange<3>(*space,false))
		{
			sx = PMLx(v.dcd1(0),0,0,0); tx = PMLx(v.dcd1(0),0,0,1);
			sy = PMLy(v.dcd2(0),0,0,0); ty = PMLy(v.dcd2(0),0,0,1);
			#pragma omp simd
			for (tw::Int k=1;k<=zDim;k++)
				A(v,k,0) = sy*A(v,k,0) + ty*(A.dfwd(v,k,10,2) + A.dfwd(v,k,11,2) - 0.5*j4(v,k,1));
			#pragma omp simd
			for (tw::Int k=1;k<=zDim;k++)
				A(v,k,1) = sz[k]*A(v,k,1) - tz[k]*(A.dfwd(v,k,8,3) + A.dfwd(v,k,9,3) + 0.5*j4(v,k,1));
			#pragma omp simd
			for (tw::Int k=1;k<=zDim;k++)
				A(v,k,2) = sz[k]*A(v,k,2) + tz[k]*(A.dfwd(v,k,6,3) + A.dfwd(v,k,7,3) - 0.5*j4(v,k,2));
			#pragma omp simd
			for (tw::Int k=1;k<=zDim;k++)
				A(v,k,3) = sx*A(v,k,3) - tx*(A.dfwd(v,k,10,1) + A.dfwd(v,k,11,1) + 0.5*j4(v,k,2));
			#pragma omp simd
			for (tw::Int k=1;k<=zDim;k++)
				A(v,k,4) = sx*A(v,k,4) + tx*(A.dfwd(v,k,8,1) + A.dfwd(v,k,9,1) - 0.5*j4(v,k,3));
			#pragma omp simd
			for (tw::Int k=1;k<=zDim;k++)
				A(v,k,5) = sy*A(v,k,5) - ty*(A.dfwd(v,k,6,2) + A.dfwd(v,k,7,2) + 0.5*j4(v,k,3));
		}
	}
	UpdateExteriorBoundary(A,PMLx,PMLy,PMLz); // must precede MPI
	A.CopyFromNeighbors(Element(0,5));
}

void YeePropagatorPML::AdvanceB(Field& A,Field& PMLx,Field& PMLy,Field& PMLz)
{
	const tw::Int xN1 = A.UNG(1);
	const tw::Int yN1 = A.UNG(2);
	const tw::Int zN1 = A.UNG(3);

	tw::Float sx,tx,sy,ty;
	std::valarray<tw::Float> sz(space->Num(3)),tz(space->Num(3));
	PMLz.GetStrip(sz,tw::strip(1,PMLz,0,0,0),3);
	PMLz.GetStrip(tz,tw::strip(1,PMLz,0,0,0),4);

	#pragma omp parallel for private(sx,sy,tx,ty) collapse(2) schedule(static)
	for (tw::Int i=1;i<=xN1;i++)
		for (tw::Int j=1;j<=yN1;j++)
		{
			tw::xstrip<3> v(*space,i,j,0);
			sx = PMLx(i,0,0,3); tx = PMLx(i,0,0,4);
			sy = PMLy(j,0,0,3); ty = PMLy(j,0,0,4);
			#pragma omp simd
			for (tw::Int k=1;k<=zN1;k++)
				A(v,k,6) = sy*A(v,k,6) - ty*(A.dbak(v,k,4,2) + A.dbak(v,k,5,2));
			#pragma omp simd
			for (tw::Int k=1;k<=zN1;k++)
				A(v,k,7) = sz[k]*A(v,k,7) + tz[k]*(A.dbak(v,k,2,3) + A.dbak(v,k,3,3));
			#pragma omp simd
			for (tw::Int k=1;k<=zN1;k++)
				A(v,k,8) = sz[k]*A(v,k,8) - tz[k]*(A.dbak(v,k,0,3) + A.dbak(v,k,1,3));
			#pragma omp simd
			for (tw::Int k=1;k<=zN1;k++)
				A(v,k,9) = sx*A(v,k,9) + tx*(A.dbak(v,k,4,1) + A.dbak(v,k,5,1));
			#pragma omp simd
			for (tw::Int k=1;k<=zN1;k++)
				A(v,k,10) = sx*A(v,k,10) - tx*(A.dbak(v,k,2,1) + A.dbak(v,k,3,1));
			#pragma omp simd
			for (tw::Int k=1;k<=zN1;k++)
				A(v,k,11) = sy*A(v,k,11) + ty*(A.dbak(v,k,0,2) + A.dbak(v,k,1,2));
		}
	A.UpwardCopy(tw::grid::x,Element(6,11),1);
	A.UpwardCopy(tw::grid::y,Element(6,11),1);
	A.UpwardCopy(tw::grid::z,Element(6,11),1);
}

void YeePropagatorPML::PrepCenteredFields(Field& F,Field& A,tw::vec3& E0,tw::vec3& B0)
{
	F = 0.0;
	add_const_vec<0,1,2>(F,E0);
	add_const_vec<3,4,5>(F,B0);
	AddMulFieldData(F,Element(3),A,Element(6),0.5);
	AddMulFieldData(F,Element(3),A,Element(7),0.5);
	AddMulFieldData(F,Element(4),A,Element(8),0.5);
	AddMulFieldData(F,Element(4),A,Element(9),0.5);
	AddMulFieldData(F,Element(5),A,Element(10),0.5);
	AddMulFieldData(F,Element(5),A,Element(11),0.5);
}

void YeePropagatorPML::CenteredFields(Field& F,Field& A)
{
	AddFieldData(F,Element(0),A,Element(0));
	AddFieldData(F,Element(0),A,Element(1));
	AddFieldData(F,Element(1),A,Element(2));
	AddFieldData(F,Element(1),A,Element(3));
	AddFieldData(F,Element(2),A,Element(4));
	AddFieldData(F,Element(2),A,Element(5));

	AddMulFieldData(F,Element(3),A,Element(6),0.5);
	AddMulFieldData(F,Element(3),A,Element(7),0.5);
	AddMulFieldData(F,Element(4),A,Element(8),0.5);
	AddMulFieldData(F,Element(4),A,Element(9),0.5);
	AddMulFieldData(F,Element(5),A,Element(10),0.5);
	AddMulFieldData(F,Element(5),A,Element(11),0.5);
}

#endif

void YeePropagatorPML::UpdateExteriorBoundary(Field& A,Field& PMLx,Field& PMLy,Field& PMLz)
{
	const tw::Int xDim = A.Dim(1);
	const tw::Int yDim = A.Dim(2);
	const tw::Int zDim = A.Dim(3);

	const tw::vec3 freq(space->dk(1),space->dk(2),space->dk(3));

	tw::Int i,j,k;

	i = A.UNG(1);
	if (task->n1[1]==MPI_PROC_NULL)
		for (j=1;j<=yDim;j++)
			for (k=1;k<=zDim;k++)
			{
				A(i,j,k,0) = PMLy(j,0,0,0)*A(i,j,k,0) + PMLy(j,0,0,1)*freq.y*(A(i,j+1,k,10)+A(i,j+1,k,11)-A(i,j,k,10)-A(i,j,k,11));
				A(i,j,k,1) = PMLz(k,0,0,0)*A(i,j,k,1) - PMLz(k,0,0,1)*freq.z*(A(i,j,k+1,8)+A(i,j,k+1,9)-A(i,j,k,8)-A(i,j,k,9));
			}
	j = A.UNG(2);
	if (task->n1[2]==MPI_PROC_NULL)
		for (i=1;i<=xDim;i++)
			for (k=1;k<=zDim;k++)
			{
				A(i,j,k,2) = PMLz(k,0,0,0)*A(i,j,k,2) + PMLz(k,0,0,1)*freq.z*(A(i,j,k+1,6)+A(i,j,k+1,7)-A(i,j,k,6)-A(i,j,k,7));
				A(i,j,k,3) = PMLx(i,0,0,0)*A(i,j,k,3) - PMLx(i,0,0,1)*freq.x*(A(i+1,j,k,10)+A(i+1,j,k,11)-A(i,j,k,10)-A(i,j,k,11));
			}
	k = A.UNG(3);
	if (task->n1[3]==MPI_PROC_NULL)
		for (i=1;i<=xDim;i++)
			for (j=1;j<=yDim;j++)
			{
				A(i,j,k,4) = PMLx(i,0,0,0)*A(i,j,k,4) + PMLx(i,0,0,1)*freq.x*(A(i+1,j,k,8)+A(i+1,j,k,9)-A(i,j,k,8)-A(i,j,k,9));
				A(i,j,k,5) = PMLy(j,0,0,0)*A(i,j,k,5) - PMLy(j,0,0,1)*freq.y*(A(i,j+1,k,6)+A(i,j+1,k,7)-A(i,j,k,6)-A(i,j,k,7));
			}
}

void YeePropagatorPML::UpdateInteriorBoundaryE(Field& A,const ScalarField& conductor)
{
	// Conductor fills cells shifted back by 1/2
	// Strategy is to zero components along all 12 edges of every hexahedral cell
	// This is highly redundant but very simple
	const tw::Int xN1 = A.UNG(1);
	const tw::Int yN1 = A.UNG(2);
	const tw::Int zN1 = A.UNG(3);
	#pragma omp parallel for collapse(2) schedule(static)
	for (tw::Int i=1;i<=xN1;i++)
		for (tw::Int j=1;j<=yN1;j++)
		{
			tw::xstrip<3> v(*space,i,j,0);
			tw::xstrip<3> vi(*space,i-1,j,0);
			tw::xstrip<3> vj(*space,i,j-1,0);
			tw::xstrip<3> vij(*space,i-1,j-1,0);
			#pragma omp simd
			for (tw::Int k=1;k<=zN1;k++)
			{
				A(v,k,0) *= conductor(v,k,0);
				A(vj,k,0) *= conductor(v,k,0);
				A(v,k-1,0) *= conductor(v,k,0);
				A(vj,k-1,0) *= conductor(v,k,0);
				A(v,k,1) *= conductor(v,k,0);
				A(vj,k,1) *= conductor(v,k,0);
				A(v,k-1,1) *= conductor(v,k,0);
				A(vj,k-1,1) *= conductor(v,k,0);

				A(v,k,2) *= conductor(v,k,0);
				A(vi,k,2) *= conductor(v,k,0);
				A(v,k-1,2) *= conductor(v,k,0);
				A(vi,k-1,2) *= conductor(v,k,0);
				A(v,k,3) *= conductor(v,k,0);
				A(vi,k,3) *= conductor(v,k,0);
				A(v,k-1,3) *= conductor(v,k,0);
				A(vi,k-1,3) *= conductor(v,k,0);

				A(v,k,4) *= conductor(v,k,0);
				A(vi,k,4) *= conductor(v,k,0);
				A(vj,k,4) *= conductor(v,k,0);
				A(vij,k,4) *= conductor(v,k,0);
				A(v,k,5) *= conductor(v,k,0);
				A(vi,k,5) *= conductor(v,k,0);
				A(vj,k,5) *= conductor(v,k,0);
				A(vij,k,5) *= conductor(v,k,0);
			}
		}
}

void YeePropagatorPML::UpdateInteriorBoundaryB(Field& A,const ScalarField& conductor)
{
	// Conductor fills cells shifted back by 1/2
	// Strategy is to zero components normal to all 6 walls of every hexahedral cell
	// This is highly redundant but very simple
	const tw::Int xN1 = A.UNG(1);
	const tw::Int yN1 = A.UNG(2);
	const tw::Int zN1 = A.UNG(3);
	#pragma omp parallel for collapse(2) schedule(static)
	for (tw::Int i=1;i<=xN1;i++)
		for (tw::Int j=1;j<=yN1;j++)
		{
			tw::xstrip<3> v(*space,i,j,0);
			tw::xstrip<3> vi(*space,i-1,j,0);
			tw::xstrip<3> vj(*space,i,j-1,0);
			#pragma omp simd
			for (tw::Int k=1;k<=zN1;k++)
			{
				A(v,k,6) *= conductor(v,k,0);
				A(vi,k,6) *= conductor(v,k,0);
				A(v,k,7) *= conductor(v,k,0);
				A(vi,k,7) *= conductor(v,k,0);

				A(v,k,8) *= conductor(v,k,0);
				A(vj,k,8) *= conductor(v,k,0);
				A(v,k,9) *= conductor(v,k,0);
				A(vj,k,9) *= conductor(v,k,0);

				A(v,k,10) *= conductor(v,k,0);
				A(v,k-1,10) *= conductor(v,k,0);
				A(v,k,11) *= conductor(v,k,0);
				A(v,k-1,11) *= conductor(v,k,0);
			}
		}
}


//////////////////////////////////////////////
// Centered Lorentz Gauge Potentials
// Assumes x-packed fields
//////////////////////////////////////////////


LorentzPropagator::LorentzPropagator(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	InitializeCLProgram("hyperbolic.cl");
}

#ifdef USE_OPENCL

LorentzPropagator::~LorentzPropagator()
{
	clReleaseKernel(k_advance);
	clReleaseKernel(k_swap);
	clReleaseKernel(k_midstep);
	clReleaseKernel(k_undoMidstep);
}

void LorentzPropagator::SetupComputeKernels(Field& A4,Field& Ao4,Field& j4)
{
	cl_int err;
	k_advance = clCreateKernel(program,"LorentzAdvance",&err);
	k_swap = clCreateKernel(program,"LorentzSwap",&err);
	k_midstep = clCreateKernel(program,"LorentzMidstep",&err);
	k_undoMidstep = clCreateKernel(program,"LorentzUndoMidstep",&err);

	clSetKernelArg(k_advance,0,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(k_advance,1,sizeof(cl_mem),&Ao4.computeBuffer);
	clSetKernelArg(k_advance,2,sizeof(cl_mem),&j4.computeBuffer);
	clSetKernelArg(k_advance,3,sizeof(cl_mem),&space->metricsBuffer);

	clSetKernelArg(k_swap,0,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(k_swap,1,sizeof(cl_mem),&Ao4.computeBuffer);

	clSetKernelArg(k_midstep,0,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(k_midstep,1,sizeof(cl_mem),&Ao4.computeBuffer);

	clSetKernelArg(k_undoMidstep,0,sizeof(cl_mem),&A4.computeBuffer);
	clSetKernelArg(k_undoMidstep,1,sizeof(cl_mem),&Ao4.computeBuffer);
}

void LorentzPropagator::Advance(Field& A4,Field& Ao4,Field& j4,const tw::Float mult,const tw::Float dt)
{
	// A4 is the current 4-potential (phi,Ax,Ay,Az)
	// Ao4 is the 4-potential from the previous step
	// J4 is the current 4-current (rho,Jx,Jy,Jz)
	// At present do not update scalar potential
	clSetKernelArg(k_advance,4,sizeof(tw::Float),&mult);
	space->LocalUpdateProtocol(k_advance,task->commandQueue);
	space->PointUpdateProtocol(k_swap,task->commandQueue);
	A4.UpdateGhostCellsInComputeBuffer(Element(1,3));
}

void LorentzPropagator::MidstepEstimate(Field& A4,Field& Ao4)
{
	space->PointUpdateProtocol(k_midstep,task->commandQueue);
}

void LorentzPropagator::UndoMidstepEstimate(Field& A4,Field& Ao4)
{
	space->PointUpdateProtocol(k_undoMidstep,task->commandQueue);
}

#else

void LorentzPropagator::Advance(Field& A4,Field& Ao4,Field& j4,const tw::Float mult,const tw::Float dt)
{
	// A4 is the current 4-potential (phi,Ax,Ay,Az)
	// Ao4 is the 4-potential from the previous step
	// J4 is the current 4-current (rho,Jx,Jy,Jz)
	// At present do not update scalar potential
	for (tw::Int c=1;c<=3;c++)
	{
		#pragma omp parallel firstprivate(c,mult,dt)
		{
			for (auto v : VectorStripRange<1>(*space,false))
				for (tw::Int i=1;i<=A4.Dim(1);i++)
					Ao4(v,i,c) = 2*A4(v,i,c) - Ao4(v,i,c) + dt*dt*(mult*j4(v,i,c) + A4.d2(v,i,c,1) + A4.d2(v,i,c,2) + A4.d2(v,i,c,3));
		}
		#pragma omp parallel firstprivate(c)
		{
			for (auto v : VectorStripRange<1>(*space,true))
				for (tw::Int i=A4.LNG(1);i<=A4.UNG(1);i++)
					std::swap(A4(v,i,c),Ao4(v,i,c));
		}
	}
	A4.CopyFromNeighbors(Element(1,3));
}

void LorentzPropagator::MidstepEstimate(Field& A4,Field& Ao4)
{
	// Put an estimate of the mid-step potential into the current potential
	for (tw::Int c=1;c<=3;c++)
	{
		#pragma omp parallel firstprivate(c)
		{
			for (auto v : VectorStripRange<1>(*space,true))
				for (tw::Int i=A4.LNG(1);i<=A4.UNG(1);i++)
					A4(v,i,c) = 0.5*(Ao4(v,i,c) + A4(v,i,c));
		}
	}
}

void LorentzPropagator::UndoMidstepEstimate(Field& A4,Field& Ao4)
{
	// Take the mid-step estimate back to the full step
	for (tw::Int c=1;c<=3;c++)
	{
		#pragma omp parallel firstprivate(c)
		{
			for (auto v : VectorStripRange<1>(*space,true))
				for (tw::Int i=A4.LNG(1);i<=A4.UNG(1);i++)
					A4(v,i,c) = 2.0*A4(v,i,c) - Ao4(v,i,c);
		}
	}
}

#endif
