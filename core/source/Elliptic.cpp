#include "definitions.h"
#include "tasks.h"
#include "ctools.h"
#include "3dmath.h"
#include "metricSpace.h"
#include "3dfields.h"
#include "region.h"
#include "numerics.h"
#include "computeTool.h"
#include "elliptic.h"

//////////////////////////////////
//                              //
// ELLIPTICAL SOLVER BASE CLASS //
//                              //
//////////////////////////////////


EllipticSolver::EllipticSolver(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	coeff = NULL;
	x0 = x1 = y0 = y1 = z0 = z1 = none;
	gammaBeam = 1.0;
}

void EllipticSolver::Initialize()
{
	ComputeTool::Initialize();

	lbc.resize(space->Num(1));
	rbc.resize(space->Num(1));
	lbc_t.resize(space->Num(1));
	rbc_t.resize(space->Num(1));
	lbc = 0.0;
	rbc = 0.0;
	lbc_t = 0.0;
	rbc_t = 0.0;
}

void EllipticSolver::FormOperatorStencil(std::valarray<tw::Float>& D,tw::Int i,tw::Int j,tw::Int k)
{
	tw::Float kx = space->Dim(1)>1 ? 1.0 : 0.0;
	tw::Float ky = space->Dim(2)>1 ? 1.0 : 0.0;
	tw::Float kz = space->Dim(3)>1 ? 1.0 : 0.0;
	tw::Float Vol = space->dS(i,j,k,0);
	D[1] = kx*(space->dS(i,j,k,1)/Vol) / space->dl(i,j,k,1);
	D[2] = kx*(space->dS(i+1,j,k,1)/Vol) / space->dl(i+1,j,k,1);
	D[3] = ky*(space->dS(i,j,k,2)/Vol) / space->dl(i,j,k,2);
	D[4] = ky*(space->dS(i,j+1,k,2)/Vol) / space->dl(i,j+1,k,2);
	D[5] = kz*(space->dS(i,j,k,3)/Vol) / space->dl(i,j,k,3);
	D[6] = kz*(space->dS(i,j,k+1,3)/Vol) / space->dl(i,j,k+1,3);
	if (coeff!=NULL)
	{
		D[1] *= 0.5*((*coeff)(i-1,j,k) + (*coeff)(i,j,k));
		D[2] *= 0.5*((*coeff)(i,j,k) + (*coeff)(i+1,j,k));
		D[3] *= 0.5*((*coeff)(i,j-1,k) + (*coeff)(i,j,k));
		D[4] *= 0.5*((*coeff)(i,j,k) + (*coeff)(i,j+1,k));
		D[5] *= 0.5*((*coeff)(i,j,k-1) + (*coeff)(i,j,k));
		D[6] *= 0.5*((*coeff)(i,j,k) + (*coeff)(i,j,k+1));
	}
	D[0] = -(D[1] + D[2] + D[3] + D[4] + D[5] + D[6]);
}

void EllipticSolver::FixPotential(ScalarField& phi,Region* theRegion,const tw::Float& thePotential)
{
}

void EllipticSolver::SetCoefficients(ScalarField *coefficients)
{
	coeff = coefficients;
}

void EllipticSolver::SetBoundaryConditions(ScalarField& phi)
{
	phi.SetBoundaryConditions(xAxis,x0,x1);
	phi.SetBoundaryConditions(yAxis,y0,y1);
	phi.SetBoundaryConditions(zAxis,z0,z1);
}

void EllipticSolver::SetBoundaryConditions(boundarySpec x0,boundarySpec x1,boundarySpec y0,boundarySpec y1,boundarySpec z0,boundarySpec z1)
{
	this->x0 = x0;
	this->y0 = y0;
	this->z0 = z0;
	this->x1 = x1;
	this->y1 = y1;
	this->z1 = z1;
}

void EllipticSolver::SaveBoundaryConditions()
{
	x0s = x0;
	x1s = x1;
	y0s = y0;
	y1s = y1;
	z0s = z0;
	z1s = z1;
}

void EllipticSolver::RestoreBoundaryConditions()
{
	x0 = x0s;
	x1 = x1s;
	y0 = y0s;
	y1 = y1s;
	z0 = z0s;
	z1 = z1s;
}

void EllipticSolver::ZeroModeGhostCellValues(tw::Float *phi0,tw::Float *phiN1,ScalarField& source,tw::Float mul)
{
	// specialized for uniform grid in z
	// assumes source has already been resolved into transverse eigenmodes, and k=0 mode has i=j=1
	// computes ghost cell potentials for k=0 mode by integrating against discretized Green's function
	// the assumption is that either z0 or z1 or both are "natural"
	// if either z0 or z1 (but not both) are dirichletCell, first assume both are natural, then shift result appropriately

	tw::Int k,kg,kN1;

	*phi0 = 0.0;
	*phiN1 = 0.0;
	kN1 = task->globalCells[3]+1;

	for (k=1;k<=source.Dim(3);k++)
	{
		kg = task->GlobalCellIndex(k,3);
		*phi0 += 0.5*mul*sqr(dz(source))*source(1,1,k)*fabs(tw::Float(kg));
		*phiN1 += 0.5*mul*sqr(dz(source))*source(1,1,k)*fabs(tw::Float(kN1-kg));
	}

	task->strip[3].AllSum(phi0,phi0,sizeof(tw::Float),0);
	task->strip[3].AllSum(phiN1,phiN1,sizeof(tw::Float),0);

	if (z0==dirichletCell)
	{
		*phiN1 -= *phi0;
		*phi0 = 0.0;
	}

	if (z1==dirichletCell)
	{
		*phi0 -= *phiN1;
		*phiN1 = 0.0;
	}
}

void EllipticSolver::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;
	ComputeTool::ReadInputFileDirective(inputString,command);
	if (command=="poisson") // eg, poisson boundary condition z = (open,open)
	{
		boundarySpec lowBC,highBC;
		std::string whichAxis,theBC;
		inputString >> word >> word >> whichAxis;
		if (whichAxis=="=")
			whichAxis = "z";
		else
			inputString >> word;

		inputString >> word;
		if (word=="none")
			lowBC = none;
		if (word=="open")
			lowBC = natural;
		if (word=="dirichlet")
			lowBC = dirichletCell;
		if (word=="neumann")
			lowBC = neumannWall;

		inputString >> word;
		if (word=="none")
			highBC = none;
		if (word=="open")
			highBC = natural;
		if (word=="dirichlet")
			highBC = dirichletCell;
		if (word=="neumann")
			highBC = neumannWall;

		if (whichAxis=="x") { x0 = lowBC; x1 = highBC; }
		if (whichAxis=="y") { y0 = lowBC; y1 = highBC; }
		if (whichAxis=="z") { z0 = lowBC; z1 = highBC; }
	}
}

void EllipticSolver::ReadData(std::ifstream& inFile)
{
	ComputeTool::ReadData(inFile);
	inFile.read((char*)&x0,sizeof(x0));
	inFile.read((char*)&x1,sizeof(x1));
	inFile.read((char*)&y0,sizeof(y0));
	inFile.read((char*)&y1,sizeof(y1));
	inFile.read((char*)&z0,sizeof(z0));
	inFile.read((char*)&z1,sizeof(z1));
	inFile.read((char*)&gammaBeam,sizeof(gammaBeam));
}

void EllipticSolver::WriteData(std::ofstream& outFile)
{
	ComputeTool::WriteData(outFile);
	outFile.write((char*)&x0,sizeof(x0));
	outFile.write((char*)&x1,sizeof(x1));
	outFile.write((char*)&y0,sizeof(y0));
	outFile.write((char*)&y1,sizeof(y1));
	outFile.write((char*)&z0,sizeof(z0));
	outFile.write((char*)&z1,sizeof(z1));
	outFile.write((char*)&gammaBeam,sizeof(gammaBeam));
}


/////////////////////////////////
//                             //
//    1D ELLIPTICAL SOLVER     //
//                             //
/////////////////////////////////



EllipticSolver1D::EllipticSolver1D(const std::string& name,MetricSpace *m,Task *tsk) : EllipticSolver(name,m,tsk)
{
	typeCode = tw::tool_type::ellipticSolver1D;
	globalIntegrator = NULL;
}

EllipticSolver1D::~EllipticSolver1D()
{
	if (globalIntegrator!=NULL)
		delete globalIntegrator;
}

void EllipticSolver1D::Initialize()
{
	EllipticSolver::Initialize();

	// SET UP GLOBAL INTEGRATOR

	if (globalIntegrator!=NULL)
		delete globalIntegrator;
	if (task->globalCells[1]>1)
		globalIntegrator = new GlobalIntegrator<tw::Float>(&task->strip[1],space->Dim(3)*space->Dim(2),space->Dim(1));
	if (task->globalCells[2]>1)
		globalIntegrator = new GlobalIntegrator<tw::Float>(&task->strip[2],space->Dim(1)*space->Dim(3),space->Dim(2));
	if (task->globalCells[3]>1)
		globalIntegrator = new GlobalIntegrator<tw::Float>(&task->strip[3],space->Dim(1)*space->Dim(2),space->Dim(3));
}

void EllipticSolver1D::Solve(ScalarField& phi,ScalarField& source,tw::Float mul)
{
	// solve div(coeff*grad(phi)) = mul*source
	// requires 1D grid

	axisSpec axis;
	tw::Int s,sDim,ax,di,dj,dk;
	std::valarray<tw::Float> D(7);
	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	di = dj = dk = 0;
	if (task->globalCells[1]>1)
	{
		axis = xAxis;
		sDim = xDim;
		di = 1;
	}
	if (task->globalCells[2]>1)
	{
		axis = yAxis;
		sDim = yDim;
		dj = 1;
	}
	if (task->globalCells[3]>1)
	{
		axis = zAxis;
		sDim = zDim;
		dk = 1;
	}
	ax = naxis(axis);

	std::valarray<tw::Float> T1(sDim),T2(sDim),T3(sDim),src(sDim),ans(sDim);

	for (s=1;s<=sDim;s++)
	{
		FormOperatorStencil(D,1-di+di*s,1-dj+dj*s,1-dk+dk*s);
		T1[s-1] = D[1]+D[3]+D[5]; // ignorable contributions go away as required
		T2[s-1] = D[0];
		T3[s-1] = D[2]+D[4]+D[6];
		src[s-1] = mul*source(s*di,s*dj,s*dk);
	}
	if (task->n0[ax]==MPI_PROC_NULL)
		phi.AdjustTridiagonalForBoundaries(axis,lowSide,&T1[0],&T2[0],&T3[0],&src[0],0.0);
	if (task->n1[ax]==MPI_PROC_NULL)
		phi.AdjustTridiagonalForBoundaries(axis,highSide,&T1[sDim-1],&T2[sDim-1],&T3[sDim-1],&src[sDim-1],0.0);
	TriDiagonal<tw::Float,tw::Float>(ans,src,T1,T2,T3);
	for (s=1;s<=sDim;s++)
		phi(s*di,s*dj,s*dk) = ans[s-1];
	globalIntegrator->SetMatrix(0,T1,T2,T3);
	globalIntegrator->SetData(0,&phi(0,0,0),1);

	globalIntegrator->Parallelize();

	phi.CopyFromNeighbors();
	phi.ApplyBoundaryCondition();
	/*if (task->n0[ax]==MPI_PROC_NULL)
		phi(0,0,0) = 0.0;
	if (task->n1[ax]==MPI_PROC_NULL)
		phi(di*(sDim+1),dj*(sDim+1),dk*(sDim+1)) = 0.0;*/
}



/////////////////////////////////
//                             //
// ITERATIVE ELLIPTICAL SOLVER //
//                             //
/////////////////////////////////



IterativePoissonSolver::IterativePoissonSolver(const std::string& name,MetricSpace *m,Task *tsk) : EllipticSolver(name,m,tsk)
{
	typeCode = tw::tool_type::iterativePoissonSolver;
	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);
	mask1 = new char[xDim*yDim*zDim];
	mask2 = new char[xDim*yDim*zDim];
	maxIterations = 1000;
	tolerance = 1e-8;
	overrelaxation = 1.9;
	minimumNorm = tw::small_pos;
	iterationsPerformed = 0;
	normSource = 0.0;
	normResidualAchieved = 0.0;
	overrelaxationChange = 0.001;
	InitializeCLProgram("elliptic.cl");
}

IterativePoissonSolver::~IterativePoissonSolver()
{
	delete mask1;
	delete mask2;
}

void IterativePoissonSolver::Initialize()
{
	EllipticSolver::Initialize();

	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	tw::Int i,j,k,n;
	bool redSquare = true;

	for (k=1;k<=zDim;k++)
	{
		for (j=1;j<=yDim;j++)
		{
			for (i=1;i<=xDim;i++)
			{
				n = (i-1) + (j-1)*xDim + (k-1)*xDim*yDim;
				mask1[n] = redSquare ? 1 : 0;
				mask2[n] = redSquare ? 0 : 1;
				if (xDim>1)
					redSquare = !redSquare;
			}
			if (yDim>1)
				redSquare = !redSquare;
		}
		if (zDim>1)
			redSquare = !redSquare;
	}
}

void IterativePoissonSolver::FixPotential(ScalarField& phi,Region* theRegion,const tw::Float& thePotential)
{
	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	tw::Int i,j,k,n;

	for (k=1;k<=zDim;k++)
		for (j=1;j<=yDim;j++)
			for (i=1;i<=xDim;i++)
			{
				n = (i-1) + (j-1)*xDim + (k-1)*xDim*yDim;
				if (theRegion->Inside(space->Pos(i,j,k),*space))
				{
					mask1[n] = 0;
					mask2[n] = 0;
					phi(i,j,k) = thePotential;
				}
			}
	phi.CopyFromNeighbors();
	phi.ApplyBoundaryCondition();
}

#ifdef USE_OPENCL

void IterativePoissonSolver::Solve(ScalarField& phi,ScalarField& source,tw::Float mul)
{
	// solve div(coeff*grad(phi)) = mul*source

	cl_int err;
	tw::Int iter,i,j,k;
	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	tw::Float normResidual;
	ScalarField residual;
	residual.Initialize(*space,task);
	residual.InitializeComputeBuffer();

	normSource = 0.0;
	for (k=1;k<=zDim;k++)
		for (j=1;j<=yDim;j++)
			for (i=1;i<=xDim;i++)
				normSource += fabs(mul*source(i,j,k));
	task->strip[0].AllSum(&normSource,&normSource,sizeof(tw::Float),0);
	normSource += minimumNorm;

	cl_kernel SORIterationKernel = clCreateKernel(program,coeff==NULL ? "SORIterationPoisson" : "SORIterationGeneral",&err);

	cl_mem mask1_buffer = clCreateBuffer(task->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(char)*xDim*yDim*zDim,mask1,&err);
	cl_mem mask2_buffer = clCreateBuffer(task->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,sizeof(char)*xDim*yDim*zDim,mask2,&err);

	clSetKernelArg(SORIterationKernel,0,sizeof(cl_mem),&phi.computeBuffer);
	clSetKernelArg(SORIterationKernel,1,sizeof(cl_mem),&source.computeBuffer);
	clSetKernelArg(SORIterationKernel,2,sizeof(cl_mem),&residual.computeBuffer);
	clSetKernelArg(SORIterationKernel,3,sizeof(cl_mem),&mask1_buffer);
	clSetKernelArg(SORIterationKernel,4,sizeof(mul),&mul);
	clSetKernelArg(SORIterationKernel,5,sizeof(overrelaxation),&overrelaxation);
	clSetKernelArg(SORIterationKernel,6,sizeof(cl_mem),&space->metricsBuffer);
	if (coeff!=NULL) clSetKernelArg(SORIterationKernel,7,sizeof(cl_mem),&coeff->computeBuffer);

	for (iter=0;iter<maxIterations;iter++)
	{
		// Update the interior (red squares)
		clSetKernelArg(SORIterationKernel,3,sizeof(cl_mem),&mask1_buffer);
		space->LocalUpdateProtocol(SORIterationKernel,task->commandQueue);

		phi.UpdateGhostCellsInComputeBuffer();

		// Update the interior (black squares)
		clSetKernelArg(SORIterationKernel,3,sizeof(cl_mem),&mask2_buffer);
		space->LocalUpdateProtocol(SORIterationKernel,task->commandQueue);

		phi.UpdateGhostCellsInComputeBuffer();

		// Compute norm of the residual

		normResidual = residual.DestructiveNorm1ComputeBuffer();
		task->strip[0].AllSum(&normResidual,&normResidual,sizeof(tw::Float),0);

		if (normResidual <= tolerance*normSource)
			break;
	}
	phi.ReceiveFromComputeBuffer();
	normResidualAchieved = normResidual;
	iterationsPerformed = iter;

	clReleaseMemObject(mask1_buffer);
	clReleaseMemObject(mask2_buffer);
	clReleaseKernel(SORIterationKernel);
}

#else

void IterativePoissonSolver::Solve(ScalarField& phi,ScalarField& source,tw::Float mul)
{
	// solve div(coeff*grad(phi)) = mul*source

	char *maskNow;
	tw::Int iter,i,j,k,ipass;
	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	tw::Float residual,normResidual,normResidual0=0.0,normResidual1=0.0,normResidual2=0.0;
	tw::Float domega,rp1,rp2;
	std::valarray<tw::Float> D(7);

	normSource = 0.0;
	for (k=1;k<=zDim;k++)
		for (j=1;j<=yDim;j++)
			for (i=1;i<=xDim;i++)
				normSource += fabs(mul*source(i,j,k));
	task->strip[0].AllSum(&normSource,&normSource,sizeof(tw::Float),0);
	normSource += minimumNorm;

	for (iter=0;iter<maxIterations;iter++)
	{
		normResidual = 0.0;
		for (ipass=1;ipass<=2;ipass++)
		{
			maskNow = ipass==1 ? mask1 : mask2;

			for (k=1;k<=zDim;k++)
				for (j=1;j<=yDim;j++)
					for (i=1;i<=xDim;i++)
					{
						if (maskNow[(i-1) + (j-1)*xDim + (k-1)*xDim*yDim]==1)
						{
							FormOperatorStencil(D,i,j,k);
							residual = D[1]*phi(i-1,j,k) + D[2]*phi(i+1,j,k) + D[3]*phi(i,j-1,k) + D[4]*phi(i,j+1,k) + D[5]*phi(i,j,k-1) + D[6]*phi(i,j,k+1) + D[0]*phi(i,j,k) - mul*source(i,j,k);
							phi(i,j,k) -= overrelaxation*residual/D[0];
							normResidual += fabs(residual);
						}
					}

			phi.CopyFromNeighbors();
			phi.ApplyBoundaryCondition();

			// Handle non-zero dirichlet BC's
			// Note: specifying cell wall potential appears unstable

			if (task->n0[3]==MPI_PROC_NULL)
				for (i=phi.N0(1);i<=phi.N1(1);i++)
					for (j=phi.N0(2);j<=phi.N1(2);j++)
						phi(i,j,0) = lbc[i-phi.N0(1)];
			if (task->n1[3]==MPI_PROC_NULL)
				for (i=phi.N0(1);i<=phi.N1(1);i++)
					for (j=phi.N0(2);j<=phi.N1(2);j++)
						phi(i,j,zDim+1) = rbc[i-phi.N0(1)];
		}

		task->strip[0].AllSum(&normResidual,&normResidual,sizeof(tw::Float),0);
		normResidual0 = normResidual1;
		normResidual1 = normResidual2;
		normResidual2 = normResidual;
		if (overrelaxation!=1.0) // If Gauss-Seidel requested, don't try to over-relax
		{
			overrelaxation += overrelaxationChange;
			if (iter>0 && iter%3==0)
			{
				rp1 = (normResidual1-normResidual0)/overrelaxationChange;
				rp2 = (normResidual2-normResidual1)/overrelaxationChange;
				domega = rp2*overrelaxationChange/(rp2-rp1);
				if (domega>0.05) domega=0.05;
				if (domega<-0.05) domega=-0.05;
				overrelaxation += domega - overrelaxationChange;
				if (overrelaxation>1.99) overrelaxation=1.99;
				if (overrelaxation<1.01) overrelaxation=1.01;
				overrelaxationChange *= -1.0;
			}
		}
		if (normResidual <= tolerance*normSource)
			break;
	}
	normResidualAchieved = normResidual;
	iterationsPerformed = iter;
}

#endif

void IterativePoissonSolver::StatusMessage(std::ostream *theStream)
{
	*theStream << "Elliptic Solver Status:" << std::endl;
	*theStream << "   Iterations = " << iterationsPerformed << std::endl;
	*theStream << "   Overrelaxation = " << overrelaxation << std::endl;
	*theStream << "   Norm[source] = " << normSource << std::endl;
	*theStream << "   Norm[residual] = " << normResidualAchieved << std::endl;
}

void IterativePoissonSolver::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	std::string word;
	EllipticSolver::ReadInputFileDirective(inputString,command);
	if (command=="tolerance")
	{
		inputString >> word >> tolerance;
	}
	if (command=="overrelaxation")
	{
		inputString >> word >> overrelaxation;
	}
}

void IterativePoissonSolver::ReadData(std::ifstream& inFile)
{
	EllipticSolver::ReadData(inFile);
	inFile.read((char*)&tolerance,sizeof(tolerance));
	inFile.read((char*)&overrelaxation,sizeof(overrelaxation));
}

void IterativePoissonSolver::WriteData(std::ofstream& outFile)
{
	EllipticSolver::WriteData(outFile);
	outFile.write((char*)&tolerance,sizeof(tolerance));
	outFile.write((char*)&overrelaxation,sizeof(overrelaxation));
}


///////////////////////////////
//                           //
// WAVE STYLE POISSON SOLVER //
//                           //
///////////////////////////////


PoissonSolver::PoissonSolver(const std::string& name,MetricSpace *m,Task *tsk) : EllipticSolver(name,m,tsk)
{
	globalIntegrator = NULL;
	z0 = dirichletWall;
	z1 = neumannWall;
	x0 = x1 = y0 = y1 = periodic;
	typeCode = tw::tool_type::facrPoissonSolver;
}

PoissonSolver::~PoissonSolver()
{
	if (globalIntegrator!=NULL)
		delete globalIntegrator;
}

void PoissonSolver::Initialize()
{
	if (!space->TransversePowersOfTwo())
		throw tw::FatalError("FACR elliptical solver requires all transverse dimensions be powers of two.");

	EllipticSolver::Initialize();

	// SET UP GLOBAL INTEGRATOR

	if (globalIntegrator!=NULL)
		delete globalIntegrator;
	globalIntegrator = new GlobalIntegrator<tw::Float>(&task->strip[3],space->Dim(1)*space->Dim(2),space->Dim(3));
}

void PoissonSolver::Solve(ScalarField& phi,ScalarField& source,tw::Float mul)
{
	tw::Int i,j,k;
	const tw::Int xDim = phi.Dim(1);
	const tw::Int yDim = phi.Dim(2);
	const tw::Int zDim = phi.Dim(3);

	// Transform to frequency space in the transverse direction

	switch (x0)
	{
    	case none:
    	case natural:
		case normalFluxFixed:
        	break;
		case periodic:
			source.TransverseFFT();
			break;
		case neumannWall:
			source.TransverseCosineTransform();
			break;
		case dirichletWall:
			source.TransverseSineTransform();
			break;
		case dirichletCell:
			source.TransverseSineTransform();
			break;
	}

	// Solve on single node

	tw::Float temp,eigenvalue,a,b,c,phi0,phiN1;
	tw::Float dz2 = sqr( gammaBeam * dz(phi) );
	a = 1.0/dz2;
	b = -2.0/dz2;
	c = 1.0/dz2;

	if (z0==natural || z1==natural)
		ZeroModeGhostCellValues(&phi0,&phiN1,source,mul);

	#pragma omp parallel for private(i,j,k,eigenvalue,temp) collapse(2) schedule(static)
	for (j=1;j<=yDim;j++)
		for (i=1;i<=xDim;i++)
		{
			std::valarray<tw::Float> s(zDim),u(zDim),T1(zDim),T2(zDim),T3(zDim);
			if (x0==periodic)
				eigenvalue = phi.CyclicEigenvalue(i,j);
			else
				eigenvalue = phi.Eigenvalue(i,j);
			for (k=1;k<=zDim;k++)
			{
				s[k-1] = mul*source(i,j,k);
				T1[k-1] = a;
				T2[k-1] = b + eigenvalue;
				T3[k-1] = c;
			}
			if (task->n0[3]==MPI_PROC_NULL)
			{
				if (z0==natural)
				{
					if (eigenvalue==0.0)
						s[0] -= T1[0]*phi0;
					else
					{
						temp = QuadraticRoot1(T1[0],T2[0],T3[0]);
						if (temp >= 1.0)
							temp = QuadraticRoot2(T1[0],T2[0],T3[0]);
						T2[0] += T1[0]*temp;
					}
				}
				else
					phi.AdjustTridiagonalForBoundaries(zAxis,lowSide,&T1[0],&T2[0],&T3[0],&s[0],0.0);
			}
			if (task->n1[3]==MPI_PROC_NULL)
			{
				if (z1==natural)
				{
					if (eigenvalue==0.0)
						s[zDim-1] -= T3[zDim-1]*phiN1;
					else
					{
						temp = QuadraticRoot1(T1[zDim-1],T2[zDim-1],T3[zDim-1]);
						if (temp >= 1.0)
							temp = QuadraticRoot2(T1[zDim-1], T2[zDim-1], T3[zDim-1]);
						T2[zDim-1] += T3[zDim-1]*temp;
					}
				}
				else
					phi.AdjustTridiagonalForBoundaries(zAxis,highSide,&T1[zDim-1],&T2[zDim-1],&T3[zDim-1],&s[zDim-1],0.0);
			}
			TriDiagonal<tw::Float,tw::Float>(u,s,T1,T2,T3);
			for (k=1;k<=zDim;k++)
				phi(i,j,k) = u[k-1];

			globalIntegrator->SetMatrix( (j-1)*xDim + (i-1) ,T1,T2,T3);
			globalIntegrator->SetData( (j-1)*xDim + (i-1) ,&phi(i,j,0),phi.Stride(3));
		}

	// Take into account the other nodes

	globalIntegrator->Parallelize();

	// Open boundaries must be done in transformed space

	if (task->n0[3]==MPI_PROC_NULL && z0==natural)
	{
		for (j=1;j<=yDim;j++)
			for (i=1;i<=xDim;i++)
			{
				if (x0==periodic)
					eigenvalue = phi.CyclicEigenvalue(i,j);
				else
					eigenvalue = phi.Eigenvalue(i,j);
				if (eigenvalue==0.0)
					phi(i,j,0) = phi0;
				else
				{
					temp = QuadraticRoot1(a,b+eigenvalue,c);
					if (temp >= 1.0)
						temp = QuadraticRoot2(a,b+eigenvalue,c);
					phi(i,j,0) = temp*phi(i,j,1);
				}
			}
	}

	if (task->n1[3]==MPI_PROC_NULL && z1==natural)
	{
		for (j=1;j<=yDim;j++)
			for (i=1;i<=xDim;i++)
			{
				if (x0==periodic)
					eigenvalue = phi.CyclicEigenvalue(i,j);
				else
					eigenvalue = phi.Eigenvalue(i,j);
				if (eigenvalue==0.0)
					phi(i,j,zDim+1) = phiN1;
				else
				{
					temp = QuadraticRoot1(a,b+eigenvalue,c);
					if (temp >= 1.0)
						temp = QuadraticRoot2(a,b+eigenvalue,c);
					phi(i,j,zDim+1) = temp*phi(i,j,zDim);
				}
			}
	}

	// Transform back to real space

	switch (x0)
	{
		case none:
    	case natural:
		case normalFluxFixed:
			break;
		case periodic:
			phi.InverseTransverseFFT();
			//source.InverseTransverseFFT();
			break;
		case neumannWall:
			phi.InverseTransverseCosineTransform();
			//source.InverseTransverseCosineTransform();
			break;
		case dirichletWall:
			phi.InverseTransverseSineTransform();
			//source.InverseTransverseSineTransform();
			break;
		case dirichletCell:
			phi.InverseTransverseSineTransform();
			//source.InverseTransverseSineTransform();
			break;
	}

	// Global boundary conditions

	phi.ApplyBoundaryCondition();
}



////////////////////////////////
//                            //
// CYLINDRICAL POISSON SOLVER //
//                            //
////////////////////////////////


EigenmodePoissonSolver::EigenmodePoissonSolver(const std::string& name,MetricSpace *m,Task *tsk) : EllipticSolver(name,m,tsk)
{
	globalIntegrator = NULL;
	z0 = dirichletWall;
	z1 = dirichletWall;
	typeCode = tw::tool_type::eigenmodePoissonSolver;
}

EigenmodePoissonSolver::~EigenmodePoissonSolver()
{
	if (globalIntegrator!=NULL)
		delete globalIntegrator;
}

void EigenmodePoissonSolver::Initialize()
{
	EllipticSolver::Initialize();

	const tw::Int rDim = space->Dim(1);
	const tw::Int zDim = space->Dim(3);

	ComputeTransformMatrices(eigenvalue,hankel,inverseHankel,space,task);

	// SET UP GLOBAL INTEGRATOR

	if (globalIntegrator!=NULL)
		delete globalIntegrator;
	globalIntegrator = new GlobalIntegrator<tw::Float>(&task->strip[3],rDim,zDim);
}

void EigenmodePoissonSolver::TransformBoundaryValues()
{
	// Presuming owner has already set lbc/rbc in real space, express lbc_t/rbc_t in eigenspace
	const tw::Int xDim = space->Dim(1);
	std::valarray<tw::Float> globalData(task->globalCells[1]+2);

	// LEFT SIDE
	task->strip[1].Gather(&lbc[space->Layers(1)],&globalData[1],xDim*sizeof(tw::Float),0);
	if (task->strip[1].Get_rank()==0)
		Transform(&globalData[1],task->globalCells[1],1,hankel);
	task->strip[1].Scatter(&lbc[space->Layers(1)],&globalData[1],xDim*sizeof(tw::Float),0);

	// RIGHT SIDE
	task->strip[1].Gather(&rbc[space->Layers(1)],&globalData[1],xDim*sizeof(tw::Float),0);
	if (task->strip[1].Get_rank()==0)
		Transform(&globalData[1],task->globalCells[1],1,hankel);
	task->strip[1].Scatter(&rbc[space->Layers(1)],&globalData[1],xDim*sizeof(tw::Float),0);
}

void EigenmodePoissonSolver::Solve(ScalarField& phi,ScalarField& source,tw::Float mul)
{
	// lbc_t and rbc_t must have been set up via TransformBoundaryValues externally
	tw::Int i,k;
	tw::Int rDim = space->Dim(1);
	tw::Int zDim = space->Dim(3);
	tw::Float temp,dz,dz1,dz2,phi0,phiN1;
	tw::Float T1_lbc,T2_lbc,T3_lbc,T1_rbc,T2_rbc,T3_rbc;
	std::valarray<tw::Float> localEig(rDim+2);

	Field T;
	for (i=1;i<=rDim;i++)
		localEig[i] = eigenvalue[i-1+task->cornerCell[1]-1];

	// Matrix elements for open boundary condition

	T1_lbc = T3_lbc = 1.0/sqr(space->dX(0,3));
	T2_lbc = -2.0/sqr(space->dX(0,3)); // add the eigenvalue later
	T1_rbc = T3_rbc = 1.0/sqr(space->dX(zDim+1,3));
	T2_rbc = -2.0/sqr(space->dX(zDim+1,3)); // add the eigenvalue later

	// Perform hankel transform
	source.Transpose(xAxis,zAxis,&T,1);
	for (i=T.N0(3);i<=T.N1(3);i++)
		Transform(&T(1,1,i,0),task->globalCells[1],T.Stride(1),hankel);
	source.Transpose(xAxis,zAxis,&T,-1);

	// Solve on single node

	if (z0==natural || z1==natural)
		ZeroModeGhostCellValues(&phi0,&phiN1,source,mul);

	#pragma omp parallel for private(i,k,temp) schedule(static)
	for (i=1;i<=rDim;i++)
	{
		std::valarray<tw::Float> s(zDim),u(zDim),T1(zDim),T2(zDim),T3(zDim);
		for (k=1;k<=zDim;k++)
		{
			s[k-1] = mul*source(i,0,k);
			dz = space->dX(k,3);
			dz1 = 0.5*(space->dX(k-1,3)+dz);
			dz2 = 0.5*(space->dX(k+1,3)+dz);
			T1[k-1] = 1.0/(dz*dz1);
			T2[k-1] = -(1.0/(dz*dz1) + 1.0/(dz*dz2)) + localEig[i];
			T3[k-1] = 1.0/(dz*dz2);
		}
		if (task->n0[3]==MPI_PROC_NULL)
		{
			if (z0==natural)
			{
				if (localEig[i]==0.0)
					s[0] -= T1[0]*phi0;
				else
				{
					temp = QuadraticRoot1(T1_lbc,T2_lbc+localEig[i],T3_lbc);
					if (temp >= 1.0)
						temp = QuadraticRoot2(T1_lbc,T2_lbc+localEig[i],T3_lbc);
					T2[0] += T1[0]*temp;
				}
			}
			else
				phi.AdjustTridiagonalForBoundaries(zAxis,lowSide,&T1[0],&T2[0],&T3[0],&s[0],lbc_t[i-phi.N0(1)]);
		}
		if (task->n1[3]==MPI_PROC_NULL)
		{
			if (z1==natural)
			{
				if (localEig[i]==0.0)
					s[zDim-1] -= T3[zDim-1]*phiN1;
				else
				{
					temp = QuadraticRoot1(T1_rbc,T2_rbc+localEig[i],T3_rbc);
					if (temp >= 1.0)
						temp = QuadraticRoot2(T1_rbc,T2_rbc+localEig[i],T3_rbc);
					T2[zDim-1] += T3[zDim-1]*temp;
				}
			}
			else
				phi.AdjustTridiagonalForBoundaries(zAxis,highSide,&T1[zDim-1],&T2[zDim-1],&T3[zDim-1],&s[zDim-1],rbc_t[i-phi.N0(1)]);
		}

		TriDiagonal<tw::Float,tw::Float>(u,s,T1,T2,T3);
		for (k=1;k<=zDim;k++)
			phi(i,0,k) = u[k-1];

		globalIntegrator->SetMatrix(i-1,T1,T2,T3);
		globalIntegrator->SetData(i-1,&phi(i,0,0),phi.Stride(3));
	}

	// Take into account the other nodes

	globalIntegrator->Parallelize();

	// Open boundaries have to be done in transformed space

	if (task->n0[3]==MPI_PROC_NULL && z0==natural)
	{
		for (i=1;i<=rDim;i++)
		{
			if (localEig[i]==0.0)
				phi(i,0,0) = phi0;
			else
			{
				temp = QuadraticRoot1(T1_lbc,T2_lbc+localEig[i],T3_lbc);
				if (temp >= 1.0)
					temp = QuadraticRoot2(T1_lbc,T2_lbc+localEig[i],T3_lbc);
				phi(i,0,0) = temp*phi(i,0,1);
			}
		}
	}

	if (task->n1[3]==MPI_PROC_NULL && z1==natural)
	{
		for (i=1;i<=rDim;i++)
		{
			if (localEig[i]==0.0)
				phi(i,0,zDim+1) = phiN1;
			else
			{
				temp = QuadraticRoot1(T1_rbc,T2_rbc+localEig[i],T3_rbc);
				if (temp >= 1.0)
					temp = QuadraticRoot2(T1_rbc,T2_rbc+localEig[i],T3_rbc);
				phi(i,0,zDim+1) = temp*phi(i,0,zDim);
			}
		}
	}

	// back to real space

	phi.Transpose(xAxis,zAxis,&T,1);
	for (i=T.N0(3);i<=T.N1(3);i++)
		Transform(&T(1,1,i,0),task->globalCells[1],T.Stride(1),inverseHankel);
	phi.Transpose(xAxis,zAxis,&T,-1);

// 	source.Transpose(xAxis,zAxis,&T,1);
// 	for (i=ZN0(T);i<=ZN1(T);i++)
// 		Transform(&T(1,1,i,0),task->globalCells[1],XStride(T),inverseHankel);
// 	source.Transpose(xAxis,zAxis,&T,-1);

	// Global boundary conditions

	phi.ApplyBoundaryCondition();
	for (i=phi.N0(1);i<=phi.N1(1);i++)
	{
		phi(i,0,0) = lbc[i-phi.N0(1)];
		phi(i,0,zDim+1) = rbc[i-phi.N0(1)];
	}
}
