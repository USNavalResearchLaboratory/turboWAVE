#include "definitions.h"
#include "tasks.h"
#include "ctools.h"
#include "3dmath.h"
#include "metricSpace.h"
#include "3dfields.h"
#include "fft.h"


/////////////////////////
//                     //
// BOUNDARY CONDITIONS //
//                     //
/////////////////////////


BoundaryCondition::BoundaryCondition()
{
	// default to none (identity matrices)
	Reset();
	// default to low side
	sgn = 1;
}

void BoundaryCondition::Reset()
{
	// Set all variables to the identity operation
	for (tw::Int i=0;i<4;i++)
	{
		for (tw::Int j=0;j<4;j++)
		{
			fold[i][j] = 0.0;
			force[i][j] = 0.0;
		}
		coeff[i] = 0.0;
		fold[i][i] = 1.0;
		force[i][i] = 1.0;
	}
}

void BoundaryCondition::Set(boundarySpec theBoundaryCondition,sideSpec whichSide)
{
	// In this routine, indexing is offset from DiscreteSpace indexing convention.
	// Namely, 0 is the outer ghost cell layer, 1 is the inner, 2 is the edge cell, etc.
	// When BC is applied, negative strides are used to produce mirror image indexing on high side

	Reset();
	switch (whichSide)
	{
		case lowSide:
			sgn = 1;
			break;
		case highSide:
			sgn = -1;
			break;
	}
	switch (theBoundaryCondition)
	{
		case normalFluxFixed:
			// quantities known on cell walls that are fixed on the wall
			if (whichSide==lowSide)
			{
				fold[3][1] = -1.0;
				force[0][0] = 0.0;
				force[1][1] = 0.0;
				force[2][2] = 0.0;
				coeff[2] = 1.0;
			}
			else
			{
				fold[2][0] = -1.0;
				force[0][0] = 0.0;
				force[1][1] = 0.0;
				coeff[1] = 1.0;
			}
			break;
		case neumannWall:
			// quantities known at cell centers that are continued outside domain
			// strip_1 = strip_2 - BC <--> strip_2 - strip_1 = BC
			force[0][0] = 0.0;
			force[1][1] = 0.0;
			force[0][3] = 1.0;
			force[1][2] = 1.0;
			coeff[0] = -1.0;
			coeff[1] = -1.0;
			break;
		case dirichletWall:
			// quantities known at cell centers which are fixed on cell walls
			// (e.g., a no-slip boundary condition for tangential fluid velocity)
			// strip_1 = 2*BC - strip_2 <--> (strip_1+strip_2)/2 = BC
			force[0][0] = 0.0;
			force[1][1] = 0.0;
			force[0][3] = -1.0;
			force[1][2] = -1.0;
			coeff[0] = 2.0;
			coeff[1] = 2.0;
			break;
		case dirichletCell:
			// quantities known at cell centers which are fixed outside domain
			fold[2][1] = 1.0;
			fold[3][0] = 1.0;
			force[0][0] = 0.0;
			force[1][1] = 0.0;
			coeff[0] = 1.0;
			coeff[1] = 1.0;
			break;
		case none:
		case periodic:
		case natural:
			// leave as the identity
			break;
	}
}



/////////////
//         //
// ELEMENT //
//         //
/////////////


Element::Element()
{
	low = 0;
	high = 0;
}

Element::Element(tw::Int l,tw::Int h)
{
	low = l;
	high = h;
}

Element::Element(tw::Int i)
{
	low = i;
	high = i;
}

tw::Int Element::Components() const
{
	return high - low + 1;
}

Element Union(const Element& e1,const Element& e2)
{
	if (e1.low < e2.low)
		return Element(e1.low,e2.high);
	else
		return Element(e2.low,e1.high);
}


///////////////////
//               //
//  FIELD CLASS  //
//               //
///////////////////


Field::Field()
{
	for (int i=0;i<4;i++)
		num[i] = 0;
	totalCells = 0;
	boundaryCells = 0;
	ghostCells = 0;
	bufferState = 0;
	packedAxis = 3;
}

Field::~Field()
{
	#ifdef USE_OPENCL
	if (bufferState==1)
	{
		clReleaseMemObject(computeBuffer);
		clReleaseMemObject(boundaryBuffer);
		clReleaseMemObject(boundaryMapBuffer);
		clReleaseMemObject(ghostBuffer);
		clReleaseMemObject(ghostMapBuffer);
	}
	#endif
}

void Field::Initialize(tw::Int components,const DiscreteSpace& ds,Task *task,const axisSpec& axis)
{
	DiscreteSpace::operator=(ds);
	this->task = task;
	packedAxis = naxis(axis);
	totalCells = num[1]*num[2]*num[3];
	num[0] = components;
	bc0.resize(4,components);
	bc1.resize(4,components);
	for (tw::Int i=0;i<4;i++)
		for (tw::Int c=0;c<components;c++)
		{
			bc0(i,c).Set(none,lowSide);
			bc1(i,c).Set(none,highSide);
		}
	if (packedAxis==0)
	{
		stride[0] = 1;
		//stride[1] = num[0];
		//stride[2] = num[0]*num[1];
		//stride[3] = num[0]*num[1]*num[2];
		stride[1] = num[0]*num[2]*num[3];
		stride[2] = num[0]*num[3];
		stride[3] = num[0];
	}
	if (packedAxis==1)
	{
		stride[0] = num[1]*num[2]*num[3];
		stride[1] = 1;
		stride[2] = num[1];
		stride[3] = num[1]*num[2];
	}
	if (packedAxis==2)
	{
		stride[0] = num[1]*num[2]*num[3];
		stride[1] = num[2];
		stride[2] = 1;
		stride[3] = num[1]*num[2];
	}
	if (packedAxis==3)
	{
		stride[0] = num[1]*num[2]*num[3];
		stride[1] = num[2]*num[3];
		stride[2] = num[3];
		stride[3] = 1;
	}
	stride[1] *= dim[1]==1 ? 0 : 1;
	stride[2] *= dim[2]==1 ? 0 : 1;
	stride[3] *= dim[3]==1 ? 0 : 1;
	// note : boundary and ghost cells have redundant information
	// (this simplifies the storage scheme, but demands operations in compute buffer be atomic)
	// boundary cells means ghost cells + as many adjacent interior cells as there are ghost cells
	tw::Int boundaryCellsVec[4],ghostCellsVec[4];
	boundaryCellsVec[1] = dim[1]>1 ? 4*layers[1]*num[2]*num[3] : 0;
	boundaryCellsVec[2] = dim[2]>1 ? 4*layers[2]*num[1]*num[3] : 0;
	boundaryCellsVec[3] = dim[3]>1 ? 4*layers[3]*num[1]*num[2] : 0;
	ghostCellsVec[1] = dim[1]>1 ? 2*layers[1]*num[2]*num[3] : 0;
	ghostCellsVec[2] = dim[2]>1 ? 2*layers[2]*num[1]*num[3] : 0;
	ghostCellsVec[3] = dim[3]>1 ? 2*layers[3]*num[1]*num[2] : 0;
	boundaryCells = boundaryCellsVec[1] + boundaryCellsVec[2] + boundaryCellsVec[3];
	ghostCells = ghostCellsVec[1] + ghostCellsVec[2] + ghostCellsVec[3];
	array.resize(totalCells*num[0]);
	boundaryData.resize(boundaryCells*num[0]);
	ghostData.resize(ghostCells*num[0]);
	boundaryDataIndexMap.resize(boundaryCells*num[0]);
	ghostDataIndexMap.resize(ghostCells*num[0]);

	// Set up the index maps going from coordinate space to memory space.
	// Ordering is not important as long as we gather all the right cells.

	tw::Int bOffset=0,gOffset=0;
	for (tw::Int c=0;c<num[0];c++)
		for (tw::Int ax=1;ax<=3;ax++)
		{
			if (dim[ax]>1)
				for (auto s : StripRange(*this,ax,strongbool::yes))
				{
					for (tw::Int i=0;i<2*layers[ax];i++)
					{
						boundaryDataIndexMap[bOffset+i] = s.Index(lb[ax]+i,c,stride);
						boundaryDataIndexMap[bOffset+2*layers[ax]+i] = s.Index(ub[ax]-2*layers[ax]+i+1,c,stride);
					}
					for (tw::Int i=0;i<layers[ax];i++)
					{
						ghostDataIndexMap[gOffset+i] = s.Index(lb[ax]+i,c,stride);
						ghostDataIndexMap[gOffset+layers[ax]+i] = s.Index(ub[ax]-layers[ax]+i+1,c,stride);
					}
					bOffset += 4*layers[ax];
					gOffset += 2*layers[ax];
				}
		}
}

void Field::MultiplyCellVolume(const MetricSpace& m)
{
	for (auto cell : EntireCellRange(*this))
		for (tw::Int c=0;c<num[0];c++)
			(*this)(cell,c) *= m.dS(cell,0);
}

void Field::DivideCellVolume(const MetricSpace& m)
{
	for (auto cell : EntireCellRange(*this))
		for (tw::Int c=0;c<num[0];c++)
			(*this)(cell,c) /= m.dS(cell,0);
}

void Field::Shift(const Element& e,const tw::strip& s,tw::Int cells,const tw::Float* incoming)
{
	// Propagate field pattern to the right if <cells> positive, to the left if negative
	// argument <incoming> is value to inject into left or right ghost cell

	tw::Int i,c,ax=s.Axis();

	if (cells>0)
	{
		for (i=N1(ax);i>=N0(ax)+cells;i--)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=N0(ax);i<N0(ax)+cells;i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = incoming[c-e.low];
	}
	if (cells<0)
	{
		for (i=N0(ax);i<=N1(ax)+cells;i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=N1(ax)+1+cells;i<=N1(ax);i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = incoming[c-e.low];
	}
}

void Field::Shift(const Element& e,const tw::strip& s,tw::Int cells,const tw::Float& incoming)
{
	// Propagate field pattern to the right if <cells> positive, to the left if negative
	// argument <incoming> is value to inject into left or right ghost cell

	tw::Int i,c,ax=s.Axis();

	if (cells>0)
	{
		for (i=N1(ax);i>=N0(ax)+cells;i--)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=N0(ax);i<N0(ax)+cells;i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = incoming;
	}
	if (cells<0)
	{
		for (i=N0(ax);i<=N1(ax)+cells;i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = (*this)(s,i-cells,c);
		for (i=N1(ax)+1+cells;i<=N1(ax);i++)
			for (c=e.low;c<=e.high;c++)
				(*this)(s,i,c) = incoming;
	}
}

//////////////////////////////////
//                              //
//  Global Boundary Conditions  //
//                              //
//////////////////////////////////


void Field::SetBoundaryConditions(const Element& e,const axisSpec& axis,boundarySpec low,boundarySpec high)
{
	tw::Int i;
	tw::Int ax = naxis(axis);

	for (i=e.low;i<=e.high;i++)
	{
		if (task->n0[ax]==MPI_PROC_NULL && dim[ax]>1)
			bc0(ax,i).Set(low,lowSide);
		if (task->n1[ax]==MPI_PROC_NULL && dim[ax]>1)
			bc1(ax,i).Set(high,highSide);
	}
}

void Field::ZeroGhostCells(const Element& e)
{
	tw::Int ax,c,s;
	for (c=e.low;c<=e.high;c++)
		for (ax=1;ax<=3;ax++)
			for (auto strip : StripRange(*this,ax,strongbool::yes))
				for (s=0;s<layers[ax];s++)
				{
					(*this)(strip,lb[ax]+s,c) = 0.0;
					(*this)(strip,dim[ax]+s+1,c) = 0.0;
				}
}

void Field::ApplyFoldingCondition(const Element& e)
{
	tw::Int ax,c;
	for (c=e.low;c<=e.high;c++)
		for (ax=1;ax<=3;ax++)
			if (num[ax]>1)
				for (auto strip : StripRange(*this,ax,strongbool::yes))
				{
					bc0(ax,c).FoldingOperation(&(*this)(strip,lb[ax],c),stride[ax]);
					bc1(ax,c).FoldingOperation(&(*this)(strip,ub[ax],c),stride[ax]);
				}
}

void Field::ApplyBoundaryCondition(const Element& e)
{
	tw::Int ax,c;
	for (c=e.low;c<=e.high;c++)
		for (ax=1;ax<=3;ax++)
			if (num[ax]>1)
				for (auto strip : StripRange(*this,ax,strongbool::yes))
				{
					bc0(ax,c).ForcingOperation(&(*this)(strip,lb[ax],c),stride[ax],0.0);
					bc1(ax,c).ForcingOperation(&(*this)(strip,ub[ax],c),stride[ax],0.0);
				}
}

void Field::AdjustTridiagonalForBoundaries(const Element& e,const axisSpec& axis,const sideSpec& side,tw::Float *T1,tw::Float *T2,tw::Float *T3,tw::Float *source,tw::Float *val)
{
	// General rule is ghostCell_i = force_ij * cell_j + coeff_i * BC
	// N.b. the indexing starts with outermost ghost cell = 0 and works inward, regardless of sideSpec
	// 2 layers is assumed in the indexing, but only 1 enters into the calculation.
	// Assumes only force_12, force_13, and coeff_1 come into the tridiagonal equation.
	tw::Int i;
	tw::Int a = naxis(axis);
	if (side==lowSide)
	{
		for (i=e.low;i<=e.high;i++)
		{
			T2[i-e.low] += T1[i-e.low] * bc0(a,i).force[1][2];
			T3[i-e.low] += T1[i-e.low] * bc0(a,i).force[1][3];
			source[i-e.low] -= T1[i-e.low] * bc0(a,i).coeff[1] * val[i-e.low];
		}
	}
	if (side==highSide)
	{
		for (i=e.low;i<=e.high;i++)
		{
			T2[i-e.low] += T3[i-e.low] * bc1(a,i).force[1][2];
			T1[i-e.low] += T3[i-e.low] * bc1(a,i).force[1][3];
			source[i-e.low] -= T3[i-e.low] * bc1(a,i).coeff[1] * val[i-e.low];
		}
	}
}

void Field::BoundaryDataToField()
{
	for (tw::Int i=0;i<boundaryData.size();i++)
		array[boundaryDataIndexMap[i]] = boundaryData[i];
}

void Field::FieldToBoundaryData()
{
	for (tw::Int i=0;i<boundaryData.size();i++)
		boundaryData[i] = array[boundaryDataIndexMap[i]];
}

void Field::FieldToGhostData()
{
	for (tw::Int i=0;i<ghostData.size();i++)
		ghostData[i] = array[ghostDataIndexMap[i]];
}

/////////////////////////////
//                         //
//  GPU Related Routines   //
//                         //
/////////////////////////////

#ifdef USE_OPENCL

void Field::InitializeComputeBuffer()
{
	cl_int err;
	if (bufferState==1)
	{
		clReleaseMemObject(computeBuffer);
		clReleaseMemObject(boundaryBuffer);
		clReleaseMemObject(boundaryMapBuffer);
		clReleaseMemObject(ghostBuffer);
		clReleaseMemObject(ghostMapBuffer);
	}
	computeBuffer = clCreateBuffer(task->context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,TotalBytes(),&array[0],&err);
	boundaryBuffer = clCreateBuffer(task->context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,BoundaryBytes(),&boundaryData[0],&err);
	boundaryMapBuffer = clCreateBuffer(task->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,BoundaryMapBytes(),&boundaryDataIndexMap[0],&err);
	ghostBuffer = clCreateBuffer(task->context,CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,GhostBytes(),&ghostData[0],&err);
	ghostMapBuffer = clCreateBuffer(task->context,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,GhostMapBytes(),&ghostDataIndexMap[0],&err);
	bufferState = 1;

	if (err)
		throw tw::FatalError("Device buffer failure in Field.");
}

void Field::SendToComputeBuffer()
{
	clEnqueueWriteBuffer(task->commandQueue,computeBuffer,CL_TRUE,0,TotalBytes(),&array[0],  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void Field::ReceiveFromComputeBuffer()
{
	clEnqueueReadBuffer(task->commandQueue,computeBuffer,CL_TRUE,0,TotalBytes(),&array[0],  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void Field::SendBoundaryCellsToComputeBuffer()
{
	FieldToBoundaryData();
	size_t elements = boundaryData.size();
	clEnqueueWriteBuffer(task->commandQueue,boundaryBuffer,CL_TRUE,0,BoundaryBytes(),&boundaryData[0]  ,0,NULL,NULL);
	clFinish(task->commandQueue);
	clSetKernelArg(task->k_boundaryToField,0,sizeof(computeBuffer),&computeBuffer);
	clSetKernelArg(task->k_boundaryToField,1,sizeof(boundaryBuffer),&boundaryBuffer);
	clSetKernelArg(task->k_boundaryToField,2,sizeof(boundaryMapBuffer),&boundaryMapBuffer);
	clEnqueueNDRangeKernel(task->commandQueue,task->k_boundaryToField,1,NULL,&elements,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void Field::ReceiveBoundaryCellsFromComputeBuffer()
{
	size_t elements = boundaryData.size();
	clSetKernelArg(task->k_fieldToBoundary,0,sizeof(computeBuffer),&computeBuffer);
	clSetKernelArg(task->k_fieldToBoundary,1,sizeof(boundaryBuffer),&boundaryBuffer);
	clSetKernelArg(task->k_fieldToBoundary,2,sizeof(boundaryMapBuffer),&boundaryMapBuffer);
	clEnqueueNDRangeKernel(task->commandQueue,task->k_fieldToBoundary,1,NULL,&elements,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);
	clEnqueueReadBuffer(task->commandQueue,boundaryBuffer,CL_TRUE,0,BoundaryBytes(),&boundaryData[0]  ,0,NULL,NULL);
	clFinish(task->commandQueue);
	BoundaryDataToField();
}

tw::Float Field::CellValueInComputeBuffer(tw::Int i,tw::Int j,tw::Int k,tw::Int c)
{
	tw::Float ans;
	tw::Int idx = c*stride[0] + (i-lb[1])*stride[1] + (j-lb[2])*stride[2] + (k-lb[3])*stride[3];
	clEnqueueReadBuffer(task->commandQueue,computeBuffer,CL_TRUE,sizeof(tw::Float)*idx,sizeof(tw::Float),&ans,  0,NULL,NULL);
	clFinish(task->commandQueue);
	return ans;
}

void Field::UpdateGhostCellsInComputeBuffer(const Element& e)
{
	ReceiveBoundaryCellsFromComputeBuffer();
	CopyFromNeighbors(e);
	ApplyBoundaryCondition(e);
	SendBoundaryCellsToComputeBuffer();
}

void Field::SendGhostCellsToComputeBuffer()
{
	FieldToGhostData();
	size_t elements = ghostData.size();
	clEnqueueWriteBuffer(task->commandQueue,ghostBuffer,CL_TRUE,0,GhostBytes(),&ghostData[0]  ,0,NULL,NULL);
	clFinish(task->commandQueue);
	clSetKernelArg(task->k_ghostToField,0,sizeof(computeBuffer),&computeBuffer);
	clSetKernelArg(task->k_ghostToField,1,sizeof(ghostBuffer),&ghostBuffer);
	clSetKernelArg(task->k_ghostToField,2,sizeof(ghostMapBuffer),&ghostMapBuffer);
	clEnqueueNDRangeKernel(task->commandQueue,task->k_ghostToField,1,NULL,&elements,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void Field::ZeroGhostCellsInComputeBuffer()
{
	size_t cells = ghostCells * num[0];
	clSetKernelArg(task->k_zeroGhostCells,0,sizeof(computeBuffer),&computeBuffer);
	clSetKernelArg(task->k_zeroGhostCells,1,sizeof(ghostMapBuffer),&ghostMapBuffer);
	clEnqueueNDRangeKernel(task->commandQueue,task->k_zeroGhostCells,1,NULL,&cells,NULL,  0,NULL,NULL);
	clFinish(task->commandQueue);
}

void SwapComputeBuffers(Field& f1,Field& f2)
{
	clSetKernelArg(f1.task->k_swapBuffers,0,sizeof(cl_mem),&f1.computeBuffer);
	clSetKernelArg(f1.task->k_swapBuffers,1,sizeof(cl_mem),&f2.computeBuffer);
	f1.ElementUpdateProtocol(f1.task->k_swapBuffers,f1.task->commandQueue);
}

void CopyComputeBuffer(Field& dst,Field& src)
{
	clEnqueueCopyBuffer(dst.task->commandQueue,src.computeBuffer,dst.computeBuffer,0,0,dst.TotalBytes(),  0,NULL,NULL);
	clFinish(dst.task->commandQueue);
}

void CopyComplexComputeBufferMod2(Field& dst,Field& src)
{
	// dst is a real scalar buffer, src is a c-packed complex scalar buffer.
	// ElementUpdateProtocol must be called on the real buffer, not the complex buffer.
	clSetKernelArg(dst.task->k_complexMod2,0,sizeof(cl_mem),&dst.computeBuffer);
	clSetKernelArg(dst.task->k_complexMod2,1,sizeof(cl_mem),&src.computeBuffer);
	dst.ElementUpdateProtocol(dst.task->k_complexMod2,dst.task->commandQueue);
}

void Field::DestructiveComplexMod2ComputeBuffer()
{
	// The field is assumed to be a c-packed complex scalar buffer.
	clSetKernelArg(task->k_destructiveComplexMod2,0,sizeof(cl_mem),&computeBuffer);
	CellUpdateProtocol(task->k_destructiveComplexMod2,task->commandQueue);
}

void Field::FillComputeBufferVec4(const Element& e,tw::vec4& A)
{
	// Overwrite elements in the set e
	// Preserve elements not in the set e
	// e must be in the range of a 4-vector
	tw::Int i;
	// Complex arguments are a device to encode which elements to preserve
	tw::Complex Ac[4];
	for (i=0;i<4;i++)
		Ac[i] = tw::Complex(0.0,1.0);
	for (i=e.low;i<=e.high;i++)
		Ac[i] = tw::Complex(A[i],0.0);
	clSetKernelArg(task->k_fillVec4Field,0,sizeof(cl_mem),&computeBuffer);
	clSetKernelArg(task->k_fillVec4Field,1,sizeof(tw::Complex),&Ac[0]);
	clSetKernelArg(task->k_fillVec4Field,2,sizeof(tw::Complex),&Ac[1]);
	clSetKernelArg(task->k_fillVec4Field,3,sizeof(tw::Complex),&Ac[2]);
	clSetKernelArg(task->k_fillVec4Field,4,sizeof(tw::Complex),&Ac[3]);
	CellUpdateProtocol(task->k_fillVec4Field,task->commandQueue);
}

tw::Float Field::DestructiveSumComputeBuffer()
{
	// add up all the components in all the cells (ghost cells included)
	size_t elements = totalCells*num[0];
	clSetKernelArg(task->k_destructiveSum,0,sizeof(cl_mem),&computeBuffer);
	while (elements>1)
	{
		clEnqueueNDRangeKernel(task->commandQueue,task->k_destructiveSum,1,NULL,&elements,NULL,0,NULL,NULL);
		clFinish(task->commandQueue);
		elements /= 2;
	}
	return CellValueInComputeBuffer(0,0,0,0);
}

tw::Float Field::DestructiveNorm1ComputeBuffer()
{
	// add up absolute value of all the components in all the cells (ghost cells included)
	size_t elements = totalCells*num[0];
	clSetKernelArg(task->k_destructiveNorm1,0,sizeof(cl_mem),&computeBuffer);
	while (elements>1)
	{
		clEnqueueNDRangeKernel(task->commandQueue,task->k_destructiveNorm1,1,NULL,&elements,NULL,0,NULL,NULL);
		clFinish(task->commandQueue);
		elements /= 2;
	}
	return CellValueInComputeBuffer(0,0,0,0);
}

void Field::WeightComputeBufferByVolume(MetricSpace& ms,tw::Float inv)
{
	clSetKernelArg(task->k_weightByVolume,0,sizeof(cl_mem),&computeBuffer);
	clSetKernelArg(task->k_weightByVolume,1,sizeof(cl_mem),&ms.metricsBuffer);
	clSetKernelArg(task->k_weightByVolume,3,sizeof(tw::Float),&inv);
	for (tw::Int c=0;c<num[0];c++)
	{
		clSetKernelArg(task->k_weightByVolume,2,sizeof(tw::Int),&c);
		LocalUpdateProtocol(task->k_weightByVolume,task->commandQueue);
	}
}

void Field::MADDComputeBuffer(tw::Float m,tw::Float a)
{
	clSetKernelArg(task->k_MADD,0,sizeof(cl_mem),&computeBuffer);
	clSetKernelArg(task->k_MADD,1,sizeof(tw::Float),&m);
	clSetKernelArg(task->k_MADD,2,sizeof(tw::Float),&a);
	ElementUpdateProtocol(task->k_MADD,task->commandQueue);
}

#endif

///////////////////////////////////////
//                                   //
//  Boundaries with Message Passing  //
//                                   //
///////////////////////////////////////


void Field::StripCopyProtocol(tw::Int axis,tw::Int shift,Slice<tw::Float> *planeIn,Slice<tw::Float> *planeOut,bool add)
{
	tw::Int src,dst;
	bool odd;

	task->strip[axis].Shift(1,shift,&src,&dst);
	odd = task->strip[axis].Get_rank() % 2;

	if (dst!=MPI_PROC_NULL)
		LoadDataIntoImage<tw::Float>(planeOut);

	if (odd)
	{
		task->strip[axis].Recv(planeIn->Buffer(),planeIn->BufferSize(),src);
		task->strip[axis].Send(planeOut->Buffer(),planeOut->BufferSize(),dst);
	}
	else
	{
		task->strip[axis].Send(planeOut->Buffer(),planeOut->BufferSize(),dst);
		task->strip[axis].Recv(planeIn->Buffer(),planeIn->BufferSize(),src);
	}

	if (src!=MPI_PROC_NULL)
	{
		if (add)
			AddDataFromImage<tw::Float>(planeIn);
		else
			SaveDataFromImage<tw::Float>(planeIn);
	}
}

void Field::DownwardCopy(const axisSpec& axis,const Element& e,tw::Int cells)
{
	Slice<tw::Float> *planeIn,*planeOut;

	switch (axis)
	{
		case tAxis:
			return;
			break;
		case xAxis:
			if (dim[1]==1) return;
			planeIn = new Slice<tw::Float>(e,dim[1]+1,dim[1]+cells, lb[2],ub[2], lb[3],ub[3]);
			planeOut = new Slice<tw::Float>(e,1,cells, lb[2],ub[2], lb[3],ub[3]);
			break;
		case yAxis:
			if (dim[2]==1) return;
			planeIn = new Slice<tw::Float>(e,lb[1],ub[1], dim[2]+1,dim[2]+cells, lb[3],ub[3]);
			planeOut = new Slice<tw::Float>(e,lb[1],ub[1], 1,cells, lb[3],ub[3]);
			break;
		case zAxis:
			if (dim[3]==1) return;
			planeIn = new Slice<tw::Float>(e,lb[1],ub[1], lb[2],ub[2], dim[3]+1,dim[3]+cells);
			planeOut = new Slice<tw::Float>(e,lb[1],ub[1], lb[2],ub[2], 1,cells);
			break;
	}

	StripCopyProtocol(naxis(axis),-1,planeIn,planeOut,false);

	delete planeOut;
	delete planeIn;
}

void Field::UpwardCopy(const axisSpec& axis,const Element& e,tw::Int cells)
{
	Slice<tw::Float> *planeIn,*planeOut;

	switch (axis)
	{
		case tAxis:
			return;
			break;
		case xAxis:
			if (dim[1]==1) return;
			planeIn = new Slice<tw::Float>(e,1-cells,0, lb[2],ub[2], lb[3],ub[3]);
			planeOut = new Slice<tw::Float>(e,dim[1]-cells+1,dim[1], lb[2],ub[2], lb[3],ub[3]);
			break;
		case yAxis:
			if (dim[2]==1) return;
			planeIn = new Slice<tw::Float>(e,lb[1],ub[1], 1-cells,0, lb[3],ub[3]);
			planeOut = new Slice<tw::Float>(e,lb[1],ub[1], dim[2]-cells+1,dim[2], lb[3],ub[3]);
			break;
		case zAxis:
			if (dim[3]==1) return;
			planeIn = new Slice<tw::Float>(e,lb[1],ub[1], lb[2],ub[2], 1-cells,0);
			planeOut = new Slice<tw::Float>(e,lb[1],ub[1], lb[2],ub[2], dim[3]-cells+1,dim[3]);
			break;
	}

	StripCopyProtocol(naxis(axis),1,planeIn,planeOut,false);

	delete planeOut;
	delete planeIn;
}

void Field::DownwardDeposit(const axisSpec& axis,const Element& e,tw::Int cells)
{
	// deposits work by sending both ghost cells and edge cells one way and adding.
	// Before accepting the return message these cells are zeroed.
	// The return message then effectively overwrites them with the correct data.
	Slice<tw::Float> *planeIn,*planeOut;

	switch (axis)
	{
		case tAxis:
			break;
		case xAxis:
			if (dim[1]==1) return;
			planeIn = new Slice<tw::Float>(e,dim[1]-cells+1,dim[1]+cells, lb[2],ub[2], lb[3],ub[3]);
			planeOut = new Slice<tw::Float>(e,1-cells,1+cells-1, lb[2],ub[2], lb[3],ub[3]);
			break;
		case yAxis:
			if (dim[2]==1) return;
			planeIn = new Slice<tw::Float>(e,lb[1],ub[1], dim[2]-cells+1,dim[2]+cells, lb[3],ub[3]);
			planeOut = new Slice<tw::Float>(e,lb[1],ub[1], 1-cells,1+cells-1, lb[3],ub[3]);
			break;
		case zAxis:
			if (dim[3]==1) return;
			planeIn = new Slice<tw::Float>(e,lb[1],ub[1], lb[2],ub[2], dim[3]-cells+1,dim[3]+cells);
			planeOut = new Slice<tw::Float>(e,lb[1],ub[1], lb[2],ub[2], 1-cells,1+cells-1);
			break;
	}

	StripCopyProtocol(naxis(axis),-1,planeIn,planeOut,true);

	delete planeOut;
	delete planeIn;
}

void Field::UpwardDeposit(const axisSpec& axis,const Element& e,tw::Int cells)
{
	// deposits work by sending both ghost cells and edge cells one way and adding.
	// Before accepting the return message these cells are zeroed.
	// The return message then effectively overwrites them with the correct data.
	Slice<tw::Float> *planeIn,*planeOut;

	switch (axis)
	{
		case tAxis:
			break;
		case xAxis:
			if (dim[1]==1) return;
			planeIn = new Slice<tw::Float>(e,1-cells,1+cells-1, lb[2],ub[2], lb[3],ub[3]);
			planeOut = new Slice<tw::Float>(e,dim[1]-cells+1,dim[1]+cells, lb[2],ub[2], lb[3],ub[3]);
			break;
		case yAxis:
			if (dim[2]==1) return;
			planeIn = new Slice<tw::Float>(e,lb[1],ub[1], 1-cells,1+cells-1, lb[3],ub[3]);
			planeOut = new Slice<tw::Float>(e,lb[1],ub[1], dim[2]-cells+1,dim[2]+cells, lb[3],ub[3]);
			break;
		case zAxis:
			if (dim[3]==1) return;
			planeIn = new Slice<tw::Float>(e,lb[1],ub[1], lb[2],ub[2], 1-cells,1+cells-1);
			planeOut = new Slice<tw::Float>(e,lb[1],ub[1], lb[2],ub[2], dim[3]-cells+1,dim[3]+cells);
			break;
	}

	StripCopyProtocol(naxis(axis),1,planeIn,planeOut,true);

	delete planeOut;
	delete planeIn;
}

void Field::CopyFromNeighbors(const Element& e)
{
	DownwardCopy(xAxis,e,1);
	UpwardCopy(xAxis,e,1);

	DownwardCopy(yAxis,e,1);
	UpwardCopy(yAxis,e,1);

	DownwardCopy(zAxis,e,1);
	UpwardCopy(zAxis,e,1);
}

void Field::DepositFromNeighbors(const Element& e)
{
	tw::Int cells = layers[0];

	if (dim[1]>1)
	{
		DownwardDeposit(xAxis,e,cells);
		if (task->n0[1]!=MPI_PROC_NULL)
		{
			Slice<tw::Float>* plane = new Slice<tw::Float>(e,1-cells,1+cells-1,lb[2],ub[2],lb[3],ub[3]);
			ZeroDataInField<tw::Float>(plane);
			delete plane;
		}
		UpwardDeposit(xAxis,e,cells);
	}

	if (dim[2]>1)
	{
		DownwardDeposit(yAxis,e,cells);
		if (task->n0[2]!=MPI_PROC_NULL)
		{
			Slice<tw::Float>* plane = new Slice<tw::Float>(e,lb[1],ub[1],1-cells,1+cells-1,lb[3],ub[3]);
			ZeroDataInField<tw::Float>(plane);
			delete plane;
		}
		UpwardDeposit(yAxis,e,cells);
	}

	if (dim[3]>1)
	{
		DownwardDeposit(zAxis,e,cells);
		if (task->n0[3]!=MPI_PROC_NULL)
		{
			Slice<tw::Float>* plane = new Slice<tw::Float>(e,lb[1],ub[1],lb[2],ub[2],1-cells,1+cells-1);
			ZeroDataInField<tw::Float>(plane);
			delete plane;
		}
		UpwardDeposit(zAxis,e,cells);
	}
}

Slice<tw::Float>* Field::FormTransposeBlock(const Element& e,const axisSpec& axis1,const axisSpec& axis2,tw::Int start1,tw::Int end1,tw::Int start2,tw::Int end2)
{
	Slice<tw::Float> *ans;
	if (axis1==xAxis)
	{
		if (axis2==yAxis)
			ans = new Slice<tw::Float>(e,start1,end1,start2,end2,lb[3],ub[3]);
		else
			ans = new Slice<tw::Float>(e,start1,end1,lb[2],ub[2],start2,end2);
	}
	if (axis1==yAxis)
	{
		if (axis2==xAxis)
			ans = new Slice<tw::Float>(e,start2,end2,start1,end1,lb[3],ub[3]);
		else
			ans = new Slice<tw::Float>(e,lb[1],ub[1],start1,end1,start2,end2);
	}
	if (axis1==zAxis)
	{
		if (axis2==xAxis)
			ans = new Slice<tw::Float>(e,start2,end2,lb[2],ub[2],start1,end1);
		else
			ans = new Slice<tw::Float>(e,lb[1],ub[1],start2,end2,start1,end1);
	}
	return ans;
}

void Field::Transpose(const Element& e,const axisSpec& axis1,const axisSpec& axis2,Field *target,tw::Int inversion)
{
	// We have in mind a matrix whose rows are the NODES arranged along axis1,
	// and whose columns are BLOCKS arranged along axis2.
	// In this picture, the global image of the data does not change.
	// Instead, the nodes are lined up along a new axis.
	// We choose the blocks so that #blocks <= #nodes.
	// 'target' is an externally owned field that will receive the transposed data
	// 'inversion' is 1 if forward transpose, -1 if reverse transpose

	// Ghost cell policy:
	// We send the ghost cells in the directions perp. to 'axis1'
	// The ghost cells are not sent in the direction parallel to 'axis1'
	// However, the ghost cells are sent in all directions if inversion=-1

	tw::Int i,j,k,l,s;
	tw::Int nodes,thisNode,dim1,dim2,dN0,dN1,cellsPerBlock,cellsRemaining,fullBlocks;
	Slice<tw::Float>* block;
	std::vector<tw::Int> blockSize,start,end;
	DiscreteSpace transposedSpace;

	dim1 = task->localCells[naxis(axis1)];
	dim2 = task->localCells[naxis(axis2)];
	nodes = task->domains[naxis(axis1)];
	thisNode = task->strip[naxis(axis1)].Get_rank(); // assumes rank=coord

	if (inversion==1)
	{
		dN0 = 1;
		dN1 = dim1;
	}
	else
	{
		dN0 = 0;
		dN1 = dim1 + 1;
	}

	cellsPerBlock = dim2/nodes + 1;
	fullBlocks = dim2/cellsPerBlock;
	cellsRemaining = dim2 - cellsPerBlock*fullBlocks;

	blockSize.resize(nodes);
	start.resize(nodes);
	end.resize(nodes);
	for (i=0;i<nodes;i++)
	{
		start[i] = i*cellsPerBlock;
		if (i<fullBlocks)
			blockSize[i] = cellsPerBlock;
		if (i==fullBlocks)
			blockSize[i] = cellsRemaining;
		if (i>fullBlocks)
			blockSize[i] = 0;
		end[i] = start[i] + blockSize[i] + 1;
	}

	if (inversion==1)
	{
		if (axis1==xAxis)
		{
			if (axis2==yAxis)
				transposedSpace.Resize(dim1*nodes,cellsPerBlock,dim[3],corner,size);
			else
				transposedSpace.Resize(dim1*nodes,dim[2],cellsPerBlock,corner,size);
		}
		if (axis1==yAxis)
		{
			if (axis2==xAxis)
				transposedSpace.Resize(cellsPerBlock,dim1*nodes,dim[3],corner,size);
			else
				transposedSpace.Resize(dim[1],dim1*nodes,cellsPerBlock,corner,size);
		}
		if (axis1==zAxis)
		{
			if (axis2==xAxis)
				transposedSpace.Resize(cellsPerBlock,dim[2],dim1*nodes,corner,size);
			else
				transposedSpace.Resize(dim[1],cellsPerBlock,dim1*nodes,corner,size);
		}
		target->Initialize(num[0],transposedSpace,task);
	}

  // The message passing pattern is to have simultaneous exchanges between pairs.
  // Each pair is made up of a sub-diagonal and corresponding super-diagonal element

	// First do the diagonal, which involves no message passing

	if (blockSize[thisNode]>0)
	{
		j = thisNode;
		block = FormTransposeBlock(e,axis1,axis2,dN0,dN1,start[j],end[j]);
        if (inversion==-1)
        {
            block->Translate(axis1,j*dim1);
            block->Translate(axis2,-j*cellsPerBlock);
            target->LoadDataIntoImage<tw::Float>(block);
            block->Translate(axis1,-j*dim1);
            block->Translate(axis2,j*cellsPerBlock);
            SaveDataFromImage<tw::Float>(block);
        }
        else
        {
            LoadDataIntoImage<tw::Float>(block);
            block->Translate(axis1,j*dim1);
            block->Translate(axis2,-j*cellsPerBlock);
            target->SaveDataFromImage<tw::Float>(block);
        }
		delete block;
	}

    // Now exchange data between sub/super diagonal pairs

	for (i=1;i<nodes;i++) // Loop over sub-diagonals
        for (k=0;k<2;k++) // each sub-diagonal has 2 independent groups
            for (l=k*i;l<nodes-i;l+=2*i) // Loop over super-blocks of this independent group
                for (s=0;s<i && s<nodes-l-i;s++) // Loop over blocks in each super-block
				{
					// sub-diagonal : (l+s , l+s+i) = (old node , old block) = (new block , new node)
					// super-diagonal : (l+s+i , l+s) = (old node , old block) = (new block , new node)
					if (thisNode==l+s+i)
					{
						// My role as super-diagonal old node
						if (blockSize[l+s]>0)
						{
							block = FormTransposeBlock(e,axis1,axis2,dN0,dN1,start[l+s],end[l+s]);
							if (inversion==-1)
							{
								task->strip[naxis(axis1)].Recv(block->Buffer(),block->BufferSize(),l+s);
								SaveDataFromImage<tw::Float>(block);
							}
							else
							{
								LoadDataIntoImage<tw::Float>(block);
								task->strip[naxis(axis1)].Send(block->Buffer(),block->BufferSize(),l+s);
							}
							delete block;
						}

						// My role as sub-diagonal new node
						if (blockSize[thisNode]>0)
						{
							block = FormTransposeBlock(e,axis1,axis2,dN0,dN1,start[thisNode],end[thisNode]);
							block->Translate(axis1,(l+s)*dim1);
							block->Translate(axis2,-start[thisNode]);
							if (inversion==-1)
							{
								target->LoadDataIntoImage<tw::Float>(block);
								task->strip[naxis(axis1)].Send(block->Buffer(),block->BufferSize(),l+s);
							}
							else
							{
								task->strip[naxis(axis1)].Recv(block->Buffer(),block->BufferSize(),l+s);
								target->SaveDataFromImage<tw::Float>(block);
							}
							delete block;
						}
					}
					if (thisNode==l+s)
					{
						// My role as super-diagonal new node
						if (blockSize[thisNode]>0)
						{
							block = FormTransposeBlock(e,axis1,axis2,dN0,dN1,start[thisNode],end[thisNode]);
							block->Translate(axis1,(l+s+i)*dim1);
							block->Translate(axis2,-start[thisNode]);
							if (inversion==-1)
							{
								target->LoadDataIntoImage<tw::Float>(block);
								task->strip[naxis(axis1)].Send(block->Buffer(),block->BufferSize(),l+s+i);
							}
							else
							{
								task->strip[naxis(axis1)].Recv(block->Buffer(),block->BufferSize(),l+s+i);
								target->SaveDataFromImage<tw::Float>(block);
							}
							delete block;
						}

						// My role as sub-diagonal old node
						if (blockSize[l+s+i]>0)
						{
							block = FormTransposeBlock(e,axis1,axis2,dN0,dN1,start[l+s+i],end[l+s+i]);
							if (inversion==-1)
							{
								task->strip[naxis(axis1)].Recv(block->Buffer(),block->BufferSize(),l+s+i);
								SaveDataFromImage<tw::Float>(block);
							}
							else
							{
								LoadDataIntoImage<tw::Float>(block);
								task->strip[naxis(axis1)].Send(block->Buffer(),block->BufferSize(),l+s+i);
							}
							delete block;
						}
					}
				}
}

/////////////////
//             //
//  Smoothing  //
//             //
/////////////////


void Field::SmoothingPass(const Element& e,const MetricSpace& ds,const tw::Float& X0,const tw::Float& X1,const tw::Float& X2)
{
	tw::Int ax,i,c;
	tw::Float ansNow,temp;
	axisSpec axs[4] = { tAxis, xAxis, yAxis, zAxis };

	for (c=e.low;c<=e.high;c++)
		for (ax=1;ax<=3;ax++)
			if (dim[ax]>1)
				for (auto s : StripRange(*this,ax,strongbool::yes))
				{
					temp = (*this)(s,0,c);
					for (i=1;i<=dim[ax];i++)
					{
						ansNow = X0*temp + X1*(*this)(s,i,c) + X2*(*this)(s,i+1,c);
						temp = (*this)(s,i,c);
						(*this)(s,i,c) = ansNow;
					}
					bc0(ax,c).ForcingOperation(&(*this)(s,lb[ax],c),stride[ax],0.0);
					bc1(ax,c).ForcingOperation(&(*this)(s,ub[ax],c),stride[ax],0.0);
				}

	for (ax=1;ax<=3;ax++)
	{
		UpwardCopy(axs[ax],e,1);
		DownwardCopy(axs[ax],e,1);
	}
}

void Field::Smooth(const Element& e,const MetricSpace& ds,tw::Int smoothPasses,tw::Int compPasses)
{
	tw::Int ipass;
	for (ipass=0;ipass<smoothPasses;ipass++)
		SmoothingPass(e,ds,0.25,0.5,0.25);
	for (ipass=0;ipass<compPasses;ipass++)
		SmoothingPass(e,ds,-1.25,3.5,-1.25);
}

void Field::ReadData(std::ifstream& inFile)
{
	DiscreteSpace ds;
	ds.ReadData(inFile);
	inFile.read((char*)&packedAxis,sizeof(tw::Int));
	Initialize(num[0],ds,task,enumaxis(packedAxis));
	inFile.read((char*)&bc0(0,0),sizeof(BoundaryCondition)*num[0]*4);
	inFile.read((char*)&bc1(0,0),sizeof(BoundaryCondition)*num[0]*4);
	inFile.read((char *)&array[0],sizeof(tw::Float)*totalCells*num[0]);
}

void Field::WriteData(std::ofstream& outFile)
{
	DiscreteSpace::WriteData(outFile);
	outFile.write((char *)&packedAxis,sizeof(tw::Int));
	outFile.write((char*)&bc0(0,0),sizeof(BoundaryCondition)*num[0]*4);
	outFile.write((char*)&bc1(0,0),sizeof(BoundaryCondition)*num[0]*4);
	outFile.write((char *)&array[0],sizeof(tw::Float)*totalCells*num[0]);
}

void CopyBoundaryConditions(Field& dst,const Element& dstElement,Field& src,const Element& srcElement)
{
	tw::Int n;
	n = dstElement.Components();
	if (n!=srcElement.Components())
		throw tw::FatalError("Incompatible copy operation (boundary conditions)");
	dst.bc0 = src.bc0;
	dst.bc1 = src.bc1;
}

void Swap(Field& f1,Field& f2)
{
	tw::Int i;
	tw::Float temp;
	for (i=0;i<f1.array.size();i++)
	{
		temp = f1.array[i];
		f1.array[i] = f2.array[i];
		f2.array[i] = temp;
	}
}

void Field::Swap(const Element& e1,const Element& e2)
{
	tw::Int c;
	tw::Float temp;
	for (c=0;c<e1.Components();c++)
	{

		for(auto cell : EntireCellRange(*this))
		{
			temp = (*this)(cell,e1.low+c);
			(*this)(cell,e1.low+c) = (*this)(cell,e2.low+c);
			(*this)(cell,e2.low+c) = temp;
		}
	}
}

void CopyFieldData(Field& dst,const Element& e_dst,Field& src,const Element& e_src)
{
	tw::Int c;
	for (c=0;c<e_dst.Components();c++)
		for (auto cell : CellRange(dst,true))
			dst(cell,e_dst.low+c) = src(cell,e_src.low+c);
}

void AddFieldData(Field& dst,const Element& e_dst,Field& src,const Element& e_src)
{
	tw::Int c;
	for (c=0;c<e_dst.Components();c++)
		for (auto cell : CellRange(dst,true))
			dst(cell,e_dst.low+c) += src(cell,e_src.low+c);
}

void AddMulFieldData(Field& dst,const Element& e_dst,Field& src,const Element& e_src,tw::Float mul)
{
	tw::Int c;
	for (c=0;c<e_dst.Components();c++)
		for (auto cell : CellRange(dst,true))
			dst(cell,e_dst.low+c) += mul*src(cell,e_src.low+c);
}


//////////////////////////
//                      //
// SCALAR FIELD METHODS //
//                      //
//////////////////////////


tw::Float ScalarField::AxialEigenvalue(tw::Int z)
{
	// eigenvalues for sine and cosine transforms are the same
	return eigenvalue_FST(task->GlobalCellIndex(z,3)-1,task->globalCells[3],freq.z);
}

tw::Float ScalarField::Eigenvalue(tw::Int x,tw::Int y)
{
	// eigenvalues for sine and cosine transforms are the same
	return eigenvalue_FST(task->GlobalCellIndex(x,1)-1,task->globalCells[1],freq.x) + eigenvalue_FST(task->GlobalCellIndex(y,2)-1,task->globalCells[2],freq.y);
}

tw::Float ScalarField::CyclicEigenvalue(tw::Int x,tw::Int y)
{
	tw::Float ans;
	x = task->GlobalCellIndex(x,1);
	y = task->GlobalCellIndex(y,2);

	ans = eigenvalue_RFFT(x-1,task->globalCells[1],freq.x);

	if (x==1 || x==2)
		ans += eigenvalue_RFFT(y-1,task->globalCells[2],freq.y);
	else
		ans += eigenvalue_CFFT(y-1,task->globalCells[2],freq.y);

	return ans;
}

tw::Float ScalarField::CyclicEigenvalue(tw::Int x,tw::Int y,tw::Int z)
{
	tw::Float ans;
	x = task->GlobalCellIndex(x,1);
	y = task->GlobalCellIndex(y,2);
	z = task->GlobalCellIndex(z,3);

	ans = eigenvalue_RFFT(x-1,task->globalCells[1],freq.x);

	if (x==1 || x==2)
		ans += eigenvalue_RFFT(y-1,task->globalCells[2],freq.y);
	else
		ans += eigenvalue_CFFT(y-1,task->globalCells[2],freq.y);

	if ((x==1 || x==2) && (y==1 || y==2))
		ans += eigenvalue_RFFT(z-1,task->globalCells[3],freq.z);
	else
		ans += eigenvalue_CFFT(z-1,task->globalCells[3],freq.z);

	return ans;
}

void ScalarField::AxialSineTransform()
{
	tw::Int i,j;
	Field T;

	if (task->globalCells[3]>1)
	{
		Transpose(zAxis,xAxis,&T,1);
		for (j=T.N0(2);j<=T.N1(2);j++)
			for (i=T.N0(1);i<=T.N1(1);i++)
				SineTransform( &T(i,j,1,0), task->globalCells[3], T.Stride(3), 1 );
		Transpose(zAxis,xAxis,&T,-1);
	}
}

void ScalarField::InverseAxialSineTransform()
{
	tw::Int i,j;
	Field T;

	if (task->globalCells[3]>1)
	{
		Transpose(zAxis,xAxis,&T,1);
		for (j=T.N0(2);j<=T.N1(2);j++)
			for (i=T.N0(1);i<=T.N1(1);i++)
			{
				SineTransform( &T(i,j,1,0), task->globalCells[3], T.Stride(3), -1 );
				T(i,j,0,0) = -T(i,j,2,0);
				T(i,j,task->globalCells[3]+1,0) = 0.0;
			}
		Transpose(zAxis,xAxis,&T,-1);
	}
}

void ScalarField::TransverseCosineTransform()
{
	tw::Int i,j,k;
	Field T;

	if (task->globalCells[1]>1)
	{
		Transpose(xAxis,zAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (j=T.N0(2);j<=T.N1(2);j++)
				CosineTransform( &T(1,j,k,0), task->globalCells[1], T.Stride(1), 1 );
		Transpose(xAxis,zAxis,&T,-1);
	}

	if (task->globalCells[2]>1)
	{
		Transpose(yAxis,zAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (i=T.N0(1);i<=T.N1(1);i++)
				CosineTransform( &T(i,1,k,0), task->globalCells[2], T.Stride(2), 1 );
		Transpose(yAxis,zAxis,&T,-1);
	}
}

void ScalarField::InverseTransverseCosineTransform()
{
	tw::Int i,j,k;
	Field T;

	if (task->globalCells[2]>1)
	{
		Transpose(yAxis,zAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (i=T.N0(1);i<=T.N1(1);i++)
			{
				CosineTransform( &T(i,1,k,0), task->globalCells[2], T.Stride(2), -1 );
				T(i,0,k,0) = T(i,1,k,0);
				T(i,task->globalCells[2]+1,k,0) = T(i,task->globalCells[2],k,0);
			}
		Transpose(yAxis,zAxis,&T,-1);
	}

	if (task->globalCells[1]>1)
	{
		Transpose(xAxis,zAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (j=T.N0(2);j<=T.N1(2);j++)
			{
				CosineTransform( &T(1,j,k,0), task->globalCells[1], T.Stride(1), -1 );
				T(0,j,k,0) = T(1,j,k,0);
				T(task->globalCells[1]+1,j,k,0) = T(task->globalCells[1],j,k,0);
			}
		Transpose(xAxis,zAxis,&T,-1);
	}
}

void ScalarField::TransverseSineTransform()
{
	tw::Int i,j,k;
	Field T;

	if (task->globalCells[1]>1)
	{
		Transpose(xAxis,zAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (j=T.N0(2);j<=T.N1(2);j++)
				SineTransform( &T(1,j,k,0), task->globalCells[1], T.Stride(1), 1 );
		Transpose(xAxis,zAxis,&T,-1);
	}

	if (task->globalCells[2]>1)
	{
		Transpose(yAxis,zAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (i=T.N0(1);i<=T.N1(1);i++)
				SineTransform( &T(i,1,k,0), task->globalCells[2], T.Stride(2), 1 );
		Transpose(yAxis,zAxis,&T,-1);
	}
}

void ScalarField::InverseTransverseSineTransform()
{
	tw::Int i,j,k;
	Field T;

	if (task->globalCells[2]>1)
	{
		Transpose(yAxis,zAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (i=T.N0(1);i<=T.N1(1);i++)
			{
				SineTransform( &T(i,1,k,0), task->globalCells[2], T.Stride(2), -1 );
				T(i,0,k,0) = -T(i,2,k,0);
				T(i,task->globalCells[2]+1,k,0) = 0.0;
			}
		Transpose(yAxis,zAxis,&T,-1);
	}

	if (task->globalCells[1]>1)
	{
		Transpose(xAxis,zAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (j=T.N0(2);j<=T.N1(2);j++)
			{
				SineTransform( &T(1,j,k,0), task->globalCells[1], T.Stride(1), -1 );
				T(0,j,k,0) = -T(2,j,k,0);
				T(task->globalCells[1]+1,j,k,0) = 0.0;
			}
		Transpose(xAxis,zAxis,&T,-1);
	}
}

void ScalarField::TransverseFFT()
{
	tw::Int i,j,k;
	Field T;

	if (task->globalCells[1]>1)
	{
		Transpose(xAxis,zAxis,&T,1);
		const tw::Int yN0=T.N0(2),yN1=T.N1(2),zN0=T.N0(3),zN1=T.N1(3);
		#pragma omp parallel for private(j,k) collapse(2) schedule(static)
		for (j=yN0;j<=yN1;j++)
			for (k=zN0;k<=zN1;k++)
				RealFFT( &T(1,j,k,0), task->globalCells[1], T.Stride(1), 1);
		Transpose(xAxis,zAxis,&T,-1);
	}

	if (task->globalCells[2]>1)
	{
		Transpose(yAxis,zAxis,&T,1);
		const tw::Int xDim=T.Dim(1),zN0=T.N0(3),zN1=T.N1(3);
		#pragma omp parallel for private(i,k) collapse(2) schedule(static)
		for (i=1;i<=xDim;i+=2) // can't include ghost cells due to complex numbers; instead do copy ops below
			for (k=zN0;k<=zN1;k++)
			{
				if (task->GlobalCellIndex(i,1)==1)
				{
					RealFFT( &T(i,1,k,0), task->globalCells[2], T.Stride(2), 1);
					RealFFT( &T(i+1,1,k,0), task->globalCells[2], T.Stride(2), 1);
				}
				else
					ComplexFFT( &T(i,1,k,0), &T(i+1,1,k,0), task->globalCells[2], T.Stride(2), 1.0);
			}
		T.UpwardCopy(xAxis,1);
		T.DownwardCopy(xAxis,1);
		Transpose(yAxis,zAxis,&T,-1);
	}
}

void ScalarField::InverseTransverseFFT()
{
	tw::Int i,j,k;
	Field T;

	if (task->globalCells[2]>1)
	{
		Transpose(yAxis,zAxis,&T,1);
		const tw::Int xDim=T.Dim(1),zN0=T.N0(3),zN1=T.N1(3);
		#pragma omp parallel for private(i,k) collapse(2) schedule(static)
		for (i=1;i<=xDim;i+=2) // can't include ghost cells due to complex numbers; instead do copy ops below
			for (k=zN0;k<=zN1;k++)
			{
				if (task->GlobalCellIndex(i,1)==1)
				{
					RealFFT( &T(i,1,k,0), task->globalCells[2], T.Stride(2), -1);
					T(i,0,k,0) = T(i,task->globalCells[2],k,0);
					T(i,task->globalCells[2]+1,k,0) = T(i,1,k,0);
					RealFFT( &T(i+1,1,k,0), task->globalCells[2], T.Stride(2), -1);
					T(i+1,0,k,0) = T(i+1,task->globalCells[2],k,0);
					T(i+1,task->globalCells[2]+1,k,0) = T(i+1,1,k,0);
				}
				else
				{
					ComplexFFT( &T(i,1,k,0), &T(i+1,1,k,0), task->globalCells[2], T.Stride(2), -1.0);
					T(i,0,k,0) = T(i,task->globalCells[2],k,0);
					T(i,task->globalCells[2]+1,k,0) = T(i,1,k,0);
					T(i+1,0,k,0) = T(i+1,task->globalCells[2],k,0);
					T(i+1,task->globalCells[2]+1,k,0) = T(i+1,1,k,0);
				}
			}
		T.UpwardCopy(xAxis,1);
		T.DownwardCopy(xAxis,1);
		Transpose(yAxis,zAxis,&T,-1);
	}

	if (task->globalCells[1]>1)
	{
		Transpose(xAxis,zAxis,&T,1);
		const tw::Int yN0=T.N0(2),yN1=T.N1(2),zN0=T.N0(3),zN1=T.N1(3);
		#pragma omp parallel for private(j,k) collapse(2) schedule(static)
		for (j=yN0;j<=yN1;j++)
			for (k=zN0;k<=zN1;k++)
			{
				RealFFT( &T(1,j,k,0), task->globalCells[1], T.Stride(1), -1);
				T(0,j,k,0) = T(task->globalCells[1],j,k,0);
				T(task->globalCells[1]+1,j,k,0) = T(1,j,k,0);
			}
		Transpose(xAxis,zAxis,&T,-1);
	}
}


///////////////////////////
//                       //
// COMPLEX FIELD METHODS //
//                       //
///////////////////////////


tw::Float ComplexField::CyclicEigenvalue(tw::Int x,tw::Int y)
{
	x = task->GlobalCellIndex(x,1);
	y = task->GlobalCellIndex(y,2);
	return eigenvalue_CFFT(x-1,task->globalCells[1],freq.x) + eigenvalue_CFFT(y-1,task->globalCells[2],freq.y);
}

tw::Float ComplexField::CyclicEigenvalue(tw::Int x,tw::Int y,tw::Int z)
{
	x = task->GlobalCellIndex(x,1);
	y = task->GlobalCellIndex(y,2);
	z = task->GlobalCellIndex(z,3);
	return eigenvalue_CFFT(x-1,task->globalCells[1],freq.x) + eigenvalue_CFFT(y-1,task->globalCells[2],freq.y) + eigenvalue_CFFT(z-1,task->globalCells[3],freq.z);
}

void ComplexField::TransverseFFT()
{
	tw::Int i,j,k;
	Field T;

	if (task->globalCells[1]>1)
	{
		Transpose(xAxis,zAxis,&T,1);
		const tw::Int yN0=T.N0(2),yN1=T.N1(2),zN0=T.N0(3),zN1=T.N1(3);
		#pragma omp parallel for private(j,k) collapse(2) schedule(static)
		for (j=yN0;j<=yN1;j++)
			for (k=zN0;k<=zN1;k++)
				ComplexFFT( &T(1,j,k,0), &T(1,j,k,1), task->globalCells[1], T.Stride(1), 1.0 );
		Transpose(xAxis,zAxis,&T,-1);
	}

	if (task->globalCells[2]>1)
	{
		Transpose(yAxis,zAxis,&T,1);
		const tw::Int xN0=T.N0(1),xN1=T.N1(1),zN0=T.N0(3),zN1=T.N1(3);
		#pragma omp parallel for private(i,k) collapse(2) schedule(static)
		for (i=xN0;i<=xN1;i++)
			for (k=zN0;k<=zN1;k++)
				ComplexFFT( &T(i,1,k,0), &T(i,1,k,1), task->globalCells[2], T.Stride(2), 1.0 );
		Transpose(yAxis,zAxis,&T,-1);
	}
}

void ComplexField::InverseTransverseFFT()
{
	tw::Int i,j,k;
	Field T;

	if (task->globalCells[2]>1)
	{
		Transpose(yAxis,zAxis,&T,1);
		const tw::Int xN0=T.N0(1),xN1=T.N1(1),zN0=T.N0(3),zN1=T.N1(3);
		#pragma omp parallel for private(i,k) collapse(2) schedule(static)
		for (i=xN0;i<=xN1;i++)
			for (k=zN0;k<=zN1;k++)
			{
				ComplexFFT( &T(i,1,k,0), &T(i,1,k,1), task->globalCells[2], T.Stride(2), -1.0 );
				T(i,0,k,0) = T(i,task->globalCells[2],k,0);
				T(i,0,k,1) = T(i,task->globalCells[2],k,1);
				T(i,task->globalCells[2]+1,k,0) = T(i,1,k,0);
				T(i,task->globalCells[2]+1,k,1) = T(i,1,k,1);
			}
		Transpose(yAxis,zAxis,&T,-1);
	}

	if (task->globalCells[1]>1)
	{
		Transpose(xAxis,zAxis,&T,1);
		const tw::Int yN0=T.N0(2),yN1=T.N1(2),zN0=T.N0(3),zN1=T.N1(3);
		#pragma omp parallel for private(j,k) collapse(2) schedule(static)
		for (j=yN0;j<=yN1;j++)
			for (k=zN0;k<=zN1;k++)
			{
				ComplexFFT( &T(1,j,k,0), &T(1,j,k,1), task->globalCells[1], T.Stride(1), -1.0 );
				T(0,j,k,0) = T(task->globalCells[1],j,k,0);
				T(0,j,k,1) = T(task->globalCells[1],j,k,1);
				T(task->globalCells[1]+1,j,k,0) = T(1,j,k,0);
				T(task->globalCells[1]+1,j,k,1) = T(1,j,k,1);
			}
		Transpose(xAxis,zAxis,&T,-1);
	}
}

void ComplexField::FFT()
{
	tw::Int i,j,k;
	Field T;
	std::slice s;
	std::valarray<tw::Float> temp;

	if (task->globalCells[1]>1)
	{
		Transpose(xAxis,yAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (j=T.N0(2);j<=T.N1(2);j++)
			{
				ComplexFFT( &T(1,j,k,0), &T(1,j,k,1), task->globalCells[1], T.Stride(1), 1.0 );
			}
		Transpose(xAxis,yAxis,&T,-1);
	}

	if (task->globalCells[2]>1)
	{
		Transpose(yAxis,xAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (i=T.N0(1);i<=T.N1(1);i++)
			{
				ComplexFFT( &T(i,1,k,0), &T(i,1,k,1), task->globalCells[2], T.Stride(2), 1.0 );
			}
		Transpose(yAxis,xAxis,&T,-1);
	}

	if (task->globalCells[3]>1)
	{
		Transpose(zAxis,xAxis,&T,1);
		for (j=T.N0(2);j<=T.N1(2);j++)
			for (i=T.N0(1);i<=T.N1(1);i++)
			{
				ComplexFFT( &T(i,j,1,0), &T(i,j,1,1), task->globalCells[3], T.Stride(3), 1.0 );
			}
		Transpose(zAxis,xAxis,&T,-1);
	}
}

void ComplexField::InverseFFT()
{
	tw::Int i,j,k;
	Field T;

	if (task->globalCells[3]>1)
	{
		Transpose(zAxis,xAxis,&T,1);
		for (j=T.N0(2);j<=T.N1(2);j++)
			for (i=T.N0(1);i<=T.N1(1);i++)
			{
				ComplexFFT( &T(i,j,1,0), &T(i,j,1,1), task->globalCells[3], T.Stride(3), -1.0);
				T(i,j,0,0) = T(i,j,task->globalCells[3],0);
				T(i,j,0,1) = T(i,j,task->globalCells[3],1);
				T(i,j,task->globalCells[3]+1,0) = T(i,j,1,0);
				T(i,j,task->globalCells[3]+1,1) = T(i,j,1,1);
			}
		Transpose(zAxis,xAxis,&T,-1);
	}

	if (task->globalCells[2]>1)
	{
		Transpose(yAxis,xAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (i=T.N0(1);i<=T.N1(1);i++)
			{
				ComplexFFT( &T(i,1,k,0), &T(i,1,k,1), task->globalCells[2], T.Stride(2), -1.0 );
				T(i,0,k,0) = T(i,task->globalCells[2],k,0);
				T(i,0,k,1) = T(i,task->globalCells[2],k,1);
				T(i,task->globalCells[2]+1,k,0) = T(i,1,k,0);
				T(i,task->globalCells[2]+1,k,1) = T(i,1,k,1);
			}
		Transpose(yAxis,xAxis,&T,-1);
	}

	if (task->globalCells[1]>1)
	{
		Transpose(xAxis,yAxis,&T,1);
		for (k=T.N0(3);k<=T.N1(3);k++)
			for (j=T.N0(2);j<=T.N1(2);j++)
			{
				ComplexFFT( &T(1,j,k,0), &T(1,j,k,1), task->globalCells[1], T.Stride(1), -1.0 );
				T(1,j,k,0) = T(task->globalCells[1],j,k,0);
				T(1,j,k,1) = T(task->globalCells[1],j,k,1);
				T(task->globalCells[1]+1,j,k,0) = T(1,j,k,0);
				T(task->globalCells[1]+1,j,k,1) = T(1,j,k,1);
			}
		Transpose(xAxis,yAxis,&T,-1);
	}
}
