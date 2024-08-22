#include "meta_base.h"

DiscreteSpace::DiscreteSpace()
{
	ignorable[0] = 0;
	ignorable[1] = 0;
	ignorable[2] = 0;
	ignorable[3] = 0;
}

DiscreteSpace::DiscreteSpace(tw::Int xDim,tw::Int yDim,tw::Int zDim,const tw::vec4& corner,const tw::vec4& size,tw::Int ghostCellLayers)
{
	const tw::Int ldim[4] = { 1, xDim, yDim, zDim };
	const tw::Int gdim[4] = { 1, xDim, yDim, zDim };
	const tw::Int dom[4] = { 0, 0, 0, 0 };
	Resize(ldim,gdim,dom,corner,size,ghostCellLayers);
}

// DiscreteSpace::DiscreteSpace(tw::Float dt0,Task& task,const tw::vec3& gcorner,const tw::vec3& gsize,tw::Int ghostCellLayers)
// {
// 	SetupTimeInfo(dt0);
// 	Resize(task,gcorner,gsize,ghostCellLayers);
// }

void DiscreteSpace::Resize(Task& task,const tw::vec4& gcorner,const tw::vec4& gsize,tw::Int ghostCellLayers)
{
	Resize(task.localCells,task.globalCells,task.domainIndex,gcorner,gsize,ghostCellLayers);
}

void DiscreteSpace::Resize(const tw::Int dim[4],const tw::Int gdim[4],const tw::Int dom[4],const tw::vec4& gcorner,const tw::vec4& gsize,tw::Int ghostCellLayers)
{
	this->dim[0] = dim[0];
	this->dim[1] = dim[1];
	this->dim[2] = dim[2];
	this->dim[3] = dim[3];

	for (tw::Int i=0;i<4;i++)
	{
		if (dim[i]==1)
		{
			layers[i] = 0;
			lfg[i] = 1;
			ufg[i] = 1;
			lng[i] = 1;
			ung[i] = 1;
		}
		else
		{
			layers[i] = ghostCellLayers;
			lfg[i] = 1 - layers[i];
			ufg[i] = dim[i] + layers[i];
			lng[i] = 0;
			ung[i] = dim[i] + 1;
		}
		num[i] = ufg[i] - lfg[i] + 1;
	}

	decodingStride[0] = num[1]*num[2]*num[3];
	decodingStride[1] = num[2]*num[3];
	decodingStride[2] = num[3];
	decodingStride[3] = 1;

	for (tw::Int i=0;i<4;i++)
	{
		encodingStride[i] = (dim[i]==1 ? 0 : decodingStride[i]);
		ignorable[i] = (dim[i]==1 ? 1 : 0);
	}

	globalCorner = gcorner;
	globalSize = gsize;
	for (tw::Int i=0;i<4;i++)
	{
		spacing[i] = globalSize[i]/gdim[i];
		freq[i] = 1/spacing[i];
		size[i] = dim[i]*spacing[i];
		corner[i] = gcorner[i] + dom[i]*size[i];
	}
}

void DiscreteSpace::ReadCheckpoint(std::ifstream& inFile)
{
	inFile.read((char*)&corner,sizeof(tw::vec4));
	inFile.read((char*)&size,sizeof(tw::vec4));
	inFile.read((char*)&globalCorner,sizeof(tw::vec4));
	inFile.read((char*)&globalSize,sizeof(tw::vec4));
	inFile.read((char*)&spacing,sizeof(tw::vec4));
	inFile.read((char*)&freq,sizeof(tw::vec4));
	inFile.read((char*)num,sizeof(num));
	inFile.read((char*)dim,sizeof(dim));
	inFile.read((char*)lfg,sizeof(lfg));
	inFile.read((char*)ufg,sizeof(ufg));
	inFile.read((char*)lng,sizeof(lng));
	inFile.read((char*)ung,sizeof(ung));
	inFile.read((char*)ignorable,sizeof(ignorable));
	inFile.read((char*)encodingStride,sizeof(encodingStride));
	inFile.read((char*)decodingStride,sizeof(decodingStride));
	inFile.read((char*)layers,sizeof(layers));
}

void DiscreteSpace::WriteCheckpoint(std::ofstream& outFile)
{
	outFile.write((char*)&corner,sizeof(tw::vec4));
	outFile.write((char*)&size,sizeof(tw::vec4));
	outFile.write((char*)&globalCorner,sizeof(tw::vec4));
	outFile.write((char*)&globalSize,sizeof(tw::vec4));
	outFile.write((char*)&spacing,sizeof(tw::vec4));
	outFile.write((char*)&freq,sizeof(tw::vec4));
	outFile.write((char*)num,sizeof(num));
	outFile.write((char*)dim,sizeof(dim));
	outFile.write((char*)lfg,sizeof(lfg));
	outFile.write((char*)ufg,sizeof(ufg));
	outFile.write((char*)lng,sizeof(lng));
	outFile.write((char*)ung,sizeof(ung));
	outFile.write((char*)ignorable,sizeof(ignorable));
	outFile.write((char*)encodingStride,sizeof(encodingStride));
	outFile.write((char*)decodingStride,sizeof(decodingStride));
	outFile.write((char*)layers,sizeof(layers));
}

#ifdef USE_OPENCL

void DiscreteSpace::CellUpdateProtocol(cl_kernel k,cl_command_queue q)
{
	size_t cells = num[1]*num[2]*num[3];
	clEnqueueNDRangeKernel(q,k,1,NULL,&cells,NULL,0,NULL,NULL);
	clFinish(q);
}

void DiscreteSpace::ElementUpdateProtocol(cl_kernel k,cl_command_queue q)
{
	size_t elements = num[0]*num[1]*num[2]*num[3];
	clEnqueueNDRangeKernel(q,k,1,NULL,&elements,NULL,0,NULL,NULL);
	clFinish(q);
}

void DiscreteSpace::LocalUpdateProtocol(cl_kernel k,cl_command_queue q)
{
 	size_t global_offset[3] = { (size_t)layers[1],(size_t)layers[2],(size_t)layers[3] };
	size_t global_work_range[3] = {(size_t)dim[1],(size_t)dim[2],(size_t)dim[3]};
	clEnqueueNDRangeKernel(q,k,3,global_offset,global_work_range,NULL,0,NULL,NULL);
	clFinish(q);
}

void DiscreteSpace::PointUpdateProtocol(cl_kernel k,cl_command_queue q)
{
	size_t global_work_range[3] = { (size_t)num[1],(size_t)num[2],(size_t)num[3] };
	clEnqueueNDRangeKernel(q,k,3,NULL,global_work_range,NULL,0,NULL,NULL);
	clFinish(q);
}

#endif
