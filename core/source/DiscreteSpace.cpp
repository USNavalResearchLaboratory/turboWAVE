#include "meta_base.h"

DiscreteSpace::DiscreteSpace()
{
	ignorable[1] = 0;
	ignorable[2] = 0;
	ignorable[3] = 0;
}

DiscreteSpace::DiscreteSpace(tw::Float dt0,tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size,tw::Int ghostCellLayers)
{
	SetupTimeInfo(dt0);
	ignorable[1] = 0;
	ignorable[2] = 0;
	ignorable[3] = 0;
	Resize(x,y,z,corner,size,ghostCellLayers);
}

void DiscreteSpace::Resize(tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size,tw::Int ghostCellLayers)
{
	// init unused elements so later assignments are memcheck clean.
	dim[0] = num[0] = 0;
	ignorable[0] = 0;
	lfg[0] = ufg[0] = lng[0] = ung[0] = 0;
	decodingStride[0] = encodingStride[0] = 0;

	layers[0] = ghostCellLayers;
	dim[1] = x;
	dim[2] = y;
	dim[3] = z;

	for (tw::Int i=1;i<=3;i++)
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

	decodingStride[1] = encodingStride[1] = num[2]*num[3];
	decodingStride[2] = encodingStride[2] = num[3];
	decodingStride[3] = encodingStride[3] = 1;

	for (tw::Int i=1;i<=3;i++)
	{
		if (dim[i]==1)
		{
			encodingStride[i] = 0;
			ignorable[i] = 1;
		}
	}

	this->corner = corner;
	this->size = size;
	spacing.x = size.x/dim[1];
	spacing.y = size.y/dim[2];
	spacing.z = size.z/dim[3];
	freq.x = 1/spacing.x;
	freq.y = 1/spacing.y;
	freq.z = 1/spacing.z;
}

void DiscreteSpace::Resize(tw::Int x,tw::Int y,tw::Int z,const tw::vec3& corner,const tw::vec3& size)
{
	Resize(x,y,z,corner,size,2); // default to 2 ghost cell layers
}

void DiscreteSpace::ReadCheckpoint(std::ifstream& inFile)
{
	inFile.read((char*)&dt,sizeof(dt));
	inFile.read((char*)&dth,sizeof(dth));
	inFile.read((char*)&dti,sizeof(dti));
	inFile.read((char*)&corner,sizeof(tw::vec3));
	inFile.read((char*)&size,sizeof(tw::vec3));
	inFile.read((char*)&globalCorner,sizeof(tw::vec3));
	inFile.read((char*)&globalSize,sizeof(tw::vec3));
	inFile.read((char*)&spacing,sizeof(tw::vec3));
	inFile.read((char*)&freq,sizeof(tw::vec3));
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
	outFile.write((char*)&dt,sizeof(dt));
	outFile.write((char*)&dth,sizeof(dth));
	outFile.write((char*)&dti,sizeof(dti));
	outFile.write((char*)&corner,sizeof(tw::vec3));
	outFile.write((char*)&size,sizeof(tw::vec3));
	outFile.write((char*)&globalCorner,sizeof(tw::vec3));
	outFile.write((char*)&globalSize,sizeof(tw::vec3));
	outFile.write((char*)&spacing,sizeof(tw::vec3));
	outFile.write((char*)&freq,sizeof(tw::vec3));
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
