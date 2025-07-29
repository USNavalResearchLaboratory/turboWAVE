module;

#include "tw_includes.h"

export module fields:transform;
import :base;

import base;
import dyn_space;
import numerics;
import fft;
import logger;
#include "tw_logger.h"

tw::Float Field::RealAxialEigenvalue(tw::Int z,const DynSpace& ds)
{
	// eigenvalues for sine and cosine transforms are the same
	return fft::eigenvalue_FST(ds.GlobalCellIndex(z,3)-1,ds.GlobalDim(3),freq[3]);
}

tw::Float Field::RealEigenvalue(tw::Int x,tw::Int y,const DynSpace& ds)
{
	// eigenvalues for sine and cosine transforms are the same
	return fft::eigenvalue_FST(ds.GlobalCellIndex(x,1)-1,ds.GlobalDim(1),freq[1]) + fft::eigenvalue_FST(ds.GlobalCellIndex(y,2)-1,ds.GlobalDim(2),freq[2]);
}

tw::Float Field::RealCyclicEigenvalue(tw::Int x,tw::Int y,const DynSpace& ds)
{
	tw::Float ans;
	x = ds.GlobalCellIndex(x,1);
	y = ds.GlobalCellIndex(y,2);

	ans = fft::eigenvalue_RFFT(x-1,ds.GlobalDim(1),freq[1]);

	if (x==1 || x==2)
		ans += fft::eigenvalue_RFFT(y-1,ds.GlobalDim(2),freq[2]);
	else
		ans += fft::eigenvalue_CFFT(y-1,ds.GlobalDim(2),freq[2]);

	return ans;
}

tw::Float Field::RealCyclicEigenvalue(tw::Int x,tw::Int y,tw::Int z,const DynSpace& ds)
{
	tw::Float ans;
	x = ds.GlobalCellIndex(x,1);
	y = ds.GlobalCellIndex(y,2);
	z = ds.GlobalCellIndex(z,3);

	ans = fft::eigenvalue_RFFT(x-1,ds.GlobalDim(1),freq[1]);

	if (x==1 || x==2)
		ans += fft::eigenvalue_RFFT(y-1,ds.GlobalDim(2),freq[2]);
	else
		ans += fft::eigenvalue_CFFT(y-1,ds.GlobalDim(2),freq[2]);

	if ((x==1 || x==2) && (y==1 || y==2))
		ans += fft::eigenvalue_RFFT(z-1,ds.GlobalDim(3),freq[3]);
	else
		ans += fft::eigenvalue_CFFT(z-1,ds.GlobalDim(3),freq[3]);

	return ans;
}


Slice<tw::Float>* Field::FormTransposeBlock(const Rng04& r,const tw::grid::axis& axis1,const tw::grid::axis& axis2,tw::Int start1,tw::Int stop1,tw::Int start2,tw::Int stop2)
{
	tw::node5 beg,end;
	beg[0] = r.b0; end[0] = r.e0;
	beg[4] = r.b4; end[4] = r.e4;
	Slice<tw::Float> *ans;
	if (axis1==tw::grid::x)
	{
		if (axis2==tw::grid::y) {
			beg[1] = start1;
			end[1] = stop1 + 1;
			beg[2] = start2;
			end[2] = stop2 + 1;
			beg[3] = lfg[3];
			end[3] = ufg[3] + 1;
		} else {
			beg[1] = start1;
			end[1] = stop1 + 1;
			beg[2] = lfg[2];
			end[2] = ufg[2] + 1;
			beg[3] = start2;
			end[3] = stop2 + 1;
		}
	}
	if (axis1==tw::grid::y)
	{
		if (axis2==tw::grid::x) {
			beg[1] = start2;
			end[1] = stop2 + 1;
			beg[2] = start1;
			end[2] = stop1 + 1;
			beg[3] = lfg[3];
			end[3] = ufg[3] + 1;
		} else {
			beg[1] = lfg[1];
			end[1] = ufg[1] + 1;
			beg[3] = start1;
			end[3] = stop1 + 1;
			beg[2] = start2;
			end[2] = stop2 + 1;
		}
	}
	if (axis1==tw::grid::z)
	{
		if (axis2==tw::grid::x) {
			beg[1] = start2;
			end[1] = stop2 + 1;
			beg[2] = lfg[2];
			end[2] = ufg[2] + 1;
			beg[3] = start1;
			end[3] = stop1 + 1;
		} else {
			beg[1] = lfg[1];
			end[1] = ufg[1] + 1;
			beg[3] = start2;
			end[3] = stop2 + 1;
			beg[2] = start1;
			end[2] = stop1 + 1;
		}
	}
	return new Slice<tw::Float>(beg,end);
}


/// We have in mind a matrix whose rows consist of NODES arranged along axis1,
/// and whose columns consist of BLOCKS arranged along axis2.
/// In this picture, the global image of the data does not change.
/// Instead, the nodes are lined up along a new axis.
/// We choose the blocks so that #blocks <= #nodes.
/// 'target' is an externally owned field that will receive the transposed data
/// This routine resizes target when the forward transpose is invoked.
/// Upon reverse transposing, the same target should be passed in.
/// 'inversion' is 1 if forward transpose, -1 if reverse transpose
///
/// Ghost cell policy:
/// We send the ghost cells in the directions perp. to 'axis1'
/// The ghost cells are not sent in the direction parallel to 'axis1'
/// However, the ghost cells are sent in all directions if inversion=-1
void Field::Transpose(const Rng04& r,const tw::grid::axis& axis1,const tw::grid::axis& axis2,Field *target,tw::Int inversion)
{
	tw::Int dN0,dN1;
	Slice<tw::Float>* block;

	const tw::Int ax1 = tw::grid::naxis(axis1);
	const tw::Int ax2 = tw::grid::naxis(axis2);
	const tw::Int nodes = task->domains[ax1];
	const tw::Int thisNode = task->strip[ax1].Get_rank(); // assumes rank=coord
	const tw::Int cellsPerBlock = num[ax2]/nodes + 1;
	const tw::Int interiorCellsPerBlock = cellsPerBlock - 2*layers[ax2];
	const tw::Int fullBlocks = num[ax2]/cellsPerBlock;
	const tw::Int cellsRemaining = num[ax2] - cellsPerBlock*fullBlocks;
	std::vector<tw::Int> blockSize(nodes),start(nodes),stop(nodes),offset(nodes);

	if (inversion==1)
	{
		dN0 = 1;
		dN1 = dim[ax1];
	}
	else
	{
		dN0 = 0;
		dN1 = dim[ax1] + 1;
	}

	for (tw::Int i=0;i<nodes;i++)
	{
		start[i] = i*cellsPerBlock + 1 - layers[ax2];
		if (i<fullBlocks)
			blockSize[i] = cellsPerBlock;
		if (i==fullBlocks)
			blockSize[i] = cellsRemaining;
		if (i>fullBlocks)
			blockSize[i] = 0;
		stop[i] = start[i] + blockSize[i] - 1;
		offset[i] = start[i] + layers[ax2] - 1;
	}

	if (inversion==1)
	{
		tw::Int ax3 = 1;
		while (ax3==ax1 || ax3==ax2) ax3++;
		tw::node5 dim_T { r.Times(),0,0,0,r.Components() };
		dim_T[ax1] = dim[ax1]*nodes;
		dim_T[ax2] = interiorCellsPerBlock;
		dim_T[ax3] = dim[ax3];
		target->Initialize(StaticSpace(dim_T,tw::vec4(1,1,1,1),packing,layers),task);
	}

  // The message passing pattern is to have simultaneous exchanges between pairs.
  // Each pair is made up of a sub-diagonal and corresponding super-diagonal element

	// First do the diagonal, which involves no message passing

	if (blockSize[thisNode]>0)
	{
		const tw::Int j = thisNode;
		block = FormTransposeBlock(r,axis1,axis2,dN0,dN1,start[j],stop[j]);
		if (inversion==-1)
		{
	    block->Translate(axis1,j*dim[ax1]);
	    block->Translate(axis2,-offset[j]);
	    target->LoadDataIntoSlice<tw::Float>(block);
	    block->Translate(axis1,-j*dim[ax1]);
	    block->Translate(axis2,offset[j]);
	    SaveDataFromSlice<tw::Float>(block);
		}
		else
		{
	    LoadDataIntoSlice<tw::Float>(block);
	    block->Translate(axis1,j*dim[ax1]);
	    block->Translate(axis2,-offset[j]);
	    target->SaveDataFromSlice<tw::Float>(block);
		}
		delete block;
	}

  // Now exchange data between sub/super diagonal pairs

	for (tw::Int i=1;i<nodes;i++) // Loop over sub-diagonals
		for (tw::Int k=0;k<2;k++) // each sub-diagonal has 2 independent groups
	    for (tw::Int l=k*i;l<nodes-i;l+=2*i) // Loop over super-blocks of this independent group
	      for (tw::Int s=0;s<i && s<nodes-l-i;s++) // Loop over blocks in each super-block
				{
					// sub-diagonal : (l+s , l+s+i) = (old node , old block) = (new block , new node)
					// super-diagonal : (l+s+i , l+s) = (old node , old block) = (new block , new node)
					if (thisNode==l+s+i)
					{
						// My role as super-diagonal old node
						if (blockSize[l+s]>0)
						{
							block = FormTransposeBlock(r,axis1,axis2,dN0,dN1,start[l+s],stop[l+s]);
							if (inversion==-1)
							{
								task->strip[ax1].Recv(block->Buffer(),block->BufferSize(),l+s);
								SaveDataFromSlice<tw::Float>(block);
							}
							else
							{
								LoadDataIntoSlice<tw::Float>(block);
								task->strip[ax1].Send(block->Buffer(),block->BufferSize(),l+s);
							}
							delete block;
						}

						// My role as sub-diagonal new node
						if (blockSize[thisNode]>0)
						{
							block = FormTransposeBlock(r,axis1,axis2,dN0,dN1,start[thisNode],stop[thisNode]);
							block->Translate(axis1,(l+s)*dim[ax1]);
							block->Translate(axis2,-offset[thisNode]);
							if (inversion==-1)
							{
								target->LoadDataIntoSlice<tw::Float>(block);
								task->strip[ax1].Send(block->Buffer(),block->BufferSize(),l+s);
							}
							else
							{
								task->strip[ax1].Recv(block->Buffer(),block->BufferSize(),l+s);
								target->SaveDataFromSlice<tw::Float>(block);
							}
							delete block;
						}
					}
					if (thisNode==l+s)
					{
						// My role as super-diagonal new node
						if (blockSize[thisNode]>0)
						{
							block = FormTransposeBlock(r,axis1,axis2,dN0,dN1,start[thisNode],stop[thisNode]);
							block->Translate(axis1,(l+s+i)*dim[ax1]);
							block->Translate(axis2,-offset[thisNode]);
							if (inversion==-1)
							{
								target->LoadDataIntoSlice<tw::Float>(block);
								task->strip[ax1].Send(block->Buffer(),block->BufferSize(),l+s+i);
							}
							else
							{
								task->strip[ax1].Recv(block->Buffer(),block->BufferSize(),l+s+i);
								target->SaveDataFromSlice<tw::Float>(block);
							}
							delete block;
						}

						// My role as sub-diagonal old node
						if (blockSize[l+s+i]>0)
						{
							block = FormTransposeBlock(r,axis1,axis2,dN0,dN1,start[l+s+i],stop[l+s+i]);
							if (inversion==-1)
							{
								task->strip[ax1].Recv(block->Buffer(),block->BufferSize(),l+s+i);
								SaveDataFromSlice<tw::Float>(block);
							}
							else
							{
								LoadDataIntoSlice<tw::Float>(block);
								task->strip[ax1].Send(block->Buffer(),block->BufferSize(),l+s+i);
							}
							delete block;
						}
					}
				}
}

void Field::RealAxialSineTransform(const Rng04& r,const DynSpace& ds)
{
	Field T;
	if (ds.GlobalDim(3)>1)
	{
		Transpose(r,tw::grid::z,tw::grid::x,&T,1);
		for (auto n=r.b0; n<r.e0; n++) {
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,3,0,n,strongbool::yes))
					fft::SineTransform( &T(strip,1,0), T.Dim(3), T.Stride(3), 1 );
			}
		}
		Transpose(r,tw::grid::z,tw::grid::x,&T,-1);
	}
}

void Field::RealInverseAxialSineTransform(const Rng04& r,const DynSpace& ds)
{
	Field T;
	if (ds.GlobalDim(3)>1)
	{
		Transpose(r,tw::grid::z,tw::grid::x,&T,1);
		for (auto n=r.b0; n<r.e0; n++) {
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,3,0,n,strongbool::yes))
				{
					fft::SineTransform( &T(strip,1,0), T.Dim(3), T.Stride(3), -1 );
					T(strip,T.LNG(3),0) = -T(strip,2,0);
					T(strip,T.UNG(3),0) = 0.0;
				}
			}
		}
		Transpose(r,tw::grid::z,tw::grid::x,&T,-1);
	}
}

void Field::RealTransverseCosineTransform(const Rng04& r,const DynSpace& ds)
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=1;ax<=2;ax++)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (ds.GlobalDim(ax)>1)
		{
			Transpose(r,axis1,axis2,&T,1);
			for (auto n=r.b0; n<r.e0; n++) {
				#pragma omp parallel
				{
					for (auto strip : StripRange(T,ax,0,n,strongbool::yes))
						fft::CosineTransform( &T(strip,1,0), T.Dim(ax), T.Stride(ax), 1 );
				}
			}
			Transpose(r,axis1,axis2,&T,-1);
		}
	}
}

void Field::RealInverseTransverseCosineTransform(const Rng04& r,const DynSpace& ds)
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=2;ax>=1;ax--)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (ds.GlobalDim(ax)>1)
		{
			Transpose(r,axis1,axis2,&T,1);
			for (auto n=r.b0; n<r.e0; n++) {
				#pragma omp parallel
				{
					for (auto strip : StripRange(T,ax,0,n,strongbool::yes))
					{
						fft::CosineTransform( &T(strip,1,0), T.Dim(ax), T.Stride(ax), -1 );
						T(strip,T.LNG(ax),0) = T(strip,1,0);
						T(strip,T.UNG(ax),0) = T(strip,T.Dim(ax),0);
					}
				}
			}
			Transpose(r,axis1,axis2,&T,-1);
		}
	}
}

void Field::RealTransverseSineTransform(const Rng04& r,const DynSpace& ds)
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=1;ax<=2;ax++)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (ds.GlobalDim(ax)>1)
		{
			Transpose(r,axis1,axis2,&T,1);
			for (auto n=r.b0; n<r.e0; n++) {
				#pragma omp parallel
				{
					for (auto strip : StripRange(T,ax,0,n,strongbool::yes))
						fft::SineTransform( &T(strip,1,0), T.Dim(ax), T.Stride(ax), 1 );
				}
			}
			Transpose(r,axis1,axis2,&T,-1);
		}
	}
}

void Field::RealInverseTransverseSineTransform(const Rng04& r,const DynSpace& ds)
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=2;ax>=1;ax--)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (ds.GlobalDim(ax)>1)
		{
			Transpose(r,axis1,axis2,&T,1);
			for (auto n=r.b0; n<r.e0; n++) {
				#pragma omp parallel
				{
					for (auto strip : StripRange(T,ax,0,n,strongbool::yes))
					{
						fft::SineTransform( &T(strip,1,0), T.Dim(ax), T.Stride(ax), -1 );
						T(strip,T.LNG(ax),0) = -T(strip,2,0);
						T(strip,T.UNG(ax),0) = 0.0;
					}
				}
			}
			Transpose(r,axis1,axis2,&T,-1);
		}
	}
}

void Field::RealTransverseFFT(const Rng04& r,const DynSpace& ds)
{
	Field T;

	if (ds.GlobalDim(1)>1)
	{
		Transpose(r,tw::grid::x,tw::grid::z,&T,1);
		for (auto n=r.b0; n<r.e0; n++) {
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,1,0,n,strongbool::yes))
					fft::RealFFT( &T(strip,1,0), T.Dim(1), T.Stride(1), 1);
			}
		}
		Transpose(r,tw::grid::x,tw::grid::z,&T,-1);
	}

	if (ds.GlobalDim(2)>1)
	{
		Transpose(r,tw::grid::y,tw::grid::z,&T,1);
		for (auto n=r.b0; n<r.e0; n++) {
			const tw::Int xDim=T.Dim(1),zN0=T.LFG(3),zN1=T.UFG(3);
			#pragma omp parallel for collapse(2) schedule(static)
			for (tw::Int i=1;i<=xDim;i+=2) { // can't include ghost cells due to complex numbers; instead do copy ops below
				for (tw::Int k=zN0;k<=zN1;k++)
				{
					if (ds.GlobalCellIndex(i,1)==1)
					{
						fft::RealFFT( &T(n,i,1,k,0), T.Dim(2), T.Stride(2), 1);
						fft::RealFFT( &T(n,i+1,1,k,0), T.Dim(2), T.Stride(2), 1);
					}
					else
						fft::ComplexFFT( &T(n,i,1,k,0), &T(n,i+1,1,k,0), T.Dim(2), T.Stride(2), 1.0);
				}
			}
		}
		T.UpwardCopy(r,tw::grid::x,1);
		T.DownwardCopy(r,tw::grid::x,1);
		Transpose(r,tw::grid::y,tw::grid::z,&T,-1);
	}
}

void Field::RealInverseTransverseFFT(const Rng04& r,const DynSpace& ds)
{
	Field T;

	if (ds.GlobalDim(2)>1)
	{
		Transpose(r,tw::grid::y,tw::grid::z,&T,1);
		for (auto n=r.b0; n<r.e0; n++) {
			const tw::Int xDim=T.Dim(1),zN0=T.LFG(3),zN1=T.UFG(3);
			#pragma omp parallel for collapse(2) schedule(static)
			for (tw::Int i=1;i<=xDim;i+=2) { // can't include ghost cells due to complex numbers; instead do copy ops below
				for (tw::Int k=zN0;k<=zN1;k++)
				{
					if (ds.GlobalCellIndex(i,1)==1)
					{
						fft::RealFFT( &T(n,i,1,k,0), T.Dim(2), T.Stride(2), -1);
						T(n,i,T.LNG(2),k,0) = T(n,i,T.Dim(2),k,0);
						T(n,i,T.UNG(2),k,0) = T(n,i,1,k,0);
						fft::RealFFT( &T(n,i+1,1,k,0), T.Dim(2), T.Stride(2), -1);
						T(n,i+1,T.LNG(2),k,0) = T(n,i+1,T.Dim(2),k,0);
						T(n,i+1,T.UNG(2),k,0) = T(n,i+1,1,k,0);
					}
					else
					{
						fft::ComplexFFT( &T(n,i,1,k,0), &T(n,i+1,1,k,0), T.Dim(2), T.Stride(2), -1.0);
						T(n,i,T.LNG(2),k,0) = T(n,i,T.Dim(2),k,0);
						T(n,i,T.UNG(2),k,0) = T(n,i,1,k,0);
						T(n,i+1,T.LNG(2),k,0) = T(n,i+1,T.Dim(2),k,0);
						T(n,i+1,T.UNG(2),k,0) = T(n,i+1,1,k,0);
					}
				}
			}
		}
		T.UpwardCopy(r,tw::grid::x,1);
		T.DownwardCopy(r,tw::grid::x,1);
		Transpose(r,tw::grid::y,tw::grid::z,&T,-1);
	}

	if (ds.GlobalDim(1)>1)
	{
		Transpose(r,tw::grid::x,tw::grid::z,&T,1);
		for (auto n=r.b0; n<r.e0; n++) {
			#pragma omp parallel
			{
				for (auto strip : StripRange(T,1,0,n,strongbool::yes))
				{
					fft::RealFFT( &T(strip,1,0), T.Dim(1), T.Stride(1), -1);
					T(strip,T.LNG(1),0) = T(strip,T.Dim(1),0);
					T(strip,T.UNG(1),0) = T(strip,1,0);
				}
			}
		}
		Transpose(r,tw::grid::x,tw::grid::z,&T,-1);
	}
}

tw::Float Field::ComplexCyclicEigenvalue(tw::Int x,tw::Int y,const DynSpace& ds)
{
	x = ds.GlobalCellIndex(x,1);
	y = ds.GlobalCellIndex(y,2);
	return fft::eigenvalue_CFFT(x-1,ds.GlobalDim(1),freq[1]) + fft::eigenvalue_CFFT(y-1,ds.GlobalDim(2),freq[2]);
}

tw::Float Field::ComplexCyclicEigenvalue(tw::Int x,tw::Int y,tw::Int z,const DynSpace& ds)
{
	x = ds.GlobalCellIndex(x,1);
	y = ds.GlobalCellIndex(y,2);
	z = ds.GlobalCellIndex(z,3);
	return fft::eigenvalue_CFFT(x-1,ds.GlobalDim(1),freq[1]) + fft::eigenvalue_CFFT(y-1,ds.GlobalDim(2),freq[2]) + fft::eigenvalue_CFFT(z-1,ds.GlobalDim(3),freq[3]);
}

void Field::ComplexTransverseFFT(const Rng04& r,const DynSpace& ds)
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=1;ax<=2;ax++)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (ds.GlobalDim(ax)>1)
		{
			Transpose(r,axis1,axis2,&T,1);
			for (auto n=r.b0; n<r.e0; n++) {
				#pragma omp parallel
				{
					for (auto strip : StripRange(T,ax,0,n,strongbool::yes))
						fft::ComplexFFT( &T(strip,1,r.b4+0), &T(strip,1,r.b4+1), T.Dim(ax), T.Stride(ax), 1.0 );
				}
			}
			Transpose(r,axis1,axis2,&T,-1);
		}
	}
}

void Field::ComplexInverseTransverseFFT(const Rng04& r,const DynSpace& ds)
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=2;ax>=1;ax--)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::z;
		if (ds.GlobalDim(ax)>1)
		{
			Transpose(r,axis1,axis2,&T,1);
			for (auto n=r.b0; n<r.e0; n++) {
				#pragma omp parallel
				{
					for (auto strip : StripRange(T,ax,0,n,strongbool::yes))
					{
						fft::ComplexFFT( &T(strip,1,0), &T(strip,1,1), T.Dim(ax), T.Stride(ax), -1.0 );
						T(strip,T.LNG(ax),r.b4+0) = T(strip,T.Dim(ax),r.b4+0);
						T(strip,T.LNG(ax),r.b4+1) = T(strip,T.Dim(ax),r.b4+1);
						T(strip,T.UNG(ax),r.b4+0) = T(strip,1,r.b4+0);
						T(strip,T.UNG(ax),r.b4+1) = T(strip,1,r.b4+1);
					}
				}
			}
			Transpose(r,axis1,axis2,&T,-1);
		}
	}
}

void Field::ComplexFFT(const Rng04& r,const DynSpace& ds)
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=1;ax<=3;ax++)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::enumaxis(ax==3 ? 1 : ax+1);
		if (ds.GlobalDim(ax)>1)
		{
			Transpose(r,axis1,axis2,&T,1);
			for (auto n=r.b0; n<r.e0; n++) {
				#pragma omp parallel
				{
					for (auto strip : StripRange(T,ax,0,n,strongbool::yes))
						fft::ComplexFFT( &T(strip,1,r.b4+0), &T(strip,1,r.b4+1), T.Dim(ax), T.Stride(ax), 1.0 );
				}
			}
			Transpose(r,axis1,axis2,&T,-1);
		}
	}
}

void Field::ComplexInverseFFT(const Rng04& r,const DynSpace& ds)
{
	Field T;
	tw::grid::axis axis1,axis2;
	for (tw::Int ax=3;ax>=1;ax--)
	{
		axis1 = tw::grid::enumaxis(ax);
		axis2 = tw::grid::enumaxis(ax==3 ? 1 : ax+1);
		if (ds.GlobalDim(ax)>1)
		{
			Transpose(r,axis1,axis2,&T,1);
			for (auto n=r.b0; n<r.e0; n++) {
				#pragma omp parallel
				{
					for (auto strip : StripRange(T,ax,0,n,strongbool::yes))
					{
						fft::ComplexFFT( &T(strip,1,0), &T(strip,1,1), T.Dim(ax), T.Stride(ax), -1.0 );
						T(strip,T.LNG(ax),r.b4+0) = T(strip,T.Dim(ax),r.b4+0);
						T(strip,T.LNG(ax),r.b4+1) = T(strip,T.Dim(ax),r.b4+1);
						T(strip,T.UNG(ax),r.b4+0) = T(strip,1,r.b4+0);
						T(strip,T.UNG(ax),r.b4+1) = T(strip,1,r.b4+1);
					}
				}
			}
			Transpose(r,axis1,axis2,&T,-1);
		}
	}
}

void Field::Hankel(const Rng04& r,tw::Int modes,std::valarray<tw::Float>& matrix)
{
	Field T;
	Transpose(r,tw::grid::x,tw::grid::z,&T,1);
	for (auto n=r.b0; n<r.e0; n++) {
		#pragma omp parallel
		{
			for (auto strip : StripRange(T,1,0,n,strongbool::yes))
				for (auto c=r.b4; c<r.e4; c++)
					Transform(&T(strip,1,c),T.Dim(1),modes,T.Stride(1),matrix);
		}
	}
	Transpose(r,tw::grid::x,tw::grid::z,&T,-1);
}

void Field::InverseHankel(const Rng04& r,tw::Int modes,std::valarray<tw::Float>& matrix)
{
	Field T;
	Transpose(r,tw::grid::x,tw::grid::z,&T,1);
	for (auto n=r.b0; n<r.e0; n++) {
		#pragma omp parallel
		{
			for (auto strip : StripRange(T,1,0,n,strongbool::yes))
				for (auto c=r.b4; c<r.e4; c++)
					ReverseTransform(&T(strip,1,c),T.Dim(1),modes,T.Stride(1),matrix);
		}
	}
	Transpose(r,tw::grid::x,tw::grid::z,&T,-1);
}
