// important: operations using index map must be atomic.
// x-packing is assumed for real buffers.
// c-packing is assumed for complex buffers.

inline void GetMetrics(tw_Float *dest,__global tw_Float *src)
{
	dest[0] = src[0];
	dest[1] = src[1];
	dest[2] = src[2];
	dest[3] = src[3];
	dest[4] = src[4];
	dest[5] = src[5];
	dest[6] = src[6];
	dest[7] = src[7];
	dest[8] = src[8];
	dest[9] = src[9];
	dest[10] = src[10];
	dest[11] = src[11];
}

inline tw_Float CellVolume(const tw_cell_id c,tw_Metrics *m)
{
	// Cartesian or Cylindrical geometry
	const tw_Float x0 = m->xo + ((tw_Float)c.s1-0.5f)*m->dx;
	const tw_Float rc = m->car + m->cyl*x0;
	return m->dx*m->dy*m->dz*rc;
}

__kernel void FieldToBoundaryData(__global tw_Float *F, __global tw_Float *B, __global tw_Int *map)
{
	// PROTOCOL : Special
	B[get_global_id(0)] = F[map[get_global_id(0)]];
}

__kernel void BoundaryDataToField(__global tw_Float *F, __global tw_Float *B, __global tw_Int *map)
{
	// PROTOCOL : Special
	F[map[get_global_id(0)]] = B[get_global_id(0)];
}

__kernel void GhostDataToField(__global tw_Float *F, __global tw_Float *B, __global tw_Int *map)
{
	// PROTOCOL : Special
	F[map[get_global_id(0)]] = B[get_global_id(0)];
}

__kernel void ZeroGhostCells(__global tw_Float *F, __global tw_Int *map)
{
	// PROTOCOL : Special
	F[map[get_global_id(0)]] = 0.0f;
}

__kernel void FillVec4Field(__global tw_Float *A4,const tw_Complex v0,const tw_Complex v1,const tw_Complex v2,const tw_Complex v3)
{
	// PROTOCOL : Cell
	// newVal = oldVal*imag(arg) + real(arg)
	const tw_Int i = get_global_id(0);
	const tw_Int n = get_global_size(0);
	A4[n*0+i] = v0.s1*A4[n*0+i] + v0.s0;
	A4[n*1+i] = v1.s1*A4[n*1+i] + v1.s0;
	A4[n*2+i] = v2.s1*A4[n*2+i] + v2.s0;
	A4[n*3+i] = v3.s1*A4[n*3+i] + v3.s0;
}

__kernel void SwapBuffers( __global tw_Float *buff1 , __global tw_Float *buff2 )
{
	// PROTOCOL : Element
	const tw_Int n = get_global_id(0);
	tw_Float temp = buff1[n];
	buff1[n] = buff2[n];
	buff2[n] = temp;
}

__kernel void DestructiveSum(__global tw_Float *data)
{
	// PROTOCOL : Special
	// perform in place (destroys data), assumes ghost cells have been zeroed
	// has to be called repeatedly
	// result winds up in ghost cell (0,0,0)
	const tw_Int ig = get_global_id(0);
	const tw_Int Ni = get_global_size(0);
	if (ig<Ni/2)
	{
		data[ig] += data[ig+Ni/2];
		if (ig==0 && Ni%2==1)
			data[ig] += data[Ni-1];
	}
}

__kernel void DestructiveNorm1(__global tw_Float *data)
{
	// PROTOCOL : Special
	// Same as "Sum" but uses absolute value
	const tw_Int ig = get_global_id(0);
	const tw_Int Ni = get_global_size(0);
	if (ig<Ni/2)
	{
		data[ig] = fabs(data[ig]) + fabs(data[ig+Ni/2]);
		if (ig==0 && Ni%2==1)
			data[ig] += fabs(data[Ni-1]);
	}
}

__kernel void WeightByVolume(__global tw_Float *data,__global tw_Float *met_g,const tw_Int comp,const tw_Float inv)
{
	LOCAL_PROTOCOL_MACRO
	tw_Float dV = CellVolume(cell,&met);
	data[n + cs*comp] *= (1.0f-inv)*dV + inv/dV;
}

__kernel void ComplexMod2(__global tw_Float *dst,__global tw_Complex *src)
{
	// PROTOCOL : Element (must be with respect to dst)
	const tw_Int ig = get_global_id(0);
	dst[ig] = src[ig].s0*src[ig].s0 + src[ig].s1*src[ig].s1;
}

__kernel void DestructiveComplexMod2(__global tw_Complex *data)
{
	// PROTOCOL : Cell
	const tw_Int ig = get_global_id(0);
	data[ig].s0 = data[ig].s0*data[ig].s0 + data[ig].s1*data[ig].s1;
	data[ig].s1 = 0.0f;
}

__kernel void MADD(__global tw_Float *data,const tw_Float m,const tw_Float a)
{
	// PROTOCOL : Element
	const tw_Int ig = get_global_id(0);
	data[ig] = m*data[ig] + a;
}
