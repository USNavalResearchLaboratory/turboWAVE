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

inline void LaplacianCoefficients(tw_Float *D,const tw_cell_id x,tw_Metrics *m)
{
	// Set D[1] through D[6], D[0] can be used externally for the sum if desired
	const tw_Float r = (tw_Float)x.s1;
	const tw_Float r1 = m->car + m->cyl*(m->xo + (r-1.0)*m->dx);
	const tw_Float rc = m->car + m->cyl*(m->xo + (r-0.5)*m->dx);
	const tw_Float r2 = m->car + m->cyl*(m->xo + (r-0.0)*m->dx);
	const tw_Float dV = m->dx*m->dy*m->dz*rc;
	D[1] = m->dy*m->dz*r1 / (dV*m->dx);
	D[2] = m->dy*m->dz*r2 / (dV*m->dx);
	D[3] = D[4] = m->dx*m->dz / (dV*m->dy*rc);
	D[5] = D[6] = m->dx*m->dy*rc / (dV*m->dz);
}

__kernel void SORIterationPoisson(	__global tw_Float *phi,
				__global tw_Float *source,
				__global tw_Float *residual,
				__global char *mask,
				const tw_Float mul,
				const tw_Float overrelaxation,
				__global tw_Float *met_g)
{
	tw_Metrics met;
	tw_Float res,f00;
	tw_Float D[7];
	
	GetMetrics((tw_Float*)&met,met_g);
	const tw_Int xDim = get_global_size(0);
	const tw_Int yDim = get_global_size(1);
	const tw_Int xo = get_global_offset(0);
	const tw_Int yo = get_global_offset(1);
	const tw_Int zo = get_global_offset(2);
	const tw_Int xs = xo;
	const tw_Int ys = yo*(xDim + 2*xo);
	const tw_Int zs = zo*(xDim + 2*xo)*(yDim + 2*yo);
	const tw_Int i = get_global_id(0) + 1 - xo;
	const tw_Int j = get_global_id(1) + 1 - yo;
	const tw_Int k = get_global_id(2) + 1 - zo;
	const tw_Int n = i*xs + j*ys + k*zs;

	// unconditionally clear the residual so L1Norm works properly
	residual[n] = 0.0;

	if (mask[(i-1) + (j-1)*xDim + (k-1)*xDim*yDim]==0)
		return;
	
	f00 = phi[n];
	
	LaplacianCoefficients(D,(tw_cell_id)(0,i,j,k),&met);
	D[1] *= (tw_Float)xo;
	D[2] *= (tw_Float)xo;
	D[3] *= (tw_Float)yo;
	D[4] *= (tw_Float)yo;
	D[5] *= (tw_Float)zo;
	D[6] *= (tw_Float)zo;
	D[0] = -(D[1] + D[2] + D[3] + D[4] + D[5] + D[6]);
	res = D[1]*(phi[n-xs]-f00);
	res += D[2]*(phi[n+xs]-f00);
	res += D[3]*(phi[n-ys]-f00);
	res += D[4]*(phi[n+ys]-f00);
	res += D[5]*(phi[n-zs]-f00);
	res += D[6]*(phi[n+zs]-f00);
	res -= mul*source[n];
	phi[n] -= overrelaxation*res/D[0];
	residual[n] = res;
}

__kernel void SORIterationGeneral(	__global tw_Float *phi,
									__global tw_Float *source,
									__global tw_Float *residual,
									__global char *mask,
									const tw_Float mul,
									const tw_Float overrelaxation,
									__global tw_Float *met_g,
									__global tw_Float *coeff)
{
	tw_Metrics met;
	tw_Float res,f00,c00;
	tw_Float D[7];
	
	GetMetrics((tw_Float*)&met,met_g);
	const tw_Int xDim = get_global_size(0);
	const tw_Int yDim = get_global_size(1);
	const tw_Int xo = get_global_offset(0);
	const tw_Int yo = get_global_offset(1);
	const tw_Int zo = get_global_offset(2);
	const tw_Int xs = xo;
	const tw_Int ys = yo*(xDim + 2*xo);
	const tw_Int zs = zo*(xDim + 2*xo)*(yDim + 2*yo);
	const tw_Int i = get_global_id(0) + 1 - xo;
	const tw_Int j = get_global_id(1) + 1 - yo;
	const tw_Int k = get_global_id(2) + 1 - zo;
	const tw_Int n = i*xs + j*ys + k*zs;

	// unconditionally clear the residual so L1Norm works properly
	residual[n] = 0.0;
	
	if (mask[(i-1) + (j-1)*xDim + (k-1)*xDim*yDim]==0)
		return;
	
	f00 = phi[n];
	c00 = coeff[n];
	
	LaplacianCoefficients(D,(tw_cell_id)(0,i,j,k),&met);
	D[1] *= (tw_Float)xo*0.5*(coeff[n-xs]+c00);
	D[2] *= (tw_Float)xo*0.5*(coeff[n+xs]+c00);
	D[3] *= (tw_Float)yo*0.5*(coeff[n-ys]+c00);
	D[4] *= (tw_Float)yo*0.5*(coeff[n+ys]+c00);
	D[5] *= (tw_Float)zo*0.5*(coeff[n-zs]+c00);
	D[6] *= (tw_Float)zo*0.5*(coeff[n+zs]+c00);
	D[0] = -(D[1] + D[2] + D[3] + D[4] + D[5] + D[6]);
	res = D[1]*(phi[n-xs]-f00);
	res += D[2]*(phi[n+xs]-f00);
	res += D[3]*(phi[n-ys]-f00);
	res += D[4]*(phi[n+ys]-f00);
	res += D[5]*(phi[n-zs]-f00);
	res += D[6]*(phi[n+zs]-f00);
	res -= mul*source[n];
	
	phi[n] -= overrelaxation*res/D[0];
	residual[n] = res;
}
