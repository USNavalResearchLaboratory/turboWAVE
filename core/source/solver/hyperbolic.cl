inline tw_Float diff2(__global tw_Float *f,tw_Int n,tw_Int s)
{
	return f[n+s] - 2*f[n] + f[n-s];
}

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

__kernel void YeeUpdateE(	__global tw_Float *A,
				__global tw_Float *J4,
				__global tw_Float *PMLx,
				__global tw_Float *PMLy,
				__global tw_Float *PMLz,
				__global tw_Float *met_g)
{
	// PROTOCOL : Local
	tw_Metrics met;
	GetMetrics((tw_Float*)&met,met_g);
	const tw_Int i = get_global_id(0);
	const tw_Int j = get_global_id(1);
	const tw_Int k = get_global_id(2);
	const tw_Int xs = get_global_offset(0)*12;
	const tw_Int ys = get_global_offset(1)*12*(get_global_size(0) + 2*get_global_offset(0));
	const tw_Int zs = get_global_offset(2)*12*(get_global_size(0) + 2*get_global_offset(0))*(get_global_size(1) + 2*get_global_offset(1));
	const tw_Int n = i*xs + j*ys + k*zs;
	const tw_Int m = n/3;

	A[n]   = PMLy[j*6]*A[n]   + PMLy[j*6+1]*(A[n+ys+10]+A[n+ys+11]-A[n+10]-A[n+11])/met.dy - 0.25*PMLy[j*6+1]*(J4[m-xs/3+1]+J4[m+1]);
	A[n+1] = PMLz[k*6]*A[n+1] - PMLz[k*6+1]*(A[n+zs+8]+A[n+zs+9]-A[n+8]-A[n+9])/met.dz - 0.25*PMLz[k*6+1]*(J4[m-xs/3+1]+J4[m+1]);
	A[n+2] = PMLz[k*6]*A[n+2] + PMLz[k*6+1]*(A[n+zs+6]+A[n+zs+7]-A[n+6]-A[n+7])/met.dz - 0.25*PMLz[k*6+1]*(J4[m-ys/3+2]+J4[m+2]);
	A[n+3] = PMLx[i*6]*A[n+3] - PMLx[i*6+1]*(A[n+xs+10]+A[n+xs+11]-A[n+10]-A[n+11])/met.dx - 0.25*PMLx[i*6+1]*(J4[m-ys/3+2]+J4[m+2]);
	A[n+4] = PMLx[i*6]*A[n+4] + PMLx[i*6+1]*(A[n+xs+8]+A[n+xs+9]-A[n+8]-A[n+9])/met.dx - 0.25*PMLx[i*6+1]*(J4[m-zs/3+3]+J4[m+3]);
	A[n+5] = PMLy[j*6]*A[n+5] - PMLy[j*6+1]*(A[n+ys+6]+A[n+ys+7]-A[n+6]-A[n+7])/met.dy - 0.25*PMLy[j*6+1]*(J4[m-zs/3+3]+J4[m+3]);
	// 66 FLOPs
}

__kernel void YeeUpdateB(	__global tw_Float *A,
				__global tw_Float *PMLx,
				__global tw_Float *PMLy,
				__global tw_Float *PMLz,
				__global tw_Float *met_g)
{
	tw_Metrics met;
	GetMetrics((tw_Float*)&met,met_g);
	// indexing scheme : [1...dim+1] , or [0] for dim=1
	// we will zero out the stride in directions where dim=1
	const tw_Int i = get_global_id(0);
	const tw_Int j = get_global_id(1);
	const tw_Int k = get_global_id(2);
	const tw_Int X = get_global_size(0)!=1;
	const tw_Int Y = get_global_size(1)!=1;
	const tw_Int Z = get_global_size(2)!=1;
	const tw_Int xs = X*12;
	const tw_Int ys = Y*12*(get_global_size(0)+X);
	const tw_Int zs = Z*12*(get_global_size(0)+X)*(get_global_size(1)+Y);
	const tw_Int n = i*xs + j*ys + k*zs;

	A[n+6] = PMLy[j*6+3]*A[n+6] - PMLy[j*6+4]*(A[n+4]+A[n+5]-A[n-ys+4]-A[n-ys+5])/met.dy;
	A[n+7] = PMLz[k*6+3]*A[n+7] + PMLz[k*6+4]*(A[n+2]+A[n+3]-A[n-zs+2]-A[n-zs+3])/met.dz;
	A[n+8] = PMLz[k*6+3]*A[n+8] - PMLz[k*6+4]*(A[n+0]+A[n+1]-A[n-zs+0]-A[n-zs+1])/met.dz;
	A[n+9] = PMLx[i*6+3]*A[n+9] + PMLx[i*6+4]*(A[n+4]+A[n+5]-A[n-xs+4]-A[n-xs+5])/met.dx;
	A[n+10] = PMLx[i*6+3]*A[n+10] - PMLx[i*6+4]*(A[n+2]+A[n+3]-A[n-xs+2]-A[n-xs+3])/met.dx;
	A[n+11] = PMLy[j*6+3]*A[n+11] + PMLy[j*6+4]*(A[n+0]+A[n+1]-A[n-ys+0]-A[n-ys+1])/met.dy;
	// 42 FLOPs
}

__kernel void PrepCenteredFields(__global tw_Float *F,__global tw_Float *A)
{
	// indexing scheme : [0...dim] , or [0] for dim=1
	// we will zero out the stride in directions where dim=1
	const tw_Int i = get_global_id(0);
	const tw_Int j = get_global_id(1);
	const tw_Int k = get_global_id(2);
	const tw_Int X = get_global_size(0)!=1;
	const tw_Int Y = get_global_size(1)!=1;
	const tw_Int Z = get_global_size(2)!=1;
	const tw_Int xs = X*12;
	const tw_Int ys = Y*12*(get_global_size(0)+X);
	const tw_Int zs = Z*12*(get_global_size(0)+X)*(get_global_size(1)+Y);
	const tw_Int n = i*xs + j*ys + k*zs;
	const tw_Int m = n/2;

	F[m+3] = 0.125*(A[n+6] + A[n+7] + A[n+ys+6] + A[n+ys+7] + A[n+zs+6] + A[n+zs+7] + A[n+ys+zs+6] + A[n+ys+zs+7]);
	F[m+4] = 0.125*(A[n+8] + A[n+9] + A[n+xs+8] + A[n+xs+9] + A[n+zs+8] + A[n+zs+9] + A[n+xs+zs+8] + A[n+xs+zs+9]);
	F[m+5] = 0.125*(A[n+10] + A[n+11] + A[n+xs+10] + A[n+xs+11] + A[n+ys+10] + A[n+ys+11] + A[n+xs+ys+10] + A[n+xs+ys+11]);
	// 24 FLOPs
}

__kernel void CenteredFields(__global tw_Float *F,__global tw_Float *A)
{
	// indexing scheme : [0...dim] , or [0] for dim=1
	// we will zero out the stride in directions where dim=1
	const tw_Int i = get_global_id(0);
	const tw_Int j = get_global_id(1);
	const tw_Int k = get_global_id(2);
	const tw_Int X = get_global_size(0)!=1;
	const tw_Int Y = get_global_size(1)!=1;
	const tw_Int Z = get_global_size(2)!=1;
	const tw_Int xs = X*12;
	const tw_Int ys = Y*12*(get_global_size(0)+X);
	const tw_Int zs = Z*12*(get_global_size(0)+X)*(get_global_size(1)+Y);
	const tw_Int n = i*xs + j*ys + k*zs;
	const tw_Int m = n/2;

	F[m+0] = 0.5*(A[n+0] + A[n+1] + A[n+xs+0] + A[n+xs+1]);
	F[m+1] = 0.5*(A[n+2] + A[n+3] + A[n+ys+2] + A[n+ys+3]);
	F[m+2] = 0.5*(A[n+4] + A[n+5] + A[n+zs+4] + A[n+zs+5]);

	F[m+3] += 0.125*(A[n+6] + A[n+7] + A[n+ys+6] + A[n+ys+7] + A[n+zs+6] + A[n+zs+7] + A[n+ys+zs+6] + A[n+ys+zs+7]);
	F[m+4] += 0.125*(A[n+8] + A[n+9] + A[n+xs+8] + A[n+xs+9] + A[n+zs+8] + A[n+zs+9] + A[n+xs+zs+8] + A[n+xs+zs+9]);
	F[m+5] += 0.125*(A[n+10] + A[n+11] + A[n+xs+10] + A[n+xs+11] + A[n+ys+10] + A[n+ys+11] + A[n+xs+ys+10] + A[n+xs+ys+11]);
	// 36 FLOPs
}

__kernel void LorentzAdvance(
	__global tw_Float *A4,
	__global tw_Float *Ao4,
	__global tw_Float *J4,
	__global tw_Float *met_g,
	const tw_Float mult)
{
	LOCAL_PROTOCOL_MACRO

	// Do not update scalar potential for now
	// Reverses role of old and new data in order to perform in place

	const tw_Float dt = met.dt;
	const tw_Float dth = 0.5*met.dt;
	const tw_Float dxi = 1/met.dx;
	const tw_Float dyi = 1/met.dy;
	const tw_Float dzi = 1/met.dz;
	const tw_Int n1 = n + 1*cs;
	const tw_Int n2 = n + 2*cs;
	const tw_Int n3 = n + 3*cs;

	//Ao4[n] = 2*A4[n] - Ao4[n] + dt*dt*mult*J4[n];
	//Ao4[n] += dt*dt*(dxi*dxi*diff2(A4,n,xs) + dyi*dyi*diff2(A4,n,ys) + dzi*dzi*diff2(A4,n,zs);

	Ao4[n1] = 2*A4[n1] - Ao4[n1] + dt*dt*mult*J4[n1];
	Ao4[n1] += dt*dt*(dxi*dxi*diff2(A4,n1,xs) + dyi*dyi*diff2(A4,n1,ys) + dzi*dzi*diff2(A4,n1,zs));

	Ao4[n2] = 2*A4[n2] - Ao4[n2] + dt*dt*mult*J4[n2];
	Ao4[n2] += dt*dt*(dxi*dxi*diff2(A4,n2,xs) + dyi*dyi*diff2(A4,n2,ys) + dzi*dzi*diff2(A4,n2,zs));

	Ao4[n3] = 2*A4[n3] - Ao4[n3] + dt*dt*mult*J4[n3];
	Ao4[n3] += dt*dt*(dxi*dxi*diff2(A4,n3,xs) + dyi*dyi*diff2(A4,n3,ys) + dzi*dzi*diff2(A4,n3,zs));
}

__kernel void LorentzSwap(
	__global tw_Float *A4,
	__global tw_Float *Ao4)
{
	POINT_PROTOCOL_MACRO

	const tw_Int n1 = n + 1*cs;
	const tw_Int n2 = n + 2*cs;
	const tw_Int n3 = n + 3*cs;
	tw_Float temp;

	temp = A4[n];
	A4[n] = Ao4[n];
	Ao4[n] = temp;

	temp = A4[n1];
	A4[n1] = Ao4[n1];
	Ao4[n1] = temp;

	temp = A4[n2];
	A4[n2] = Ao4[n2];
	Ao4[n2] = temp;

	temp = A4[n3];
	A4[n3] = Ao4[n3];
	Ao4[n3] = temp;
}

__kernel void LorentzMidstep(
	__global tw_Float *A4,
	__global tw_Float *Ao4)
{
	POINT_PROTOCOL_MACRO

	const tw_Int n1 = n + 1*cs;
	const tw_Int n2 = n + 2*cs;
	const tw_Int n3 = n + 3*cs;

	A4[n] = 0.5*(Ao4[n] + A4[n]);
	A4[n1] = 0.5*(Ao4[n1] + A4[n1]);
	A4[n2] = 0.5*(Ao4[n2] + A4[n2]);
	A4[n3] = 0.5*(Ao4[n3] + A4[n3]);
}

__kernel void LorentzUndoMidstep(
	__global tw_Float *A4,
	__global tw_Float *Ao4)
{
	POINT_PROTOCOL_MACRO

	const tw_Int n1 = n + 1*cs;
	const tw_Int n2 = n + 2*cs;
	const tw_Int n3 = n + 3*cs;

	A4[n] = 2.0*A4[n] - Ao4[n];
	A4[n1] = 2.0*A4[n1] - Ao4[n1];
	A4[n2] = 2.0*A4[n2] - Ao4[n2];
	A4[n3] = 2.0*A4[n3] - Ao4[n3];
}
