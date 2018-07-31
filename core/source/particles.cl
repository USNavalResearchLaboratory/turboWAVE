// important: operations using index map must be atomic

struct Primitive
{
	float x[3];	// x[0],x[1],x[2] = relative position in cell using interval [-0.5,0.5); linearly mapped to curvilinear coordinates in cell
	tw_Int cell; // index in memory of the reference cell
};

struct Particle
{
	Primitive q;
	tw_Float p[3];
	float n,aux1,aux2;
};

inline void RefCellIndices(Primitive *q,tw_Int *rc,tw_DiscreteSpace *ds)
{
	tw_Int X = ds->xN1 - ds->xN0 + 1;
	tw_Int Y = ds->yN1 - ds->yN0 + 1;
	tw_Int Z = ds->zN1 - ds->zN0 + 1;
	rc[2] = q.cell / (X*Y);
	rc[1] = (q.cell - rc[2]*X*Y) / X;
	rc[0] = q.cell - rc[1]*X - rc[2]*X*Y;
}

inline bool PrimitiveInDomain(Primitive *q,tw_DiscreteSpace *ds,tw_Metrics *m)
{
	tw_Int rc[3];
	tw_Float r[3];
	RefCellIndices(q,rc,ds);
	r[0] = (tw_Float)q.x[0] + (tw_Float)rc[0] - 0.5;
	r[1] = (tw_Float)q.x[1] + (tw_Float)rc[1] - 0.5;
	r[2] = (tw_Float)q.x[2] + (tw_Float)rc[2] - 0.5;
	return	((r[0]>=0.0 && r[0]<(tw_Float)ds->xDim) || ds->xDim==1) &&
			((r[1]>=0.0 && r[1]<(tw_Float)ds->yDim) || ds->yDim==1) &&
			((r[2]>=0.0 && r[2]<(tw_Float)ds->zDim) || ds->zDim==1);
}

inline bool PrimitiveInCell(Primitive *q)
{
	return q->x[0]>=-0.5 && q->x[0]<0.5 && q->x[1]>=-0.5 && q->x[1]<0.5 && q->x[2]>=-0.5 && q->x[2]<0.5;
}

inline tw_Float vi_new_basis(tw_Int i,__global tw_Float *mat,__global tw_Float *vec)
{
	// basis vectors are layed out in mat as ux,uy,uz,vx,vy,vz,wx,vy,wz
	return mat[i*3+0]*vec[0] + mat[i*3+1]*vec[1] + mat[i*3+2]*vec[2];
}

inline tw_Float vi_std_basis(tw_Int i,__global tw_Float *mat,tw_Float *vec)
{
	// basis vectors are layed out in mat as ux,uy,uz,vx,vy,vz,wx,vy,wz
	return mat[i+0]*vec[0] + mat[i+3]*vec[1] + mat[i+6]*vec[2];
}

inline void GetWeights(float *w,float *x)
{
	w[0] = 0.125f - 0.5f*x[0] + 0.5f*x[0]*x[0];
	w[1] = 0.125f - 0.5f*x[1] + 0.5f*x[1]*x[1];
	w[2] = 0.125f - 0.5f*x[2] + 0.5f*x[2]*x[2];
	w[3] = 0.75f - x[0]*x[0];
	w[4] = 0.75f - x[1]*x[1];
	w[5] = 0.75f - x[2]*x[2];
	w[6] = 0.125f + 0.5f*x[0] + 0.5f*x[0]*x[0];
	w[7] = 0.125f + 0.5f*x[1] + 0.5f*x[1]*x[1];
	w[8] = 0.125f + 0.5f*x[2] + 0.5f*x[2]*x[2];
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

inline tw_Float CellVolume(const tw_cell_id c,tw_Metrics *m)
{
	const tw_Float x0 = m->xo + ((tw_Float)c.s1-0.5)*m->dx;
	const tw_Float y0 = m->yo + ((tw_Float)c.s2-0.5)*m->dy;
	const tw_Float rc = m->car + m->cyl*x0 + m->par*0.25*(x0+y0);
	return m->dx*m->dy*m->dz*rc;
}

__kernel void PushBound(__global tw_Float *dens,
						__global tw_Float *R0,
						__global tw_Float *R1,
						__global tw_Float *EM,
						__global tw_Float *J4,
						__global tw_Float *packet,
						__global tw_Float *met_g,
						const tw_Float q0,const tw_Float m0,const tw_Float dt)
{
	// PROTOCOL : Local
	// packet : resFreq , dampFreq , oscStrength , a36 , b , d , T33
	// T33 layout : ux,uy,uz,vx,vy,vz,wx,wy,wz
	// a36 layout : a11...a16,a21...a26,a31...a36
	// Solve Dt^2(R) + 2*nu*Dt(R) + omega^2*R = FNL - E
	// FNL = -a*R*R + b*(R.R)*R (really F/m)
	// At start, E is known at t = 0, R0 at t = -1, R1 at t = 0
	
	tw_Metrics met;
	GetMetrics((tw_Float*)&met,met_g);
	
	const tw_Int ig = get_global_id(0);
	const tw_Int jg = get_global_id(1);
	const tw_Int kg = get_global_id(2);
	const tw_Int xs = get_global_offset(0);
	const tw_Int ys = get_global_offset(1)*(get_global_size(0) + 2*get_global_offset(0));
	const tw_Int zs = get_global_offset(2)*(get_global_size(0) + 2*get_global_offset(0))*(get_global_size(1) + 2*get_global_offset(1));
	const tw_Int n = ig*xs + jg*ys + kg*zs;
	
	tw_Float dV = CellVolume((tw_cell_id)(0,ig,jg,kg),&met);
	
	tw_Int i;
	tw_Float r2,FNL[3],r0[3],r1[3];

	for (i=0;i<3;i++)
	{
		r0[i] = R0[n*3+i];
		r1[i] = R1[n*3+i];
	}
	r2 = r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2];
	for (i=0;i<3;i++)
	{
		FNL[i] = -(packet[9+i*6+0]*r1[0]*r1[0] + packet[9+i*6+1]*r1[1]*r1[1] + packet[9+i*6+2]*r1[2]*r1[2]);
		FNL[i] -= 2.0*(packet[9+i*6+3]*r1[1]*r1[2] + packet[9+i*6+4]*r1[0]*r1[2] + packet[9+i*6+5]*r1[0]*r1[1]);
		FNL[i] += packet[27]*r2*r1[i] + packet[28]*r2*r2*r1[i];
	}
	for (i=0;i<3;i++)
	{
		r1[i] = r1[i]*(2.0 + 2.0*packet[3+i]*dt - packet[0+i]*packet[0+i]*dt*dt) - r0[i] + FNL[i]*dt*dt;
		r1[i] += (q0/m0) * vi_new_basis(i,&packet[29],&EM[n*6]) * dt * dt;
		r1[i] /= 1.0 + 2.0*packet[3+i]*dt;
		r0[i] = R0[n*3+i] = R1[n*3+i];
		R1[n*3+i] = r1[i];
		r0[i] *= packet[6+i];
		r1[i] *= packet[6+i];
		r0[i] = (r1[i] - r0[i])/dt; // put crystal basis velocity in r0
	}
	for (i=0;i<3;i++)
		J4[n*4+i+1] += q0*dens[n]*vi_std_basis(i,&packet[29],r0)*dV;
}

__kernel void Push(	__global Particle *par,
					__global tw_Float *EM,
					__global tw_Float *J4,
					__global tw_Float *met_g,
					const tw_Float q0,
					const tw_Float m0,
					tw_cell_id stride,
					const tw_Float dt)
{
	tw_Int i,j,k,n;
	Primitive q;
	tw_Float wijk,dens;
	float w[9],xh[3],vel[3];
	tw_vec4 E,B,p,t,s,impulse;
	tw_Metrics met;
	const tw_Int ig = get_global_id(0);
	GetMetrics((tw_Float*)&met,met_g);
	q = par[ig].q;
	p = (tw_vec4)(par[ig].p[0],par[ig].p[1],par[ig].p[2],0);
	dens = par[ig].n;
	
	// tw::vec4 signature during momentum update is +++-
	
	GetWeights(w,q.x);
	E = (tw_vec4)(0,0,0,0);
	B = (tw_vec4)(0,0,0,0);
	for (k=-1;k<=1;k++)
		for (j=-1;j<=1;j++)
			for (i=-1;i<=1;i++)
			{
				wijk = w[i*3+3]*w[j*3+4]*w[k*3+5];
				n = q.cell + i*stride[1] + j*stride[2] + k*stride[3];
				E.x += EM[n*6+0]*wijk;
				E.y += EM[n*6+1]*wijk;
				E.z += EM[n*6+2]*wijk;
				B.x += EM[n*6+3]*wijk;
				B.y += EM[n*6+4]*wijk;
				B.z += EM[n*6+5]*wijk;
			}
	
	impulse = 0.5*q0*E*dt;
	p += impulse;
	t = 0.5*q0*B*dt/sqrt(m0*m0+dot(p,p));
	s = 2.0*t/(1.0 + dot(t,t));
	p += cross(p + cross(p,t),s);
	p += impulse;
	
	wijk = sqrt(m0*m0 + dot(p,p));
	vel[0] = p.x/wijk;
	vel[1] = p.y/wijk;
	vel[2] = p.z/wijk;
	
	xh[0] = q.x[0] + 0.5*vel[0]*dt/met.dx;
	xh[1] = q.x[1] + 0.5*vel[1]*dt/met.dy;
	xh[2] = q.x[2] + 0.5*vel[2]*dt/met.dz;
	q.x[0] += vel[0]*dt/met.dx;
	q.x[1] += vel[1]*dt/met.dy;
	q.x[2] += vel[2]*dt/met.dz;
	
	// tw::vec4 signature during deposition is -+++
	
	s = q.n*q0*tw_vec4(1,vel.x,vel.y,vel.z);
	
	if (PrimitiveInCell(&q))
		goto deposit;
	return;
	
	deposit:
	
	GetWeights(w,xh);
	for (k=-1;k<=1;k++)
		for (j=-1;j<=1;j++)
			for (i=-1;i<=1;i++)
			{
				wijk = w[i*3+3]*w[j*3+4]*w[k*3+5];
				n = q.cell + i*stride[1] + j*stride[2] + k*stride[3];
				J4[n*4+1] += s.s1*wijk;
				J4[n*4+2] += s.s2*wijk;
				J4[n*4+3] += s.s3*wijk;
			}
	
	GetWeights(w,q.x);
	for (k=-1;k<=1;k++)
		for (j=-1;j<=1;j++)
			for (i=-1;i<=1;i++)
			{
				wijk = w[i*3+3]*w[j*3+4]*w[k*3+5];
				n = q.cell + i*stride[1] + j*stride[2] + k*stride[3];
				J4[n*4+0] += s.s0*wijk;
			}
}
