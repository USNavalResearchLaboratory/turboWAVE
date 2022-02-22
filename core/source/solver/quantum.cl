void GetUnitaryFactor(
	const tw_Int n,
	__global tw_Float *A4,
	const tw_Float partitionFactor,
	const tw_cell_id cell,
	tw_Metrics *met,
	tw_Strip *strip,
	tw_Complex *H1,
	tw_Complex *H2,
	tw_Complex *H3,
	tw_Complex dt,
	tw_Float sgn);

inline void GetWave(
	__global tw_Float *waves,
	tw_Int wave_num,
	tw_Float *a0,
	tw_Float *w0,
	tw_Float *u,
	tw_Float *w,
	tw_Float *f,
	tw_Float *tp)
{
	tw_Int n = wave_num*15;
	*a0 = waves[n+0];
	*w0 = waves[n+1];
	for (tw_Int i=0;i<3;i++)
	{
		u[i] = waves[n+2+i];
		w[i] = waves[n+5+i];
		f[i] = waves[n+8+i];
	}
	for (tw_Int i=0;i<4;i++)
		tp[i] = waves[n+11+i];
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

inline void GetStrip(tw_Int *dest,__global tw_Int *src)
{
	dest[0] = src[0];
	dest[1] = src[1];
	dest[2] = src[2];
	dest[3] = src[3];
	dest[4] = src[4];
	dest[5] = src[5];
}

inline void LaplacianParameters(const tw_cell_id cell,tw_Metrics *m,tw_Strip *strip,tw_Float *D1,tw_Float *D2,tw_Float *l1,tw_Float *l2)
{
	const tw_Float r = (tw_Float)cell.s1;
	const tw_Float rw1 = m->car + m->cyl*(m->xo + (r-1.0f)*m->dx);
	const tw_Float rc = m->car + m->cyl*(m->xo + (r-0.5f)*m->dx);
	const tw_Float rw2 = m->car + m->cyl*(m->xo + (r-0.0f)*m->dx);
	const tw_Float dV = m->dx*m->dy*m->dz*rc;
	*D1 = strip->di*m->dy*m->dz*rw1 + strip->dj*m->dx*m->dz + strip->dk*m->dx*m->dy*rc;
	*D2 = strip->di*m->dy*m->dz*rw2 + strip->dj*m->dx*m->dz + strip->dk*m->dx*m->dy*rc;
	*l1 = strip->di*m->dx + strip->dj*m->dy*rc + strip->dk*m->dz;
	*l2 = strip->di*m->dx + strip->dj*m->dy*rc + strip->dk*m->dz;
	*D1 /= (*l1)*dV;
	*D2 /= (*l2)*dV;
}

inline tw_Float QuinticRise(tw_Float x)
{
	return 10.0f*x*x*x - 15.0f*x*x*x*x + 6.0f*x*x*x*x*x;
}

inline tw_Float PulseShapeFactor(tw_Float t,tw_Float *tp)
{
	const tw_Float hold = (tw_Float)(t > tp[1] && t <= tp[2]);
	const tw_Float tau_rise = (tw_Float)(t > tp[0] && t <= tp[1]) * (t-tp[0]) / (tp[1] - tp[0]);
	const tw_Float tau_fall = (tw_Float)(t > tp[2] && t <= tp[3]) * (1.0 - (t-tp[2]) / (tp[3] - tp[2]));
	return QuinticRise(tau_rise + tau_fall + hold);
}

inline tw_Complex conj(tw_Complex z)
{
	return (tw_Complex)(z.s0,-z.s1);
}

inline tw_Complex cmul(tw_Complex a,tw_Complex b)
{
	return (tw_Complex)(a.s0*b.s0 - a.s1*b.s1 , a.s0*b.s1 + a.s1*b.s0);
}

inline tw_Complex cdiv(tw_Complex a,tw_Complex b)
{
	return (tw_Complex)(a.s0*b.s0 + a.s1*b.s1 , a.s1*b.s0 - a.s0*b.s1)/(b.s0*b.s0 + b.s1*b.s1);
}

inline tw_Float sqr(tw_Float x)
{
	return x*x;
}

inline tw_Float diff(__global tw_Float *f,tw_Int n,tw_Int s)
{
	return f[n+s] - f[n-s];
}

inline tw_Float diff2(__global tw_Float *f,tw_Int n,tw_Int s)
{
	return f[n+s] - 2*f[n] + f[n-s];
}

inline tw_Float sfwd(__global tw_Float *f,tw_Int n,tw_Int s)
{
	return 0.5*(f[n] + f[n+s]);
}

inline tw_Float sbak(__global tw_Float *f,tw_Int n,tw_Int s)
{
	return 0.5*(f[n-s] + f[n]);
}

inline void split_step_advance(tw_Float *Zr,tw_Float *Zi,const tw_Float Ur,const tw_Float Ui,const tw_Float Drh,const tw_Float Dih)
{
	*Zr += Drh;
	*Zi += Dih;
	tw_Float temp = *Zr;
	*Zr = *Zr*Ur - *Zi*Ui;
	*Zi = temp*Ui + *Zi*Ur;
	*Zr += Drh;
	*Zi += Dih;
}

void GetUnitaryFactor(
	const tw_Int n,
	__global tw_Float *A4,
	const tw_Float partitionFactor,
	const tw_cell_id cell,
	tw_Metrics *met,
	tw_Strip *strip,
	tw_Complex *U1,
	tw_Complex *U2,
	tw_Complex *U3,
	const tw_Complex dt,
	const tw_Float sgn)
{
	const tw_Complex iihalfsgn = (tw_Complex)(0.0f,sgn*0.5f);
	tw_Float l1,l2,D1,D2,A1,A2,Ueff;

	LaplacianParameters(cell,met,strip,&D1,&D2,&l1,&l2);
	A1 = 0.5f*( A4[n*4+strip->axis] + A4[(n-strip->stride)*4+strip->axis] );
	A2 = 0.5f*( A4[n*4+strip->axis] + A4[(n+strip->stride)*4+strip->axis] );
	Ueff = partitionFactor*(0.5f*(A4[n*4+1]*A4[n*4+1] + A4[n*4+2]*A4[n*4+2] + A4[n*4+3]*A4[n*4+3]) - A4[n*4]);
	*U1 = (tw_Complex)( -0.5f*D1 , 0.5f*A1*D1*l1 ); // H1
	*U2 = (tw_Complex)( 0.5f*(D1+D2) + Ueff , 0.0f ); // H2
	*U3 = (tw_Complex)( -0.5f*D2 , -0.5f*A2*D2*l2 ); // H3
	*U1 = cmul(dt,cmul(iihalfsgn,*U1));
	*U2 = (tw_Complex)(1.0f,0.0f) + cmul(dt,cmul(iihalfsgn,*U2));
	*U3 = cmul(dt,cmul(iihalfsgn,*U3));
}

__kernel void DepositCharge(__global tw_Complex *psi,__global tw_Float *J4)
{
	// PROTOCOL : Cell
	const tw_Int n = get_global_id(0);
	J4[n] += 0.5f*(psi[n].s0*psi[n].s0 + psi[n].s1*psi[n].s1);
}

__kernel void DepositCurrent(
	__global tw_Complex *psi0,
	__global tw_Complex *psi1,
	__global tw_Float *A4,
	__global tw_Float *J4,
	__global tw_Float *met_g,
	__global tw_Int *strip_g)
{
	// PROTOCOL : Strip
	tw_Int i,n;
	tw_Complex f11,f21,f12,f22,temp;
	tw_Float rc,ds;
	tw_Metrics met;
	tw_Strip strip;
	tw_cell_id cell,cell_inc;

	GetMetrics((tw_Float*)&met,met_g);
	GetStrip((tw_Int*)&strip,strip_g);
	const tw_Int ig = get_global_id(0);
	const tw_Int jg = get_global_id(1);
	const tw_Int kg = get_global_id(2);
	const tw_Int Ni = get_global_size(0);
	const tw_Int Nj = get_global_size(1);
	const tw_Int ys = (1-strip.di)*Ni + strip.di*(strip.dim+2);
	const tw_Int zs = ys*((1-strip.dj)*Nj + strip.dj*(strip.dim+2));
	n = ig + jg*ys + kg*zs;
	cell = (tw_cell_id)(0,ig,jg,kg);
	cell_inc = (tw_cell_id)(0,strip.di,strip.dj,strip.dk);

	const tw_Complex ii = (tw_Complex)(0.0f,1.0f);

	for (i=1;i<=strip.dim;i++)
	{
		n += strip.stride;
		cell += cell_inc;

		rc = met.car + met.cyl*(met.xo + ((tw_Float)cell.s1-0.5f)*met.dx);
		ds = strip.di*met.dx + strip.dj*met.dy*rc + strip.dk*met.dz;

		f11 = -0.5f*cmul(ii,psi0[n] - psi0[n-strip.stride])/ds + 0.25f*(A4[(n-strip.stride)*4+strip.axis]+A4[n*4+strip.axis])*psi0[n-strip.stride];
		f21 = -0.5f*cmul(ii,psi1[n] - psi1[n-strip.stride])/ds + 0.25f*(A4[(n-strip.stride)*4+strip.axis]+A4[n*4+strip.axis])*psi1[n-strip.stride];
		f12 = -0.5f*cmul(ii,psi0[n+strip.stride] - psi0[n])/ds + 0.25f*(A4[(n+strip.stride)*4+strip.axis]+A4[n*4+strip.axis])*psi0[n+strip.stride];
		f22 = -0.5f*cmul(ii,psi1[n+strip.stride] - psi1[n])/ds + 0.25f*(A4[(n+strip.stride)*4+strip.axis]+A4[n*4+strip.axis])*psi1[n+strip.stride];

		temp = cmul(conj(psi0[n]),f11+f21) + cmul(conj(psi1[n]),f11+f21);
		J4[n*4+strip.axis] += 0.25f*temp.s0;
		temp = cmul(conj(psi0[n]),f12+f22) + cmul(conj(psi1[n]),f12+f22);
		J4[n*4+strip.axis] += 0.25f*temp.s0;
	}
}

__kernel void UpdateSpin(
	__global tw_Complex *psi,
	__global tw_Complex *chi,
	__global tw_Float *A4,
	__global tw_Float *met_g)
{
	// PROTOCOL : Local
	tw_Metrics met;
	GetMetrics((tw_Float*)&met,met_g);

	const tw_Complex c16 = (tw_Complex)(16.0f,0.0f);

	const tw_Int i = get_global_id(0);
	const tw_Int j = get_global_id(1);
	const tw_Int k = get_global_id(2);
	const tw_Int xs = get_global_offset(0);
	const tw_Int ys = get_global_offset(1)*(get_global_size(0) + 2*get_global_offset(0));
	const tw_Int zs = get_global_offset(2)*(get_global_size(0) + 2*get_global_offset(0))*(get_global_size(1) + 2*get_global_offset(1));
	const tw_Int n = i*xs + j*ys + k*zs;
	const tw_Float adt = 0.00729735257f*met.dt;

	tw_Complex temp,U11,U12,U21,U22;
	tw_Float Bx,By,Bz,B2,denom;
	Bx = 0.5f*(A4[(n+ys)*4+3]-A4[(n-ys)*4+3])/met.dy;
	Bx -= 0.5f*(A4[(n+zs)*4+2]-A4[(n-zs)*4+2])/met.dz;
	By = (A4[(n+zs)*4+1]-A4[(n-zs)*4+1])/met.dz;
	By -= (A4[(n+xs)*4+3]-A4[(n-xs)*4+3])/met.dx;
	Bz = (A4[(n+xs)*4+2]-A4[(n-xs)*4+2])/met.dx;
	Bz -= (A4[(n+ys)*4+1]-A4[(n-ys)*4+1])/met.dy;

	B2 = Bx*Bx + By*By + Bz*Bz;
	denom = 16.0f+B2*adt*adt;
	U11 = (c16+adt*(tw_Complex)(-B2*adt,8.0f*Bz))/denom;
	U12 = 8.0f*(tw_Complex)(By,Bx)*adt/denom;
	U21 = 8.0f*(tw_Complex)(-By,Bx)*adt/denom;
	U22 = (c16+adt*(tw_Complex)(-B2*adt,-8.0f*Bz))/denom;

	temp = psi[n];
	psi[n] = cmul(U11,psi[n]) + cmul(U12,chi[n]);
	chi[n] = cmul(U21,temp) + cmul(U22,chi[n]);
}

__kernel void ApplyNumerator(
	__global tw_Complex *psi_in_out,
	__global tw_Float *A4,
	__global tw_Float *met_g,
	__global tw_Int *strip_g,
	const tw_Float partitionFactor,
	const tw_Float relax)
{
	// PROTOCOL : Strip
	tw_Int i,n;
	tw_Complex psi_prev,psi_now,T1,T2,T3;
	tw_Metrics met;
	tw_Strip strip;
	tw_cell_id cell,cell_inc;

	GetMetrics((tw_Float*)&met,met_g);
	GetStrip((tw_Int*)&strip,strip_g);
	const tw_Int ig = get_global_id(0);
	const tw_Int jg = get_global_id(1);
	const tw_Int kg = get_global_id(2);
	const tw_Int Ni = get_global_size(0);
	const tw_Int Nj = get_global_size(1);
	const tw_Int ys = (1-strip.di)*Ni + strip.di*(strip.dim+2);
	const tw_Int zs = ys*((1-strip.dj)*Nj + strip.dj*(strip.dim+2));
	n = ig + jg*ys + kg*zs;
	cell = (tw_cell_id)(0,ig,jg,kg);
	cell_inc = (tw_cell_id)(0,strip.di,strip.dj,strip.dk);

	const tw_Complex dt = (tw_Complex)(met.dt*(1.0f-relax),-met.dt*relax);

	// Zeroing ghost cells preserves Hermiticity
	// However, needs conditional in case of multi-node jobs
	//[n] = (tw_Complex)(0.0f,0.0f);
	//psi_in_out[n+strip.stride*(strip.dim+1)] = (tw_Complex)(0.0f,0.0f);

	psi_prev = psi_in_out[n];
	for (i=1;i<=strip.dim;i++)
	{
		n += strip.stride;
		cell += cell_inc;

		GetUnitaryFactor(n,A4,partitionFactor,cell,&met,&strip,&T1,&T2,&T3,dt,-1.0f);

		psi_now = psi_in_out[n];
		psi_in_out[n] = cmul(T1,psi_prev);
		psi_in_out[n] += cmul(T2,psi_now);
		psi_in_out[n] += cmul(T3,psi_in_out[n+strip.stride]);
		psi_prev = psi_now;
	}
}

__kernel void ApplyDenominator(
	__global tw_Complex *psi_in_out,
	__global tw_Float *A4,
	__global tw_Float *met_g,
	__global tw_Int *strip_g,
	const tw_Float partitionFactor,
	const tw_Float relax,
	__global tw_Complex *gam,
	__global tw_Complex *v,
	__global tw_Complex *w,
	__global tw_Complex *mpi_packet)
{
	// PROTOCOL : Strip
	tw_Int i,n,n0;
	tw_Complex bet,T1,T2,T3;
	tw_Metrics met;
	tw_Strip strip;
	tw_cell_id cell0,cell_inc,cell;

	// Each work item has to update strip.dim cells
	// In the strip direction, global_id = global_offset = 0, global_size = 1
	// Derivatives evaluated only in strip direction, so zero strides unnecessary

	GetMetrics((tw_Float*)&met,met_g);
	GetStrip((tw_Int*)&strip,strip_g);
	const tw_Int ig = get_global_id(0);
	const tw_Int jg = get_global_id(1);
	const tw_Int kg = get_global_id(2);
	const tw_Int Ni = get_global_size(0);
	const tw_Int Nj = get_global_size(1);
	const tw_Int Nk = get_global_size(2);
	const tw_Int ys = (1-strip.di)*Ni + strip.di*(strip.dim+2);
	const tw_Int zs = ys*((1-strip.dj)*Nj + strip.dj*(strip.dim+2));
	n0 = ig + jg*ys + kg*zs;
	cell0 = (tw_cell_id)(0,ig,jg,kg);
	cell_inc = (tw_cell_id)(0,strip.di,strip.dj,strip.dk);

	const tw_Int s = strip.stride;
	const tw_Complex dt = (tw_Complex)(met.dt*(1.0f-relax),-met.dt*relax);

	// UNCOUPLED NODE SOLUTION and GREEN'S FUNCTIONS (v,w)
	n = n0 + s;
	cell = cell0 + cell_inc;
	GetUnitaryFactor(n,A4,partitionFactor,cell,&met,&strip,&T1,&T2,&T3,dt,1.0f);
	bet = T2;
	psi_in_out[n] = cdiv(psi_in_out[n],T2);
	v[n] = cdiv(T1,T2); // T1 = theta
	w[n] = (tw_Complex)(0.0f,0.0f);
	for (i=2;i<=strip.dim;i++)
	{
		n = n0 + i*s;
		cell = cell0 + (int)i*cell_inc;
		gam[n] = cdiv(T3,bet);
		GetUnitaryFactor(n,A4,partitionFactor,cell,&met,&strip,&T1,&T2,&T3,dt,1.0f);
		bet = T2 - cmul(T1,gam[n]);
		psi_in_out[n] = cdiv(psi_in_out[n] - cmul(T1,psi_in_out[n-s]),bet);
		v[n] = cdiv(-cmul(T1,v[n-s]),bet);
		w[n] = (tw_Complex)(0.0f,0.0f);
	}
	w[n] = cdiv(T3 - cmul(T1,w[n-s]),bet); // T3 = eta
	for (i=strip.dim-1;i>=1;i--)
	{
		n = n0 + i*s;
		psi_in_out[n] -= cmul(gam[n+s],psi_in_out[n+s]);
		v[n] -= cmul(gam[n+s],v[n+s]);
		w[n] -= cmul(gam[n+s],w[n+s]);
	}

	// FORM MPI PACKET
	n = 8*strip.di*(jg + kg*Nj);
	n += 8*strip.dj*(kg + ig*Nk);
	n += 8*strip.dk*(ig + jg*Ni);
	mpi_packet[n+0] = psi_in_out[n0+s];
	mpi_packet[n+1] = v[n0+s*strip.dim];
	mpi_packet[n+2] = w[n0+s*strip.dim];
	mpi_packet[n+3] = v[n0+s];
	mpi_packet[n+4] = w[n0+s];
	mpi_packet[n+5] = psi_in_out[n0+s*strip.dim];
	mpi_packet[n+6] = (tw_Complex)(0.0f,0.0f);
	mpi_packet[n+7] = (tw_Complex)(0.0f,0.0f);
}

__kernel void Stitch(
	__global tw_Complex *psi_in_out,
	__global tw_Int *strip_g,
	__global tw_Complex *v,
	__global tw_Complex *w,
	__global tw_Complex *mpi_packet)
{
	// PROTOCOL : Local
	tw_Int n,m;
	tw_Strip strip;

	GetStrip((tw_Int*)&strip,strip_g);

	const tw_Int ig = get_global_id(0);
	const tw_Int jg = get_global_id(1);
	const tw_Int kg = get_global_id(2);
	const tw_Int Ni = get_global_size(0) + 2*get_global_offset(0);
	const tw_Int Nj = get_global_size(1) + 2*get_global_offset(1);
	const tw_Int Nk = get_global_size(2) + 2*get_global_offset(2);

	n = ig + jg*Ni + kg*Ni*Nj;
	m = 8*strip.di*(jg + kg*Nj);
	m += 8*strip.dj*(kg + ig*Nk);
	m += 8*strip.dk*(ig + jg*Ni);

	psi_in_out[n] -= cmul(mpi_packet[m+6],v[n]) + cmul(mpi_packet[m+7],w[n]);
}

__kernel void KGAdvance1(
	__global tw_Float *psi_r,
	__global tw_Float *psi_i,
	__global tw_Float *A4,
	__global tw_Float *met_g,
	const tw_Float m0,
	const tw_Float q0)
{
	LOCAL_PROTOCOL_MACRO
	const tw_Float dt = met.dt;
	const tw_Float dth = (tw_Float)0.5f*met.dt;
	// Update psi
	// unitary operator of time translation for diagonal part of Hamiltonian
	const tw_Float Ur = cos(dt*q0*A4[n]);
	const tw_Float Ui = -sin(dt*q0*A4[n]);
	// midpoint method for off-diagonal term
	tw_Float Dr,Di,ans_r,ans_i;
	Dr = m0*psi_i[n+cs];
	Di = -m0*psi_r[n+cs];
	ans_r = psi_r[n];
	ans_i = psi_i[n];
	split_step_advance(&ans_r,&ans_i,Ur,Ui,Dr*dth,Di*dth);
	psi_r[n] = ans_r;
	psi_i[n] = ans_i;
}

__kernel void KGAdvance2(
	__global tw_Float *psi_r,
	__global tw_Float *psi_i,
	__global tw_Float *A4,
	__global tw_Float *met_g,
	const tw_Float m0,
	const tw_Float q0)
{
	LOCAL_PROTOCOL_MACRO
	const tw_Float dt = met.dt;
	const tw_Float dth = (tw_Float)0.5f*met.dt;
	const tw_Float fx = (tw_Float)1.0f/met.dx;
	const tw_Float fy = (tw_Float)1.0f/met.dy;
	const tw_Float fz = (tw_Float)1.0f/met.dz;
	// Update chi
	// unitary operator of time translation for diagonal part of Hamiltonian
	const tw_Float Ur = cos(dt*q0*A4[n]);
	const tw_Float Ui = -sin(dt*q0*A4[n]);
	// Evaluate i(Dk)^2 assuming div(A)=0
	// This adds i*del^2(psi)
	tw_Float Dr,Di,ans_r,ans_i;
	Dr = -(diff2(psi_i,n,xs)*fx*fx + diff2(psi_i,n,ys)*fy*fy + diff2(psi_i,n,zs)*fz*fz);
	Di = (diff2(psi_r,n,xs)*fx*fx + diff2(psi_r,n,ys)*fy*fy + diff2(psi_r,n,zs)*fz*fz);
	// This adds 2q*div(A*psi)
	Dr += 2*q0*fx*(sfwd(psi_r,n,xs)*sfwd(A4,n+1*cs,xs) - sbak(psi_r,n,xs)*sbak(A4,n+1*cs,xs));
	Dr += 2*q0*fy*(sfwd(psi_r,n,ys)*sfwd(A4,n+2*cs,ys) - sbak(psi_r,n,ys)*sbak(A4,n+2*cs,ys));
	Dr += 2*q0*fz*(sfwd(psi_r,n,zs)*sfwd(A4,n+3*cs,zs) - sbak(psi_r,n,zs)*sbak(A4,n+3*cs,zs));
	Di += 2*q0*fx*(sfwd(psi_i,n,xs)*sfwd(A4,n+1*cs,xs) - sbak(psi_i,n,xs)*sbak(A4,n+1*cs,xs));
	Di += 2*q0*fy*(sfwd(psi_i,n,ys)*sfwd(A4,n+2*cs,ys) - sbak(psi_i,n,ys)*sbak(A4,n+2*cs,ys));
	Di += 2*q0*fz*(sfwd(psi_i,n,zs)*sfwd(A4,n+3*cs,zs) - sbak(psi_i,n,zs)*sbak(A4,n+3*cs,zs));
	// This adds -i*q^2*A^2*psi
	Dr += q0*q0*(sqr(A4[n+1*cs])+sqr(A4[n+2*cs])+sqr(A4[n+3*cs]))*psi_i[n];
	Di -= q0*q0*(sqr(A4[n+1*cs])+sqr(A4[n+2*cs])+sqr(A4[n+3*cs]))*psi_r[n];
	// Finish by incorporating mass
	Dr = Dr/m0 + m0*psi_i[n];
	Di = Di/m0 - m0*psi_r[n];
	// Perform the integration
	ans_r = psi_r[n+cs];
	ans_i = psi_i[n+cs];
	split_step_advance(&ans_r,&ans_i,Ur,Ui,Dr*dth,Di*dth);
	psi_r[n+cs] = ans_r;
	psi_i[n+cs] = ans_i;
}

__kernel void DiracAdvance(
	__global tw_Float *psi_r,
	__global tw_Float *psi_i,
	__global tw_Float *A4,
	__global tw_Float *met_g,
	const tw_Float m0,
	const tw_Float q0,
	const int OUT1,
	const int OUT2,
	const int IN1,
	const int IN2)
{
	// PROTOCOL : Local
	// Advances 2 components of the Dirac bispinor
	// Intended to be called twice in leap-frog fashion
	// Arguments determine which components are leap-frogged
	// Mass must change sign from one call to the next
	LOCAL_PROTOCOL_MACRO
	tw_Float Dr,Di,Ar,Ai;

	const tw_Int out1 = n + cs*OUT1;
	const tw_Int out2 = n + cs*OUT2;
	const tw_Int in1 = n + cs*IN1;
	const tw_Int in2 = n + cs*IN2;
	const tw_Float dt = met.dt;
	const tw_Float dth = 0.5*met.dt;
	const tw_Float fx = 0.5/met.dx;
	const tw_Float fy = 0.5/met.dy;
	const tw_Float fz = 0.5/met.dz;

	const tw_Float phi = A4[n];
	const tw_Float Ax = A4[n + 1*cs];
	const tw_Float Ay = A4[n + 2*cs];
	const tw_Float Az = A4[n + 3*cs];

	const tw_Float Ur = cos(dt*(m0 + q0*phi));
	const tw_Float Ui = -sin(dt*(m0 + q0*phi));

	Dr = -fz*diff(psi_r,in1,zs) - fx*diff(psi_r,in2,xs) - fy*diff(psi_i,in2,ys);
	Di = -fz*diff(psi_i,in1,zs) - fx*diff(psi_i,in2,xs) + fy*diff(psi_r,in2,ys);
	Dr += q0*(-Az*psi_i[in1] - Ax*psi_i[in2] + Ay*psi_r[in2]);
	Di += q0*(Az*psi_r[in1] + Ax*psi_r[in2] + Ay*psi_i[in2]);
	Ar = psi_r[out1];
	Ai = psi_i[out1];
	split_step_advance(&Ar,&Ai,Ur,Ui,Dr*dth,Di*dth);
	psi_r[out1] = Ar;
	psi_i[out1] = Ai;

	Dr = fz*diff(psi_r,in2,zs) - fx*diff(psi_r,in1,xs) + fy*diff(psi_i,in1,ys);
	Di = fz*diff(psi_i,in2,zs) - fx*diff(psi_i,in1,xs) - fy*diff(psi_r,in1,ys);
	Dr += q0*(Az*psi_i[in2] - Ax*psi_i[in1] - Ay*psi_r[in1]);
	Di += q0*(-Az*psi_r[in2] + Ax*psi_r[in1] - Ay*psi_i[in1]);
	Ar = psi_r[out2];
	Ai = psi_i[out2];
	split_step_advance(&Ar,&Ai,Ur,Ui,Dr*dth,Di*dth);
	psi_r[out2] = Ar;
	psi_i[out2] = Ai;
}

// __kernel void UpdatePotentials(
// 	__global tw_Float *A4,
// 	__global tw_Float *met_g,
// 	__global tw_Float *waves,
// 	const int num_waves,
// 	const tw_Float B0,
// 	const tw_Float t)
// {
// 	POINT_PROTOCOL_METRIC_MACRO
// 	tw_Float A[3],z0;
// 	tw_Float a0,w0,u[3],w[3],f[3],tp[4];
// 	const tw_Float x = met.xo + met.dx * ((tw_Float)ig-0.5);
// 	const tw_Float y = met.yo + met.dy * ((tw_Float)jg-0.5);
// 	const tw_Float z = met.zo + met.dz * ((tw_Float)kg-0.5);
//
// 	A[0] = -0.5*B0*y;
// 	A[1] = 0.5*B0*x;
// 	A[2] = 0.0;
//
// 	for (tw_Int wv=0;wv<num_waves;wv++)
// 	{
// 		GetWave(waves,wv,&a0,&w0,u,w,f,tp);
// 		z0 = 0.0;
// 		z0 += w[0]*(x - f[0]);
// 		z0 += w[1]*(y - f[1]);
// 		z0 += w[2]*(z - f[2]);
// 		a0 *= cos(w0*z0 - w0*(t-tp[0]))*PulseShapeFactor(t-z0,tp);
// 		A[0] += u[0]*a0;
// 		A[1] += u[1]*a0;
// 		A[2] += u[2]*a0;
// 	}
// 	A4[n+1*cs] = A[0];
// 	A4[n+2*cs] = A[1];
// 	A4[n+3*cs] = A[2];
// }
