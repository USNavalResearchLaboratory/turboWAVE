module;

#include "tw_includes.h"

export module parabolic;
import input;
import compute_tool;
import region;
import fields;
import numerics;

using namespace tw::bc;

export enum tw_polarization_type {linearPolarization,circularPolarization,radialPolarization};

export struct LaserPropagator:ComputeTool
{
	tw::Float w0,dt;
	tw_polarization_type polarization;
	bool movingWindow;

	LaserPropagator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~LaserPropagator();
	virtual void SetData(tw::Float w0,tw::Float dt,tw_polarization_type pol,bool mov,MetricSpace *refined = NULL);
	void SetBoundaryConditions(ComplexField& a0,ComplexField& a1,ComplexField& chi);
	virtual void Advance(ComplexField& a0,ComplexField& a1,ComplexField& chi) = 0;
};

export struct EigenmodePropagator:LaserPropagator
{
	// following selects number of radial modes to keep and absorbing layers
	tw::Int modes,layers;
	// causality = 0.0 gives centered differencing
	// causality = -1.0 uses forward (downwind) differencing
	// causality = 1.0 uses backward (upwind) differencing
	// non-integral values are legal, and give a linear combination
	tw::Float causality,dampingTime;
	std::valarray<tw::Float> eigenvalue,hankel,inverseHankel;
	GlobalIntegrator<tw::Complex>* globalIntegrator;

	EigenmodePropagator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~EigenmodePropagator();
	virtual void SetData(tw::Float w0,tw::Float dt,tw_polarization_type pol,bool mov,MetricSpace *refined = NULL);
	virtual void Advance(ComplexField& a0,ComplexField& a1,ComplexField& chi);
};

export struct ADIPropagator:LaserPropagator
{
	GlobalIntegrator<tw::Complex>* xGlobalIntegrator;
	GlobalIntegrator<tw::Complex>* yGlobalIntegrator;
	std::valarray<tw::Complex> X1,X2,X3;
	std::valarray<tw::Complex> Y1,Y2,Y3;
	bool evenTime;

	ADIPropagator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~ADIPropagator();
	virtual void SetData(tw::Float w0,tw::Float dt,tw_polarization_type pol,bool mov,MetricSpace *refined = NULL);
	virtual void Advance(ComplexField& a0,ComplexField& a1,ComplexField& chi);
};

struct ConservativePropagator:LaserPropagator
{
};

export struct SchroedingerPropagator:ComputeTool
{
	std::vector<GlobalIntegrator<tw::Complex>*> globalIntegrator;

	SchroedingerPropagator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~SchroedingerPropagator();
	virtual void DepositCurrent(const tw::grid::axis& axis,ComplexField& psi0,ComplexField& psi1,Field& A4,Field& J4,tw::Complex dt);
	virtual void ApplyNumerator(const tw::grid::axis& axis,ComplexField& psi,Field& A4,bool keepA2Term,tw::Complex dt);
	virtual void ApplyDenominator(const tw::grid::axis& axis,ComplexField& psi,Field& A4,bool keepA2Term,tw::Complex dt);
	virtual void UpdateSpin(ComplexField& psi,ComplexField& chi,Field& A4,tw::Float adt);
};

export struct ParabolicSolver:BoundedTool
{
	// inhomogeneous boundary conditions are handled as described in elliptic.h
	std::vector<GlobalIntegrator<tw::Float>*> globalIntegrator;

	ParabolicSolver(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~ParabolicSolver();

	void FixTemperature(Field& T,const Element& e,Region* theRegion,const tw::Float& T0);
	void FormOperatorStencil(tw::Float *D1,tw::Float *D2,const ScalarField& fluxMask,Field *coeff,tw::Int c,const tw::strip& s,tw::Int i);
	virtual void Advance(const tw::grid::axis& axis,ScalarField& psi,ScalarField& fluxMask,tw::Float coeff,tw::Float dt);
	virtual void Advance(ScalarField& psi,ScalarField& fluxMask,tw::Float coeff,tw::Float dt);
	virtual void Advance(	const tw::grid::axis& axis,
							Field& psi,
							tw::Int psi_idx,
							ScalarField& fluxMask,
							Field *coeff1,
							tw::Int c1_idx,
							Field *coeff2,
							tw::Int c2_idx,
							tw::Float dt);
	virtual void Advance(	Field& psi,
							tw::Int psi_idx,
							ScalarField& fluxMask,
							Field *coeff1,
							tw::Int c1_idx,
							Field *coeff2,
							tw::Int c2_idx,
							tw::Float dt);
};

export struct IsotropicPropagator:ComputeTool
{
	GlobalIntegrator<tw::Complex>* zGlobalIntegrator;
	std::valarray<tw::Complex> Z1,Z2,Z3;

	IsotropicPropagator(const std::string& name,MetricSpace *m,Task *tsk);
	virtual ~IsotropicPropagator();

	void SetupIncomingWaveLeft(const tw::strip& s,ComplexField& amplitude,tw::Complex a0,tw::Complex a1,tw::Complex w0) const
	{
		amplitude(s,0) = a1 - a0 + ii*w0*space->dl(s,1,3)*half*(a0+a1);
	}

	void SetupIncomingWaveRight(const tw::strip& s,ComplexField& amplitude,tw::Complex aN,tw::Complex aN1,tw::Complex w0) const
	{
		amplitude(s,space->Dim(s.Axis())+1) = aN1 - aN - ii*w0*space->dl(s,space->Dim(s.Axis())+1,3)*half*(aN+aN1);
	}

	virtual void Advance(ComplexField& amplitude, // this is the array to update
		ComplexField& refractiveIndex, // index of refraction array
		ScalarField& nu_e, // electron collision frequency array
		tw::Float laserFrequency,tw::Float dt);
};

/////////////////////////////
//                         //
//     Laser Propagator    //
//                         //
/////////////////////////////


LaserPropagator::LaserPropagator(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	w0 = 10.0;
	polarization = linearPolarization;
	movingWindow = true;
}

LaserPropagator::~LaserPropagator()
{
}

void LaserPropagator::SetData(tw::Float w0,tw::Float dt,tw_polarization_type pol,bool mov,MetricSpace *refined)
{
	this->w0 = w0;
	this->dt = dt;
	polarization = pol;
	movingWindow = mov;
	if (refined != NULL) {
		space = refined;
	}
}

void LaserPropagator::SetBoundaryConditions(ComplexField& a0,ComplexField& a1,ComplexField& chi)
{
	tw::bc::fld xl,xh,yl,yh;

	xl = space->cyl==1.0 && polarization==radialPolarization ? fld::dirichletWall : fld::neumannWall;
	xh = fld::neumannWall;
	yl = fld::neumannWall;
	yh = fld::neumannWall;

	a0.SetBoundaryConditions(tw::grid::x,xl,xh);
	a1.SetBoundaryConditions(tw::grid::x,xl,xh);
	a0.SetBoundaryConditions(tw::grid::y,yl,yh);
	a1.SetBoundaryConditions(tw::grid::y,yl,yh);

	xl = space->cyl==1.0 ? fld::neumannWall : fld::dirichletWall;

	chi.SetBoundaryConditions(tw::grid::x,xl,fld::dirichletCell);
	chi.SetBoundaryConditions(tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	chi.SetBoundaryConditions(tw::grid::z,fld::dirichletCell,fld::dirichletCell);
}



/////////////////////////////
//                         //
//  Eigenmode Propagator   //
//                         //
/////////////////////////////


EigenmodePropagator::EigenmodePropagator(const std::string& name,MetricSpace *m,Task *tsk) : LaserPropagator(name,m,tsk)
{
	if (space->car==1.0 && !space->TransversePowersOfTwo())
		throw tw::FatalError("eigenmode propagator in Cartesian coordinates requires transverse dimensions be powers of two.");

	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	modes = space->GlobalDim(1);
	causality = 0.0;
	dampingTime = 1e10;
	layers = 8;
	globalIntegrator = new GlobalIntegrator<tw::Complex>(&task->strip[3],xDim*yDim,zDim);

	directives.Add("modes",new tw::input::Int(&modes),false);
	directives.Add("damping time",new tw::input::Float(&dampingTime),false);
	directives.Add("absorbing layers",new tw::input::Int(&layers),false);
	directives.Add("causality",new tw::input::Float(&causality),false);
}

EigenmodePropagator::~EigenmodePropagator()
{
	delete globalIntegrator;
}

void EigenmodePropagator::SetData(tw::Float w0,tw::Float dt,tw_polarization_type pol,bool mov,MetricSpace *refined)
{
	LaserPropagator::SetData(w0,dt,pol,mov,refined);
	if (refined!=NULL) {
		delete globalIntegrator;
		globalIntegrator = new GlobalIntegrator<tw::Complex>(&task->strip[3],refined->Dim(1)*refined->Dim(2),refined->Dim(3));
	}
	if (space->car!=1.0)
		ComputeTransformMatrices(fld::dirichletWall,eigenvalue,hankel,inverseHankel,space,task);
	if (modes > space->GlobalDim(1))
		throw tw::FatalError("more modes than radial cells");
	if (space->cyl==1.0 && pol==radialPolarization)
		throw tw::FatalError("eigenmode propagator incompatible with radial polarization.");
}

void EigenmodePropagator::Advance(ComplexField& a0,ComplexField& a1,ComplexField& chi)
{
	// unlike the papers, we take a~ = (a/2)exp(i(kx-wt)) + cc
	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);
	const tw::Float dz = space->dx(3);

	std::valarray<tw::Float> localEig(xDim+2),chi_ref(zDim+2);

	const tw::Complex T1 = -(1+causality)/(2*dz*dt);
	const tw::Complex T2 = ii*w0/dt + causality/(dz*dt);
	const tw::Complex T3 = (1-causality)/(2*dz*dt);

	if (space->car!=1.0)
	{
		for (tw::Int i=1;i<=xDim;i++)
			localEig[i] = eigenvalue[space->GlobalCellIndex(i,3)-1];
	}

	// Setup reference chi and explicit current
	// Current goes into chi, destroying chi
	if (space->car==1.0)
	{
		tw::Int i = space->LocalCellIndex(space->GlobalDim(1)/2,1) + 1;
		tw::Int j = space->LocalCellIndex(space->GlobalDim(2)/2,2) + 1;
		if (i>=1 && i<=xDim && j>=1 && j<=yDim)
			for (tw::Int k=0;k<=zDim+1;k++)
				chi_ref[k] = std::real(chi(i,j,k));
		else
			chi_ref = 0.0;
		task->strip[1].AllSum(&chi_ref[0],&chi_ref[0],sizeof(tw::Float)*chi_ref.size(),0);
		task->strip[2].AllSum(&chi_ref[0],&chi_ref[0],sizeof(tw::Float)*chi_ref.size(),0);
	}
	else
	{
		for (tw::Int k=0;k<=zDim+1;k++)
			chi_ref[k] = std::real(chi(1,1,k));
		task->strip[1].Bcast(&chi_ref[0],sizeof(tw::Float)*chi_ref.size(),0);
	}
	for (auto s : StripRange(*space,3,strongbool::yes))
		for (tw::Int k=0;k<=zDim+1;k++)
			chi(s,k) = (chi(s,k) - chi_ref[k])*a1(s,k);

	// Diagonalize Laplacian and solve in laser frame

	if (space->car==1.0)
	{
		a0.TransverseFFT();
		chi.TransverseFFT();
	}
	else
	{
		a0.Hankel(modes,hankel);
		chi.Hankel(modes,hankel);
	}

	// NOTE: VectorStripRange does not work for complex fields
	#pragma omp parallel
	{
		tw::Float lambda;
		std::valarray<tw::Complex> T2z(zDim),source(zDim),ans(zDim);
		StripRange range(*space,3,strongbool::no);
		for (auto it=range.begin();it!=range.end();++it)
		{
			tw::strip s = *it;
			if (space->car==1.0)
				lambda = a0.CyclicEigenvalue(s.dcd1(0),s.dcd2(0));
			else
				lambda = localEig[s.dcd1(0)];
			for (tw::Int k=1;k<=zDim;k++)
			{
				T2z[k-1] = T2 - tw::Float(0.5)*(lambda + chi_ref[k]);
				source[k-1] = T1*a0(s,k-1) + T2z[k-1]*a0(s,k) + T3*a0(s,k+1) - chi(s,k);
				T2z[k-1] = T2 + tw::Float(0.5)*(lambda + chi_ref[k]);
			}

			TriDiagonal(ans,source,T1,T2z,T3);

			for (tw::Int k=1;k<=zDim;k++)
			{
				a0(s,k) = a1(s,k);
				a1(s,k) = ans[k-1];
			}

			globalIntegrator->SetMatrix(it.global_count(),T1,T2z,T3,T1,T3);
			globalIntegrator->SetData(it.global_count(),&a1(s,0),a1.Stride(3)/2); // need complex stride and stride[0]=1
		}
	}
	CopyGhostCellData(a0,All(a0),a1,All(a1));

	// Take into account the other nodes

	globalIntegrator->Parallelize();

	// Back to real space

	if (space->car==1.0)
	{
		a1.InverseTransverseFFT();
		chi.InverseTransverseFFT();
	}
	else
	{
		a1.InverseHankel(modes,inverseHankel);
		chi.InverseHankel(modes,inverseHankel);
	}

	// Damping layer

	tw::Int xstart = space->GlobalDim(1) - layers + 1;
	tw::Int xstart_local = space->LocalCellIndex(xstart,1);
	tw::Int first = xstart_local;
	if (first<space->LFG(1))
		first = space->LFG(1);
	for (tw::Int i=first;i<=space->UFG(1);i++)
		for (tw::Int k=0;k<=zDim+1;k++)
			a1(i,1,k) *= exp(-(i-xstart_local+1)*dt/dampingTime/tw::Float(layers));

	a1.ApplyBoundaryCondition();

	if (task->n0[3]==MPI_PROC_NULL && zDim>1)
		for (auto s : StripRange(*space,3,strongbool::yes))
			a1(s,space->LFG(3)) = a0(s,space->LFG(3));

	if (task->n1[3]==MPI_PROC_NULL && zDim>1)
		for (auto s : StripRange(*space,3,strongbool::yes))
			a1(s,space->UFG(3)) = a0(s,space->UFG(3));
}


////////////////////////
//                    //
//   ADI Propagator   //
//                    //
////////////////////////


ADIPropagator::ADIPropagator(const std::string& name,MetricSpace *m,Task *tsk) : LaserPropagator(name,m,tsk)
{
	evenTime = true;
	xGlobalIntegrator = NULL;
	yGlobalIntegrator = NULL;

	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	X1.resize(xDim);
	X2.resize(xDim);
	X3.resize(xDim);
	Y1.resize(yDim);
	Y2.resize(yDim);
	Y3.resize(yDim);

	if (xDim>1)
		xGlobalIntegrator = new GlobalIntegrator<tw::Complex>(&task->strip[1],yDim*zDim,xDim);

	if (yDim>1)
		yGlobalIntegrator = new GlobalIntegrator<tw::Complex>(&task->strip[2],xDim*zDim,yDim);
}

ADIPropagator::~ADIPropagator()
{
	if (xGlobalIntegrator!=NULL)
		delete xGlobalIntegrator;
	if (yGlobalIntegrator!=NULL)
		delete yGlobalIntegrator;
}

void ADIPropagator::SetData(tw::Float w0,tw::Float dt,tw_polarization_type pol,bool mov,MetricSpace *refined)
{
	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);
	LaserPropagator::SetData(w0,dt,pol,mov,refined);
	if (refined != NULL && xGlobalIntegrator != NULL) {
		delete xGlobalIntegrator;
		xGlobalIntegrator = new GlobalIntegrator<tw::Complex>(&task->strip[1],yDim*zDim,xDim);
	}
	if (refined != NULL && yGlobalIntegrator != NULL) {
		delete yGlobalIntegrator;
		yGlobalIntegrator = new GlobalIntegrator<tw::Complex>(&task->strip[2],xDim*zDim,yDim);
	}
}

void ADIPropagator::Advance(ComplexField& a0,ComplexField& a1,ComplexField& chi)
{
	// unlike the papers, we take a~ = (a/2)exp(i(kx-wt)) + cc
	// this solves (delperp^2 + 2iw0Dt + 2Dzt)a = -chi*a
	// where use the estimate 2Dzt(a) ~ (i/w0)*Dz[(chi + delperp^2)a]
	// here, we must remember that aNow is leapfrogged, so effective timestep is 2*dt
	// thus, we end up with (1 + i*dt*H0)a1 = (1 - i*dt*H0)a0 - 2*i*dt*H1*aNow
	// with H0 = -Dx^2/(2*w0) and H1 = -(1/w0)*(chi/2 + Dzt + Dy^2/2)

	tw::Int i,j,k;
	tw::Float Ax0,Ax1,Ay0,Ay1,Vol,radialPolarizationFactor;
	tw::Complex D1,D2,D3,D4,Dzt;

	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	radialPolarizationFactor = space->cyl==1.0 && polarization==radialPolarization ? 1.0 : 0.0;

	ComplexField aNow;
	aNow.Initialize(*space,task);
	aNow = a1;

	if (xDim==1 && yDim==1)
	{
		for (k=1;k<=zDim;k++)
		{
			a1(1,1,k) = a0(1,1,k) + (ii*dt/w0)*chi(1,1,k)*aNow(1,1,k);
			a1(1,1,k) -= (dt/sqr(w0))*(chi(1,1,k+1)*aNow(1,1,k+1) - chi(1,1,k-1)*aNow(1,1,k-1))/space->dL(1,1,k,3);
		}
	}

	if ((!evenTime && xDim>1) || (yDim==1 && xDim>1))
	{
		std::valarray<tw::Complex> src(xDim),ans(xDim);

		for (k=1;k<=zDim;k++)
			for (j=1;j<=yDim;j++)
			{
				for (i=1;i<=xDim;i++)
				{
					Vol = space->dS(i,j,k,0);
					Ax0 = space->dS(i,j,k,1);
					Ax1 = space->dS(i+1,j,k,1);
					Ay0 = space->dS(i,j,k,2);
					Ay1 = space->dS(i,j+1,k,2);
					D1 = (Ax0/Vol) / space->dl(i,j,k,1);
					D2 = (Ax1/Vol) / space->dl(i+1,j,k,1);
					D3 = (Ay0/Vol) / space->dl(i,j,k,2);
					D4 = (Ay1/Vol) / space->dl(i,j+1,k,2);
					Dzt = (ii/w0)/space->dL(i,j,k,3);
					src[i-1] = (ii*w0/dt)*a0(i,j,k) - chi(i,j,k)*aNow(i,j,k);
					src[i-1] -= half*(D1*a0(i-1,j,k) - (D1 + D2 + radialPolarizationFactor/sqr(space->X(i,1)))*a0(i,j,k) + D2*a0(i+1,j,k));
					src[i-1] -= D3*aNow(i,j-1,k) - (D3 + D4)*aNow(i,j,k) + D4*aNow(i,j+1,k);
					src[i-1] -= Dzt*(chi(i,j,k+1)*aNow(i,j,k+1) - chi(i,j,k-1)*aNow(i,j,k-1));
					src[i-1] -= Dzt*(D1*aNow(i-1,j,k+1) - (D1+D2)*aNow(i,j,k+1) + D2*aNow(i+1,j,k+1));
					src[i-1] += Dzt*(D1*aNow(i-1,j,k-1) - (D1+D2)*aNow(i,j,k-1) + D2*aNow(i+1,j,k-1));
					src[i-1] -= Dzt*(D3*aNow(i,j-1,k+1) - (D3+D4)*aNow(i,j,k+1) + D4*aNow(i,j+1,k+1));
					src[i-1] += Dzt*(D3*aNow(i,j-1,k-1) - (D3+D4)*aNow(i,j,k-1) + D4*aNow(i,j+1,k-1));
					X1[i-1] = half*D1;
					X2[i-1] = ii*w0/dt - half*(D1 + D2 + radialPolarizationFactor/sqr(space->X(i,1)));
					X3[i-1] = half*D2;
				}

				if (task->n0[1]==MPI_PROC_NULL)
					a1.AdjustTridiagonalForBoundaries(tw::grid::x,tw::grid::low,X1,X2,X3,src,tw::Complex(0.0));
				if (task->n1[1]==MPI_PROC_NULL)
					a1.AdjustTridiagonalForBoundaries(tw::grid::x,tw::grid::high,X1,X2,X3,src,tw::Complex(0.0));

				TriDiagonal(ans,src,X1,X2,X3);
				for (i=1;i<=xDim;i++)
					a1(i,j,k) = ans[i-1];

				xGlobalIntegrator->SetMatrix((k-1)*yDim + (j-1),X1,X2,X3);
				xGlobalIntegrator->SetData((k-1)*yDim + (j-1),&a1(0,j,k),a1.Stride(1)/2); // need complex stride and stride[0]=1
			}

		xGlobalIntegrator->Parallelize();
		a1.DownwardCopy(tw::grid::y,1);
		a1.UpwardCopy(tw::grid::y,1);
		a1.DownwardCopy(tw::grid::z,1);
		a1.UpwardCopy(tw::grid::z,1);
		a1.ApplyBoundaryCondition();
	}

	if ((evenTime && yDim>1) || (xDim==1 && yDim>1))
	{
		std::valarray<tw::Complex> src(yDim),ans(yDim);

		for (k=1;k<=zDim;k++)
			for (i=1;i<=xDim;i++)
			{
				for (j=1;j<=yDim;j++)
				{
					Vol = space->dS(i,j,k,0);
					Ax0 = space->dS(i,j,k,1);
					Ax1 = space->dS(i+1,j,k,1);
					Ay0 = space->dS(i,j,k,2);
					Ay1 = space->dS(i,j+1,k,2);
					D1 = (Ax0/Vol) / space->dl(i,j,k,1);
					D2 = (Ax1/Vol) / space->dl(i+1,j,k,1);
					D3 = (Ay0/Vol) / space->dl(i,j,k,2);
					D4 = (Ay1/Vol) / space->dl(i,j+1,k,2);
					Dzt = (ii/w0)/space->dL(i,j,k,3);
					src[j-1] = (ii*w0/dt)*a0(i,j,k) - chi(i,j,k)*aNow(i,j,k);
					src[j-1] -= half*(D3*a0(i,j-1,k) - (D3 + D4)*a0(i,j,k) + D4*a0(i,j+1,k));
					src[j-1] -= D1*aNow(i-1,j,k) - (D1 + D2)*aNow(i,j,k) + D2*aNow(i+1,j,k);
					src[j-1] -= Dzt*(chi(i,j,k+1)*aNow(i,j,k+1) - chi(i,j,k-1)*aNow(i,j,k-1));
					src[j-1] -= Dzt*(D1*aNow(i-1,j,k+1) - (D1+D2)*aNow(i,j,k+1) + D2*aNow(i+1,j,k+1));
					src[j-1] += Dzt*(D1*aNow(i-1,j,k-1) - (D1+D2)*aNow(i,j,k-1) + D2*aNow(i+1,j,k-1));
					src[j-1] -= Dzt*(D3*aNow(i,j-1,k+1) - (D3+D4)*aNow(i,j,k+1) + D4*aNow(i,j+1,k+1));
					src[j-1] += Dzt*(D3*aNow(i,j-1,k-1) - (D3+D4)*aNow(i,j,k-1) + D4*aNow(i,j+1,k-1));
					Y1[j-1] = half*D3;
					Y2[j-1] = ii*w0/dt - half*(D3 + D4);
					Y3[j-1] = half*D3;
				}

				if (task->n0[2]==MPI_PROC_NULL)
					a1.AdjustTridiagonalForBoundaries(tw::grid::y,tw::grid::low,Y1,Y2,Y3,src,tw::Complex(0.0));
				if (task->n1[2]==MPI_PROC_NULL)
					a1.AdjustTridiagonalForBoundaries(tw::grid::y,tw::grid::high,Y1,Y2,Y3,src,tw::Complex(0.0));

				TriDiagonal(ans,src,Y1,Y2,Y3);
				for (j=1;j<=yDim;j++)
					a1(i,j,k) = ans[j-1];

				yGlobalIntegrator->SetMatrix((k-1)*xDim + (i-1),Y1,Y2,Y3);
				yGlobalIntegrator->SetData((k-1)*xDim + (i-1),&a1(i,0,k),a1.Stride(2)/2); // need complex stride and stride[0]=1
			}

		yGlobalIntegrator->Parallelize();
		a1.DownwardCopy(tw::grid::x,1);
		a1.UpwardCopy(tw::grid::x,1);
		a1.DownwardCopy(tw::grid::z,1);
		a1.UpwardCopy(tw::grid::z,1);
		a1.ApplyBoundaryCondition();
	}

	a0 = aNow;
	if (task->n0[3]==MPI_PROC_NULL && zDim>1)
		for (j=space->LFG(2);j<=space->UFG(2);j++)
			for (i=space->LFG(1);i<=space->UFG(1);i++)
				a1(i,j,space->LFG(3)) = a0(i,j,space->LFG(3));

	if (task->n1[3]==MPI_PROC_NULL && zDim>1)
		for (j=space->LFG(2);j<=space->UFG(2);j++)
			for (i=space->LFG(1);i<=space->UFG(1);i++)
				a1(i,j,space->UFG(3)) = a0(i,j,space->UFG(3));

	evenTime = !evenTime;
}



//////////////////////////////////////////
//                                      //
// Propagator for Schroedinger Equation //
//                                      //
//////////////////////////////////////////


SchroedingerPropagator::SchroedingerPropagator(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	globalIntegrator.assign(4,NULL);

	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	if (xDim>1)
		globalIntegrator[1] = new GlobalIntegrator<tw::Complex>(&task->strip[1],yDim*zDim,xDim);

	if (yDim>1)
		globalIntegrator[2] = new GlobalIntegrator<tw::Complex>(&task->strip[2],xDim*zDim,yDim);

	if (zDim>1)
		globalIntegrator[3] = new GlobalIntegrator<tw::Complex>(&task->strip[3],xDim*yDim,zDim);
}

SchroedingerPropagator::~SchroedingerPropagator()
{
	for (auto g : globalIntegrator)
		if (g!=NULL)
			delete g;
}

void SchroedingerPropagator::DepositCurrent(const tw::grid::axis& axis,ComplexField& psi0,ComplexField& psi1,Field& A4,Field& J4,tw::Complex dt)
{
	// for J1,J2,J3, psi0 and psi1 should be the wavefunction before and after the corresponding one dimensional sweep, respectively
	// for J0, this has to be called twice, with psi1 containing the wavefunction at the beginning of the time step on the first call, and at the end on the second call
	// this is the conservative current evaluation from J. Comp. Phys. 280, 457 (2015)

	const tw::Int ax=tw::grid::naxis(axis);

	if (imag(dt)!=0.0)
		return;

	if (axis==tw::grid::t)
	{
		#pragma omp parallel
		{
			for (auto cell : InteriorCellRange(*space))
				J4(cell,0) += half*norm(psi1(cell));
		}
	}
	else
	{
		if (space->Dim(ax)>1)
		{
			#pragma omp parallel
			{
				tw::Complex f11,f12,f21,f22;
				for (auto s : StripRange(*space,ax,strongbool::no))
				{
					for (tw::Int i=1;i<=space->Dim(ax);i++)
					{
						f11 = -half*ii*(psi0(s,i)-psi0(s,i-1))/space->dl(s,i,ax) + tw::Float(0.25)*(A4(s,i-1,ax)+A4(s,i,ax))*psi0(s,i-1);
						f21 = -half*ii*(psi1(s,i)-psi1(s,i-1))/space->dl(s,i,ax) + tw::Float(0.25)*(A4(s,i-1,ax)+A4(s,i,ax))*psi1(s,i-1);
						f12 = -half*ii*(psi0(s,i+1)-psi0(s,i))/space->dl(s,i+1,ax) + tw::Float(0.25)*(A4(s,i,ax)+A4(s,i+1,ax))*psi0(s,i+1);
						f22 = -half*ii*(psi1(s,i+1)-psi1(s,i))/space->dl(s,i+1,ax) + tw::Float(0.25)*(A4(s,i,ax)+A4(s,i+1,ax))*psi1(s,i+1);
						J4(s,i,ax) += real(tw::Float(0.25)*(conj(psi0(s,i))*(f11 + f21) + conj(psi1(s,i))*(f11 + f21)));
						J4(s,i,ax) += real(tw::Float(0.25)*(conj(psi0(s,i))*(f12 + f22) + conj(psi1(s,i))*(f12 + f22)));
					}
				}
			}
		}
	}
}

void SchroedingerPropagator::ApplyNumerator(const tw::grid::axis& axis,ComplexField& psi,Field& A4,bool keepA2Term,tw::Complex dt)
{
	const tw::Int ax = tw::grid::naxis(axis);
	const tw::Float partitionFactor = 1.0 / ((space->Dim(1)>1 ? 1.0 : 0.0) + (space->Dim(2)>1 ? 1.0 : 0.0) + (space->Dim(3)>1 ? 1.0 : 0.0));
	const tw::Float A2Factor = keepA2Term ? 0.5 : 0.0;

	if (space->Dim(ax)>1)
	{
		#pragma omp parallel firstprivate(dt)
		{
			std::valarray<tw::Complex> src(space->Dim(ax));
			for (auto s : StripRange(*space,ax,strongbool::no))
			{
				for (tw::Int i=1;i<=space->Dim(ax);i++)
				{
					const tw::Float Vol = space->dS(s,i,0);
					const tw::Float S1 = space->dS(s,i,ax);
					const tw::Float S2 = space->dS(s,i+1,ax);
					const tw::Float D1 = (S1/Vol) / space->dl(s,i,ax);
					const tw::Float D2 = (S2/Vol) / space->dl(s,i+1,ax);
					const tw::Float A1 = half*(A4(s,i-1,ax) + A4(s,i,ax));
					const tw::Float A2 = half*(A4(s,i,ax) + A4(s,i+1,ax));
					const tw::Float Ueff = partitionFactor*(A2Factor*Norm(A4.Vec3(s,i,1)) - A4(s,i,0));
					const tw::Complex H1 = -half*D1 + half*ii*S1*A1/Vol;
					const tw::Complex H2 = half*(D1+D2) + Ueff;
					const tw::Complex H3 = -half*D2 - half*ii*S2*A2/Vol;

					src[i-1]  = -half*ii*dt*H1*psi(s,i-1);
					src[i-1] += (one - half*ii*dt*H2)*psi(s,i);
					src[i-1] += -half*ii*dt*H3*psi(s,i+1);
				}

				for (tw::Int i=1;i<=space->Dim(ax);i++)
					psi(s,i) = src[i-1];
			}
		}
		psi.CopyFromNeighbors();
		psi.ApplyBoundaryCondition();
	}
}

void SchroedingerPropagator::ApplyDenominator(const tw::grid::axis& axis,ComplexField& psi,Field& A4,bool keepA2Term,tw::Complex dt)
{
	const tw::Int ax = tw::grid::naxis(axis);
	const tw::Int sDim = space->Dim(ax);
	const tw::Float partitionFactor = 1.0 / ((space->Dim(1)>1 ? 1.0 : 0.0) + (space->Dim(2)>1 ? 1.0 : 0.0) + (space->Dim(3)>1 ? 1.0 : 0.0));
	const tw::Float A2Factor = keepA2Term ? 0.5 : 0.0;

	if (sDim>1)
	{
		#pragma omp parallel firstprivate(dt)
		{
			StripRange range(*space,ax,strongbool::no);
			std::valarray<tw::Complex> src(sDim),ans(sDim),T1(sDim),T2(sDim),T3(sDim);
			for (auto it=range.begin();it!=range.end();++it)
			{
				tw::strip s = *it;
				for (tw::Int i=1;i<=sDim;i++)
				{
					const tw::Float Vol = space->dS(s,i,0);
					const tw::Float S1 = space->dS(s,i,ax);
					const tw::Float S2 = space->dS(s,i+1,ax);
					const tw::Float D1 = (S1/Vol) / space->dl(s,i,ax);
					const tw::Float D2 = (S2/Vol) / space->dl(s,i+1,ax);
					const tw::Float A1 = half*(A4(s,i-1,ax) + A4(s,i,ax));
					const tw::Float A2 = half*(A4(s,i,ax) + A4(s,i+1,ax));
					const tw::Float Ueff = partitionFactor*(A2Factor*Norm(A4.Vec3(s,i,1)) - A4(s,i,0));
					const tw::Complex H1 = -half*D1 + half*ii*S1*A1/Vol;
					const tw::Complex H2 = half*(D1+D2) + Ueff;
					const tw::Complex H3 = -half*D2 - half*ii*S2*A2/Vol;

					src[i-1] = psi(s,i);
					T1[i-1] = half*ii*dt*H1;
					T2[i-1] = one + half*ii*dt*H2;
					T3[i-1] = half*ii*dt*H3;
				}

				if (task->n0[1]==MPI_PROC_NULL)
					psi.AdjustTridiagonalForBoundaries(axis,tw::grid::low,T1,T2,T3,src,tw::Complex(0.0));
				if (task->n1[1]==MPI_PROC_NULL)
					psi.AdjustTridiagonalForBoundaries(axis,tw::grid::high,T1,T2,T3,src,tw::Complex(0.0));

				TriDiagonal(ans,src,T1,T2,T3);
				for (tw::Int i=1;i<=sDim;i++)
					psi(s,i) = ans[i-1];

				globalIntegrator[ax]->SetMatrix(it.global_count(),T1,T2,T3);
				globalIntegrator[ax]->SetData(it.global_count(),&psi(s,0),psi.Stride(ax)/2); // need complex stride and stride[0]=1
			}
		}

		globalIntegrator[ax]->Parallelize();
		psi.CopyFromNeighbors();
		psi.ApplyBoundaryCondition();
	}
}

void SchroedingerPropagator::UpdateSpin(ComplexField& psi,ComplexField& chi,Field& A4,tw::Float adt)
{
	tw::Int i,j,k;
	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);
	tw::Complex temp,U11,U12,U21,U22;
	tw::Float Bx,By,Bz,B2,denom;
	for (k=1;k<=zDim;k++)
		for (j=1;j<=yDim;j++)
			for (i=1;i<=xDim;i++)
			{
				Bx = (A4(i,j+1,k,3)-A4(i,j-1,k,3))/space->dL(i,j,k,2);
				Bx -= (A4(i,j,k+1,2)-A4(i,j,k-1,2))/space->dL(i,j,k,3);
				By = (A4(i,j,k+1,1)-A4(i,j,k-1,1))/space->dL(i,j,k,3);
				By -= (A4(i+1,j,k,3)-A4(i-1,j,k,3))/space->dL(i,j,k,1);
				Bz = (A4(i+1,j,k,2)-A4(i-1,j,k,2))/space->dL(i,j,k,1);
				Bz -= (A4(i,j+1,k,1)-A4(i,j-1,k,1))/space->dL(i,j,k,2);

				B2 = Bx*Bx + By*By + Bz*Bz;
				denom = tw::Float(16)+B2*sqr(adt);
				U11 = (tw::Float(16)+adt*(tw::Float(8)*ii*Bz-B2*adt))/denom;
				U12 = tw::Float(8)*(ii*Bx+By)*adt/denom;
				U21 = tw::Float(8)*ii*(Bx+ii*By)*adt/denom;
				U22 = (tw::Float(16)+adt*(-tw::Float(8)*ii*Bz-B2*adt))/denom;

				temp = psi(i,j,k);
				psi(i,j,k) = U11*psi(i,j,k) + U12*chi(i,j,k);
				chi(i,j,k) = U21*temp + U22*chi(i,j,k);
			}
}


//////////////////////////////////////////
//                                      //
//    Solver for Parabolic Equations    //
//                                      //
//////////////////////////////////////////


ParabolicSolver::ParabolicSolver(const std::string& name,MetricSpace *m,Task *tsk) : BoundedTool(name,m,tsk)
{
	globalIntegrator.assign(4,NULL);

	if (space->Dim(1)>1)
		globalIntegrator[1] = new GlobalIntegrator<tw::Float>(&task->strip[1],space->Dim(2)*space->Dim(3),space->Dim(1));

	if (space->Dim(2)>1)
		globalIntegrator[2] = new GlobalIntegrator<tw::Float>(&task->strip[2],space->Dim(3)*space->Dim(1),space->Dim(2));

	if (space->Dim(3)>1)
		globalIntegrator[3] = new GlobalIntegrator<tw::Float>(&task->strip[3],space->Dim(1)*space->Dim(2),space->Dim(3));
}

ParabolicSolver::~ParabolicSolver()
{
	for (auto g : globalIntegrator)
		if (g!=NULL)
			delete g;
}

void ParabolicSolver::FixTemperature(Field& T,const Element& e,Region* theRegion,const tw::Float& T0)
{
	#pragma omp parallel
	{
		for (auto cell : EntireCellRange(*space))
			if (theRegion->Inside(space->Pos(cell),*space))
				for (tw::Int c=e.low;c<=e.high;c++)
					T(cell,c) = T0;
	}
}

inline void ParabolicSolver::FormOperatorStencil(tw::Float *D1,tw::Float *D2,const ScalarField& fluxMask,Field *coeff,tw::Int c,const tw::strip& s,tw::Int i)
{
	tw::Float dV = space->dS(s,i,0);
	tw::Int ax = s.Axis();
	*D1 = (space->dS(s,i,ax)/dV)/space->dl(s,i,ax);
	*D2 = (space->dS(s,i+1,ax)/dV)/space->dl(s,i+1,ax);
	if (coeff!=NULL)
	{
		*D1 *= 0.5*((*coeff)(s,i-1,c) + (*coeff)(s,i,c));
		*D2 *= 0.5*((*coeff)(s,i,c) + (*coeff)(s,i+1,c));
	}
	*D1 *= fluxMask(s,i-1,0)*fluxMask(s,i,0);
	*D2 *= fluxMask(s,i,0)*fluxMask(s,i+1,0);
}

void ParabolicSolver::Advance(const tw::grid::axis& axis,ScalarField& psi,ScalarField& fluxMask,tw::Float coeff,tw::Float dt)
{
	// solve d(psi)/dt - coeff*div(grad(psi)) = 0
	// solve on strips running parallel to axis

	const tw::Int ax=tw::grid::naxis(axis);
	const tw::Int sDim=space->Dim(ax);

	#pragma omp parallel
	{
		StripRange range(*space,ax,strongbool::no);
		tw::Float D1,D2;
		std::valarray<tw::Float> src,ans,T1,T2,T3;

		ans.resize(sDim);
		src.resize(sDim);
		T1.resize(sDim);
		T2.resize(sDim);
		T3.resize(sDim);
		for (auto it=range.begin();it!=range.end();++it)
		{
			tw::strip s = *it;
			for (tw::Int i=1;i<=sDim;i++)
			{
				FormOperatorStencil(&D1,&D2,fluxMask,NULL,0,s,i);
				src[i-1] = psi(s,i);
				T1[i-1] = -dt*coeff*D1;
				T2[i-1] = 1.0 + dt*coeff*(D1+D2);
				T3[i-1] = -dt*coeff*D1;
			}

			if (task->n0[ax]==MPI_PROC_NULL)
				psi.AdjustTridiagonalForBoundaries(axis,tw::grid::low,T1,T2,T3,src,psi(s,space->LFG(ax)));
			if (task->n1[ax]==MPI_PROC_NULL)
				psi.AdjustTridiagonalForBoundaries(axis,tw::grid::high,T1,T2,T3,src,psi(s,space->UFG(ax)));

			TriDiagonal(ans,src,T1,T2,T3);
			for (tw::Int i=1;i<=sDim;i++)
				psi(s,i) = ans[i-1];

			globalIntegrator[ax]->SetMatrix(it.global_count(),T1,T2,T3);
			globalIntegrator[ax]->SetData(it.global_count(),&psi(s,0),psi.Stride(ax));
		}
	}

	globalIntegrator[ax]->Parallelize();
	psi.ApplyBoundaryCondition(false);
	// Leaving ghost cells in orthogonal directions unspecified.
	// They will be set when those axes are propagated.
}

void ParabolicSolver::Advance(ScalarField& psi,ScalarField& fluxMask,tw::Float coeff,tw::Float dt)
{
	// solve d(psi)/dt - coeff*div(grad(psi)) = 0

	if (space->Dim(1)>1)
		Advance(tw::grid::x,psi,fluxMask,coeff,dt);

	if (space->Dim(2)>1)
		Advance(tw::grid::y,psi,fluxMask,coeff,dt);

	if (space->Dim(3)>1)
		Advance(tw::grid::z,psi,fluxMask,coeff,dt);
}

void ParabolicSolver::Advance(	const tw::grid::axis& axis,
								Field& psi,
								tw::Int psi_idx,
								ScalarField& fluxMask,
								Field *coeff1,
								tw::Int c1_idx,
								Field *coeff2,
								tw::Int c2_idx,
								tw::Float dt)
{
	// solve coeff1*d(psi)/dt - div(coeff2*grad(psi)) = 0
	// solve on strips running parallel to axis
	// in hydro operator splitting context, putting coeff1 inside time derivative would be incorrect

	const tw::Int ax=tw::grid::naxis(axis);
	const tw::Int sDim=space->Dim(ax);

	#pragma omp parallel
	{
		StripRange range(*space,ax,strongbool::no);
		tw::Float D1,D2;
		std::valarray<tw::Float> src,ans,T1,T2,T3;

		ans.resize(sDim);
		src.resize(sDim);
		T1.resize(sDim);
		T2.resize(sDim);
		T3.resize(sDim);
		for (auto it=range.begin();it!=range.end();++it)
		{
			tw::strip s = *it;
			for (tw::Int i=1;i<=sDim;i++)
			{
				FormOperatorStencil(&D1,&D2,fluxMask,coeff2,c2_idx,s,i);
				src[i-1] = (*coeff1)(s,i,c1_idx)*psi(s,i,psi_idx);
				T1[i-1] = -dt*D1;
				T2[i-1] = (*coeff1)(s,i,c1_idx) + dt*(D1+D2);
				T3[i-1] = -dt*D2;
			}

			if (task->n0[ax]==MPI_PROC_NULL)
				psi.AdjustTridiagonalForBoundaries(Element(psi_idx),axis,tw::grid::low,T1,T2,T3,src,psi(s,space->LFG(ax),psi_idx));
			if (task->n1[ax]==MPI_PROC_NULL)
				psi.AdjustTridiagonalForBoundaries(Element(psi_idx),axis,tw::grid::high,T1,T2,T3,src,psi(s,space->UFG(ax),psi_idx));

			TriDiagonal(ans,src,T1,T2,T3);
			for (tw::Int i=1;i<=sDim;i++)
				psi(s,i,psi_idx) = ans[i-1];

			globalIntegrator[ax]->SetMatrix(it.global_count(),T1,T2,T3);
			globalIntegrator[ax]->SetData(it.global_count(),&psi(s,0,psi_idx),psi.Stride(ax));
		}
	}

	globalIntegrator[ax]->Parallelize();
	psi.ApplyBoundaryCondition(Element(psi_idx),false);
	// Leaving ghost cells in orthogonal directions unspecified.
	// They will be set when those axes are propagated.
}

void ParabolicSolver::Advance(	Field& psi,
								tw::Int psi_idx,
								ScalarField& fluxMask,
								Field *coeff1,
								tw::Int c1_idx,
								Field *coeff2,
								tw::Int c2_idx,
								tw::Float dt)
{
	// solve coeff1*d(psi)/dt - div(coeff2*grad(psi)) = 0

	if (space->Dim(1)>1)
		Advance(tw::grid::x,psi,psi_idx,fluxMask,coeff1,c1_idx,coeff2,c2_idx,dt);

	if (space->Dim(2)>1)
		Advance(tw::grid::y,psi,psi_idx,fluxMask,coeff1,c1_idx,coeff2,c2_idx,dt);

	if (space->Dim(3)>1)
		Advance(tw::grid::z,psi,psi_idx,fluxMask,coeff1,c1_idx,coeff2,c2_idx,dt);
}


////////////////////////////////////////////
//                                        //
//      ISOTROPIC WAVE PROPAGATOR         //
//                                        //
////////////////////////////////////////////


IsotropicPropagator::IsotropicPropagator(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk)
{
	zGlobalIntegrator = NULL;

	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	Z1.resize(zDim);
	Z2.resize(zDim);
	Z3.resize(zDim);

	if (zDim>1)
		zGlobalIntegrator = new GlobalIntegrator<tw::Complex>(&task->strip[3],xDim*yDim,zDim);
}

IsotropicPropagator::~IsotropicPropagator()
{
	if (zGlobalIntegrator!=NULL)
		delete zGlobalIntegrator;
}

void IsotropicPropagator::Advance(
	ComplexField& amplitude, // this is the array to update
	ComplexField& refractiveIndex, // index of refraction array
	ScalarField& nu_e, // electron collision frequency array
	tw::Float laserFrequency,tw::Float dt)
{
	// Solve d/dt(a) = (0.5i/w0) ((d/dz)^2(a) + w0^2 n^2 a)
	// Source functions for incoming waves are assumed to be loaded into ghost cells of "amplitude"

	const tw::Int xDim = space->Dim(1);
	const tw::Int yDim = space->Dim(2);
	const tw::Int zDim = space->Dim(3);

	tw::Int i,j,k,idx;
	tw::Float S1,S2,Vol,D1,D2;
	tw::Complex H1,H2,H3; // related to hamiltonian
	tw::Complex nw,eta_p,eta_m;
	std::valarray<tw::Complex> src,ans;

	const tw::Float dz = space->dl(1,1,1,3); // assume uniform grid

	src.resize(zDim);
	ans.resize(zDim);

	eta_p = one+half*ii*laserFrequency*dz;
	eta_m = one-half*ii*laserFrequency*dz;

	idx = 0;
	for (j=1;j<=yDim;j++)
		for (i=1;i<=xDim;i++)
		{
			for (k=1;k<=zDim;k++)
			{
				Vol = space->dS(i,j,k,0);
				S1 = space->dS(i,j,k,3);
				S2 = space->dS(i,j,k+1,3);
				D1 = (S1/Vol) / space->dl(i,j,k,3);
				D2 = (S2/Vol) / space->dl(i,j,k+1,3);
				nw = laserFrequency*refractiveIndex(i,j,k);

				H1 = tw::Float(0.25)*ii*dt*D1/laserFrequency;
				H2 = tw::Float(0.25)*ii*dt*(nw*nw - D1 - D2)/laserFrequency;
				H3 = tw::Float(0.25)*ii*dt*D2/laserFrequency;

				// Use this block for fully implicit advance
				src[k-1] = amplitude(i,j,k);
				Z1[k-1] = -two*H1;
				Z2[k-1] = one-two*H2;
				Z3[k-1] = -two*H3;

				// Use this block for semi implicit advance
				//src[k-1] = H1*amplitude(i,j,k-1) + (one+H2)*amplitude(i,j,k) + H3*amplitude(i,j,k+1);
				//Z1[k-1] = -H1;
				//Z2[k-1] = one-H2;
				//Z3[k-1] = -H3;
			}

			if (task->n0[3]==MPI_PROC_NULL)
			{
				Z2[0] += Z1[0]*eta_p/eta_m;
				src[0] += Z1[0]*amplitude(i,j,0)/eta_m;
			}
			if (task->n1[3]==MPI_PROC_NULL)
			{
				Z2[zDim-1] += Z3[zDim-1]*eta_p/eta_m;
				src[zDim-1] -= Z3[zDim-1]*amplitude(i,j,zDim+1)/eta_m;
			}

			TriDiagonal(ans,src,Z1,Z2,Z3);
			for (k=1;k<=zDim;k++)
				amplitude(i,j,k) = ans[k-1];

			zGlobalIntegrator->SetMatrix(idx,Z1,Z2,Z3);
			zGlobalIntegrator->SetData(idx,&amplitude(i,j,0),amplitude.Stride(3)/2); // need complex stride and stride[0]=1
			idx++;
		}

	zGlobalIntegrator->Parallelize();
	amplitude.CopyFromNeighbors();
//	amplitude.ApplyBoundaryCondition();
	for (j=1;j<=yDim;j++)
		for (i=1;i<=xDim;i++)
		{
			if (task->n0[3]==MPI_PROC_NULL)
				amplitude(i,j,0) = amplitude(i,j,1)*eta_p/eta_m - amplitude(i,j,0)/eta_m;
			if (task->n1[3]==MPI_PROC_NULL)
				amplitude(i,j,zDim+1) = amplitude(i,j,zDim)*eta_p/eta_m + amplitude(i,j,zDim+1)/eta_m;
		}
}
