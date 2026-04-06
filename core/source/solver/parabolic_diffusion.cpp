module;

#include "tw_includes.h"
#include "tw_logger.h"

export module parabolic:diffusion;
import input;
import driver;
import fields;
import numerics;
import logger;

using namespace tw::bc;
using SharedRegion = std::shared_ptr<Region>;

export struct ParabolicSolver:BoundedTool
{
	// inhomogeneous boundary conditions are handled as described in elliptic.h
	std::unique_ptr<GlobalIntegrator<tw::Float>> globalIntegrator[4]; 

	ParabolicSolver(const std::string& name,MetricSpace *m,Task *tsk): BoundedTool(name,m,tsk) {
        if (space->Dim(1)>1)
            globalIntegrator[1] = std::make_unique<GlobalIntegrator<tw::Float>>(&task->strip[1],space->Dim(2)*space->Dim(3),space->Dim(1));

        if (space->Dim(2)>1)
            globalIntegrator[2] = std::make_unique<GlobalIntegrator<tw::Float>>(&task->strip[2],space->Dim(3)*space->Dim(1),space->Dim(2));

        if (space->Dim(3)>1)
            globalIntegrator[3] = std::make_unique<GlobalIntegrator<tw::Float>>(&task->strip[3],space->Dim(1)*space->Dim(2),space->Dim(3));
    }

	void FixTemperature(Field& T,const Rng04& r,SharedRegion theRegion,const tw::Float& T0) {
        #pragma omp parallel
        {
            for (auto n=r.b0; n<r.e0; n++)
            for (auto cell : EntireCellRange(*space,n))
                if (theRegion->Inside(space->Pos4(cell),0))
                    for (tw::Int c=r.b4; c<=r.e4; c++)
                        T(cell,c) = T0;
        }
    }
	void FormOperatorStencil(tw::Float *D1,tw::Float *D2,const ScalarField& fluxMask,Field *coeff,tw::Int c,const tw::strip& s,tw::Int i) {
        tw::Float dV = space->dS(s,i,0);
        tw::Int ax = s.StripAxis();
        *D1 = (space->dS(s,i,ax)/dV)/space->dl(s,i,ax);
        *D2 = (space->dS(s,i+1,ax)/dV)/space->dl(s,i+1,ax);
        if (coeff!=NULL)
        {
            *D1 *= 0.5*((*coeff)(s,i-1,c) + (*coeff)(s,i,c));
            *D2 *= 0.5*((*coeff)(s,i,c) + (*coeff)(s,i+1,c));
        }
        *D1 *= fluxMask(s,i-1)*fluxMask(s,i);
        *D2 *= fluxMask(s,i)*fluxMask(s,i+1);
    }
    /// solve d(psi)/dt - coeff*div(grad(psi)) = 0
    /// solve on strips running parallel to axis
	void Advance(const tw::grid::axis& axis,ScalarField& psi,ScalarField& fluxMask,tw::Float coeff,tw::Float dt) {

        const tw::Int ax=tw::grid::naxis(axis);
        const tw::Int sDim=space->Dim(ax);

        #pragma omp parallel
        {
            StripRange range(*space,ax,0,1,strongbool::no);
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
                globalIntegrator[ax]->SetData(it.global_count(),&psi(s,0),psi.Stride(ax),psi.Stride(4));
            }
        }

        globalIntegrator[ax]->Parallelize();
        psi.ApplyBoundaryCondition(false);
        // Leaving ghost cells in orthogonal directions unspecified.
        // They will be set when those axes are propagated.
    }
	/// solve d(psi)/dt - coeff*div(grad(psi)) = 0
	void Advance(ScalarField& psi,ScalarField& fluxMask,tw::Float coeff,tw::Float dt) {
        if (space->Dim(1)>1)
            Advance(tw::grid::x,psi,fluxMask,coeff,dt);

        if (space->Dim(2)>1)
            Advance(tw::grid::y,psi,fluxMask,coeff,dt);

        if (space->Dim(3)>1)
            Advance(tw::grid::z,psi,fluxMask,coeff,dt);
    }
	/// solve coeff1*d(psi)/dt - div(coeff2*grad(psi)) = 0
	/// solve on strips running parallel to axis
	/// in hydro operator splitting context, putting coeff1 inside time derivative would be incorrect
	void Advance(	const tw::grid::axis& axis,
							Field& psi,
							tw::Int psi_idx,
							ScalarField& fluxMask,
							Field *coeff1,
							tw::Int c1_idx,
							Field *coeff2,
							tw::Int c2_idx,
							tw::Float dt)
    {

        const tw::Int ax=tw::grid::naxis(axis);
        const tw::Int sDim=space->Dim(ax);

        #pragma omp parallel
        {
            StripRange range(*space,ax,0,1,strongbool::no);
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
                    psi.AdjustTridiagonalForBoundaries(Rng(psi_idx),axis,tw::grid::low,T1,T2,T3,src,psi(s,space->LFG(ax),psi_idx));
                if (task->n1[ax]==MPI_PROC_NULL)
                    psi.AdjustTridiagonalForBoundaries(Rng(psi_idx),axis,tw::grid::high,T1,T2,T3,src,psi(s,space->UFG(ax),psi_idx));

                TriDiagonal(ans,src,T1,T2,T3);
                for (tw::Int i=1;i<=sDim;i++)
                    psi(s,i,psi_idx) = ans[i-1];

                globalIntegrator[ax]->SetMatrix(it.global_count(),T1,T2,T3);
                globalIntegrator[ax]->SetData(it.global_count(),&psi(s,0,psi_idx),psi.Stride(ax),psi.Stride(4));
            }
        }

        globalIntegrator[ax]->Parallelize();
        psi.ApplyBoundaryCondition(Rng(psi_idx),false);
        // Leaving ghost cells in orthogonal directions unspecified.
        // They will be set when those axes are propagated.

    }
	/// solve coeff1*d(psi)/dt - div(coeff2*grad(psi)) = 0
	virtual void Advance(	Field& psi,
							tw::Int psi_idx,
							ScalarField& fluxMask,
							Field *coeff1,
							tw::Int c1_idx,
							Field *coeff2,
							tw::Int c2_idx,
							tw::Float dt)
    {
        if (space->Dim(1)>1)
            Advance(tw::grid::x,psi,psi_idx,fluxMask,coeff1,c1_idx,coeff2,c2_idx,dt);

        if (space->Dim(2)>1)
            Advance(tw::grid::y,psi,psi_idx,fluxMask,coeff1,c1_idx,coeff2,c2_idx,dt);

        if (space->Dim(3)>1)
            Advance(tw::grid::z,psi,psi_idx,fluxMask,coeff1,c1_idx,coeff2,c2_idx,dt);

        }
};
