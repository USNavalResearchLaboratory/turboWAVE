module;

#include "tw_includes.h"
#include "tw_logger.h"

export module parabolic:quantum;
import input;
import driver;
import fields;
import numerics;
import logger;

using namespace tw::bc;
using SharedRegion = std::shared_ptr<Region>;

/// Uses a probability conserving algorithm to advance the quantum mechanical
/// Schroedinger equation in configuration space.
export struct SchroedingerPropagator:ComputeTool
{
	std::unique_ptr<GlobalIntegrator<tw::Complex>> globalIntegrator[4];

	SchroedingerPropagator(const std::string& name,MetricSpace *m,Task *tsk): ComputeTool(name,m,tsk) {
        const tw::Int xDim = space->Dim(1);
        const tw::Int yDim = space->Dim(2);
        const tw::Int zDim = space->Dim(3);

        if (xDim>1)
            globalIntegrator[1] = std::make_unique<GlobalIntegrator<tw::Complex>>(&task->strip[1],yDim*zDim,xDim);

        if (yDim>1)
            globalIntegrator[2] = std::make_unique<GlobalIntegrator<tw::Complex>>(&task->strip[2],xDim*zDim,yDim);

        if (zDim>1)
            globalIntegrator[3] = std::make_unique<GlobalIntegrator<tw::Complex>>(&task->strip[3],xDim*yDim,zDim);

        logger::TRACE("end constructor");
    }
	virtual void DepositCurrent(const tw::grid::axis& axis,ComplexField& psi0,ComplexField& psi1,Field& A4,Field& J4,tw::Complex dt) {
        // for J1,J2,J3, psi0 and psi1 should be the wavefunction before and after the corresponding one dimensional sweep, respectively
        // for J0, this has to be called twice, with psi1 containing the wavefunction at the beginning of the time step on the first call, and at the end on the second call
        // this is the conservative current evaluation from J. Comp. Phys. 280, 457 (2015)

        const tw::Int ax=tw::grid::naxis(axis);
        logger::TRACE(std::format("deposit J{}",ax));
        
        if (imag(dt)!=0.0)
            return;

        if (axis==tw::grid::t)
        {
            #pragma omp parallel
            {
                for (auto cell : InteriorCellRange(*space,1))
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
                    for (auto s : StripRange(*space,ax,0,1,strongbool::no))
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
	virtual void ApplyNumerator(const tw::grid::axis& axis,ComplexField& psi,Field& A4,bool keepA2Term,tw::Complex dt) {
        const tw::Int ax = tw::grid::naxis(axis);
        const tw::Float partitionFactor = 1.0 / ((space->Dim(1)>1 ? 1.0 : 0.0) + (space->Dim(2)>1 ? 1.0 : 0.0) + (space->Dim(3)>1 ? 1.0 : 0.0));
        const tw::Float A2Factor = keepA2Term ? 0.5 : 0.0;

        logger::TRACE(std::format("apply numerator {}",ax));

        if (space->Dim(ax)>1)
        {
            #pragma omp parallel firstprivate(dt)
            {
                tw::vec<tw::Complex> src(space->Dim(ax));
                for (auto s : StripRange(*space,ax,0,1,strongbool::no))
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
                        psi.Pack(s, i, src[i-1]);
                }
            }
            psi.CopyFromNeighbors();
            psi.ApplyBoundaryCondition();
        }
	}
	virtual void ApplyDenominator(const tw::grid::axis& axis,ComplexField& psi,Field& A4,bool keepA2Term,tw::Complex dt) {
        const tw::Int ax = tw::grid::naxis(axis);
        const tw::Int sDim = space->Dim(ax);
        const tw::Float partitionFactor = 1.0 / ((space->Dim(1)>1 ? 1.0 : 0.0) + (space->Dim(2)>1 ? 1.0 : 0.0) + (space->Dim(3)>1 ? 1.0 : 0.0));
        const tw::Float A2Factor = keepA2Term ? 0.5 : 0.0;

        logger::TRACE(std::format("apply denominator {}",ax));

        if (sDim>1)
        {
            #pragma omp parallel firstprivate(dt)
            {
                StripRange range(*space,ax,0,1,strongbool::no);
                tw::vec<tw::Complex> src(sDim),ans(sDim),T1(sDim),T2(sDim),T3(sDim);
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

                    if (task->n0[ax]==MPI_PROC_NULL)
                        psi.AdjustTridiagonalForBoundaries(axis,tw::grid::low,T1,T2,T3,src,tw::Complex(0.0));
                    if (task->n1[ax]==MPI_PROC_NULL)
                        psi.AdjustTridiagonalForBoundaries(axis,tw::grid::high,T1,T2,T3,src,tw::Complex(0.0));

                    TriDiagonal(ans,src,T1,T2,T3);
                    for (tw::Int i=1;i<=sDim;i++)
                        psi.Pack(s, i, ans[i-1]);

                    globalIntegrator[ax]->SetData(it.global_count(),&psi(s,0,0),psi.Stride(ax),psi.Stride(4));
                    globalIntegrator[ax]->SetMatrix(it.global_count(),T1,T2,T3);
                }
            }

            globalIntegrator[ax]->Parallelize();
            psi.CopyFromNeighbors();
            psi.ApplyBoundaryCondition();
        }
    }
	virtual void UpdateSpin(ComplexField& psi,ComplexField& chi,Field& A4,tw::Float adt) {
        const tw::Int xDim = space->Dim(1);
        const tw::Int yDim = space->Dim(2);
        const tw::Int zDim = space->Dim(3);
        tw::Complex temp,U11,U12,U21,U22;
        tw::Float Bx,By,Bz,B2,denom;
        for (auto k=1;k<=zDim;k++) {
            for (auto j=1;j<=yDim;j++) {
                for (auto i=1;i<=xDim;i++) {
                    Bx = (A4(1,i,j+1,k,3)-A4(1,i,j-1,k,3))/space->dL(i,j,k,2);
                    Bx -= (A4(1,i,j,k+1,2)-A4(1,i,j,k-1,2))/space->dL(i,j,k,3);
                    By = (A4(1,i,j,k+1,1)-A4(1,i,j,k-1,1))/space->dL(i,j,k,3);
                    By -= (A4(1,i+1,j,k,3)-A4(1,i-1,j,k,3))/space->dL(i,j,k,1);
                    Bz = (A4(1,i+1,j,k,2)-A4(1,i-1,j,k,2))/space->dL(i,j,k,1);
                    Bz -= (A4(1,i,j+1,k,1)-A4(1,i,j-1,k,1))/space->dL(i,j,k,2);

                    B2 = Bx*Bx + By*By + Bz*Bz;
                    denom = tw::Float(16)+B2*sqr(adt);
                    U11 = (tw::Float(16)+adt*(tw::Float(8)*ii*Bz-B2*adt))/denom;
                    U12 = tw::Float(8)*(ii*Bx+By)*adt/denom;
                    U21 = tw::Float(8)*ii*(Bx+ii*By)*adt/denom;
                    U22 = (tw::Float(16)+adt*(-tw::Float(8)*ii*Bz-B2*adt))/denom;

                    temp = psi(i,j,k);
                    psi.Pack(i,j,k,U11*psi(i,j,k) + U12*chi(i,j,k));
                    chi.Pack(i,j,k,U21*temp + U22*chi(i,j,k));
                }
            }
        }
    }
};
