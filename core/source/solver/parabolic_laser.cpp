module;

#include "tw_includes.h"
#include "tw_logger.h"

export module parabolic:laser;
import input;
import driver;
import fields;
import numerics;
import logger;

using namespace tw::bc;
using SharedRegion = std::shared_ptr<Region>;

export enum tw_polarization_type {linearPolarization,circularPolarization,radialPolarization};

/// Propagator to advance coherent radiation unidirectionally.
/// Take a~ = (a/2)exp(i(kx-wt)) + cc and solve (delperp^2 + 2iw0Dt + 2Dzt)a = -chi*a.
/// The 2Dzt term uses a second order forward difference [W. Zhu et al., Phys. Plasmas 19, 033105 (2012)].
/// This is designed to work for either Cartesian or axisymmetric coordinates.
export struct ForwardPropagator:ComputeTool
{
	tw::Float w0,dt;
	tw_polarization_type polarization;
	bool evenTime,movingWindow;
	std::unique_ptr<GlobalIntegrator<tw::Complex>> globalIntegrator[4];
	ComplexField aNow;

	ForwardPropagator(const std::string& name,MetricSpace *m,Task *tsk): ComputeTool(name,m,tsk) {
		w0 = 10.0;
		polarization = linearPolarization;
		movingWindow = true;
		evenTime = true;
	}
	virtual void SetData(tw::Float w0,tw::Float dt,tw_polarization_type pol,bool mov,MetricSpace *refined = NULL) {
		this->w0 = w0;
		this->dt = dt;
		polarization = pol;
		movingWindow = mov;
		if (refined != NULL) {
			space = refined;
		}

		aNow.Initialize(*space,task);

		tw::Int systems = space->Dim(2);
		tw::Int cells = space->Dim(1);
		globalIntegrator[1] = std::make_unique<GlobalIntegrator<tw::Complex>>(&task->strip[1],systems,cells);

		systems = space->Dim(1);
		cells = space->Dim(2);
		globalIntegrator[2] = std::make_unique<GlobalIntegrator<tw::Complex>>(&task->strip[2],systems,cells);

	}
	void SetBoundaryConditions(ComplexField& a,ComplexField& chi) {
		tw::bc::fld xl,xh,yl,yh;

		xl = space->cyl==1.0 && polarization==radialPolarization ? fld::dirichletWall : fld::neumannWall;
		xh = fld::neumannWall;
		yl = fld::neumannWall;
		yh = fld::neumannWall;

		a.SetBoundaryConditions(tw::grid::x,xl,xh);
		a.SetBoundaryConditions(tw::grid::y,yl,yh);

		xl = space->cyl==1.0 ? fld::neumannWall : fld::dirichletWall;

		chi.SetBoundaryConditions(tw::grid::x,xl,fld::dirichletCell);
		chi.SetBoundaryConditions(tw::grid::y,fld::dirichletCell,fld::dirichletCell);
		chi.SetBoundaryConditions(tw::grid::z,fld::dirichletCell,fld::dirichletCell);
	}
	/// @brief ADI advance with `ua` being the implicitly treated axis (1 or 2)
	void AdvanceAxis(tw::Int ua,ComplexField& a,ComplexField& chi)
	{
		const tw::Int va = ua == 1 ? 2 : 1; // explicit axis
		const tw::Int uDim = space->Dim(ua);
		const tw::Int vDim = space->Dim(va);
		const tw::Int wDim = space->Dim(3);
		const tw::Float radialPolarizationFactor = space->cyl==1.0 && polarization==radialPolarization ? 1.0 : 0.0;
		tw::vec<tw::Complex> src(uDim),ans(uDim),T1(uDim),T2(uDim),T3(uDim);
		for (auto k=wDim;k>=-1;k--) {
			for (auto j=1;j<=vDim;j++) {
				// we have to setup a next (n=1), current (n=2), and previous (n=3) strip
				auto coord = ua == 1 ? tw::node4 {1,0,j,k} : tw::node4 {1,j,0,k};
				const auto un = tw::strip(a,ua,0,coord);
				const auto uc = tw::strip(un,2);
				const auto up = tw::strip(un,3);
				#pragma omp parallel for
				for (auto i=1;i<=uDim;i++) {
					coord = ua == 1 ? tw::node4 {1,i,0,k} : tw::node4 {1,0,i,k};
					const auto vn = tw::strip(a,va,0,coord);
					const auto vp = tw::strip(vn,3);
					coord = ua == 1 ? tw::node4 {1,i,j,0} : tw::node4 {1,j,i,0};
					const auto wn = tw::strip(a,3,0,coord);
					const auto wp = tw::strip(wn,3);

					const auto Vol = space->dS(un,i,0);
					auto A0 = space->dS(un,i,ua);
					auto A1 = space->dS(un,i+1,ua);
					const auto D1 = (A0/Vol) / space->dl(un,i,ua);
					const auto D2 = (A1/Vol) / space->dl(un,i+1,ua);
					const auto dtidzi = 1/space->dl(wn,k+1,3)/dt;
					const auto R = radialPolarizationFactor / sqr(space->X(i,ua));
					const auto T2p = (2.0*ii*w0/dt) - 3.0*dtidzi - (D1+D2+R) + chi(uc,i); 
					const auto T2m = (2.0*ii*w0/dt) - 3.0*dtidzi + (D1+D2+R) - chi(uc,i);
					const auto lookahead = dtidzi*(4.0*a(wn,k+1) - 4.0*a(wp,k+1) + a(wp,k+2) - a(wn,k+2));
					src[i-1] = T2m*a(up,i) - lookahead - D1*a(up,i-1) - D2*a(up,i+1);

					A0 = space->dS(vn,j,va);
					A1 = space->dS(vn,j+1,va);
					const auto D3 = (A0/Vol) / space->dl(vn,j,va);
					const auto D4 = (A1/Vol) / space->dl(vn,j+1,va);
					src[i-1] -= 2.0*(D3*a(vp,j-1) - (D3+D4)*a(vp,j) + D4*a(vp,j+1));

					T1[i-1] = D1;
					T2[i-1] = T2p;
					T3[i-1] = D2;
				}

				if (task->n0[ua]==MPI_PROC_NULL)
					a.AdjustTridiagonalForBoundaries(tw::grid::enumaxis(ua),tw::grid::low,T1,T2,T3,src,tw::Complex(0.0));
				if (task->n1[ua]==MPI_PROC_NULL)
					a.AdjustTridiagonalForBoundaries(tw::grid::enumaxis(ua),tw::grid::high,T1,T2,T3,src,tw::Complex(0.0));

				TriDiagonal(ans,src,T1,T2,T3);
				for (auto i=1;i<=uDim;i++) {
					a.Pack(un, i, ans[i-1]);
				}

				const auto system = j - 1;
				globalIntegrator[ua]->SetMatrix(system,T1,T2,T3);
				globalIntegrator[ua]->SetData(system,&a(un,0,0),a.Stride(ua),a.Stride(4));
			}
			globalIntegrator[ua]->Parallelize();
		}

		a.Field::UpwardCopy(Rng04(1,2,0,2),tw::grid::enumaxis(va),1);
		a.Field::DownwardCopy(Rng04(1,2,0,2),tw::grid::enumaxis(va),1);
		a.Field::ApplyBoundaryCondition(Rng04(1,2,0,2));
	}
	void Advance(ComplexField& a,ComplexField& chi)
	{
		// TODO: the shift has replaced copying current data to a scratch array (aNow),
		// need to finish implementing this change.
		a.ShiftTimeLevels();
		const tw::Int xDim = space->Dim(1);
		const tw::Int yDim = space->Dim(2);
		const tw::Int zDim = space->Dim(3);

		Slice<tw::Float> send_lookahead(
			{1,space->LFG(1),space->LFG(2),1,0},
			{2,space->UFG(1)+1,space->UFG(2)+1,3,2}
		);
		Slice<tw::Float> recv_lookahead(
			{1,space->LFG(1),space->LFG(2),zDim+1,0},
			{2,space->UFG(1)+1,space->UFG(2)+1,zDim+3,2}
		);
		
		tw::Int n0,n1;
		task->strip[3].Shift(1,1,&n0,&n1);

		// We have to handle axial strips serially.
		// Downstream nodes will have to wait for their upstream neighbor before starting.
		if (n1!=MPI_PROC_NULL) {
			task->strip[3].Recv(recv_lookahead.Buffer(),recv_lookahead.BufferSize(),n1);
			a.SaveDataFromSlice(&recv_lookahead);
		}

		if (xDim==1 && yDim==1)
		{
			for (auto k=zDim;k>=-1;k--) {
				const auto dtidzi = 1/space->dl(1,1,k+1,3)/dt;
				const auto T2p = (2.0*ii*w0/dt) - 3.0*dtidzi + chi(1,1,k); 
				const auto T2m = (2.0*ii*w0/dt) - 3.0*dtidzi - chi(1,1,k);
				const auto lookahead = dtidzi*(4.0*a(1,1,k+1) - 4.0*a(3,1,1,k+1) + a(3,1,1,k+2) - a(1,1,k+2));
				a.Pack(1,1,k, (T2m*a(3,1,1,k) - lookahead) / T2p);
			}
		} else if ((!evenTime && xDim>1) || (yDim==1 && xDim>1)) {
			AdvanceAxis(1,a,chi);
		} else if ((evenTime && yDim>1) || (xDim==1 && yDim>1)) {
			AdvanceAxis(2,a,chi);
		}

		if (n0!=MPI_PROC_NULL) {
			a.LoadDataIntoSlice(&send_lookahead);
			task->strip[3].Send(send_lookahead.Buffer(),send_lookahead.BufferSize(),n0);
		}

		evenTime = !evenTime;
	}
};


/// Propagator to advance coherent radiation bidirectionally.
/// Solve d/dt(a) = (0.5i/w0) ((d/dz)^2(a) + w0^2 n^2 a).
export struct IsotropicPropagator:ComputeTool
{
	std::unique_ptr<GlobalIntegrator<tw::Complex>> zGlobalIntegrator;
	tw::vec<tw::Complex> Z1,Z2,Z3;

	IsotropicPropagator(const std::string& name,MetricSpace *m,Task *tsk): ComputeTool(name,m,tsk) {
        zGlobalIntegrator = NULL;

        const tw::Int xDim = space->Dim(1);
        const tw::Int yDim = space->Dim(2);
        const tw::Int zDim = space->Dim(3);

        Z1.resize(zDim);
        Z2.resize(zDim);
        Z3.resize(zDim);

        if (zDim>1) {
            zGlobalIntegrator = std::make_unique<GlobalIntegrator<tw::Complex>>(&task->strip[3],xDim*yDim,zDim);
        }
    }

	void SetupIncomingWaveLeft(const tw::strip& s,ComplexField& amplitude,tw::Complex a0,tw::Complex a1,tw::Complex w0) const
	{
		amplitude.Pack(s, 0, a1 - a0 + ii*w0*space->dl(s,1,3)*half*(a0+a1));
	}

	void SetupIncomingWaveRight(const tw::strip& s,ComplexField& amplitude,tw::Complex aN,tw::Complex aN1,tw::Complex w0) const
	{
		amplitude.Pack(s,space->Dim(s.StripAxis())+1, aN1 - aN - ii*w0*space->dl(s,space->Dim(s.StripAxis())+1,3)*half*(aN+aN1));
	}
    /// Source functions for incoming waves are assumed to be loaded into ghost cells of "amplitude"
	virtual void Advance(ComplexField& amplitude, // this is the array to update
		ComplexField& refractiveIndex, // index of refraction array
		ScalarField& nu_e, // electron collision frequency array
		tw::Float laserFrequency,tw::Float dt)
    {
        const tw::Int xDim = space->Dim(1);
        const tw::Int yDim = space->Dim(2);
        const tw::Int zDim = space->Dim(3);

        tw::Float S1,S2,Vol,D1,D2;
        tw::Complex H1,H2,H3; // related to hamiltonian
        tw::Complex nw,eta_p,eta_m;
        tw::vec<tw::Complex> src,ans;

        const tw::Float dz = space->dl(1,1,1,3); // assume uniform grid

        src.resize(zDim);
        ans.resize(zDim);

        eta_p = one+half*ii*laserFrequency*dz;
        eta_m = one-half*ii*laserFrequency*dz;

        tw::Int idx = 0;
        for (auto j=1;j<=yDim;j++)
            for (auto i=1;i<=xDim;i++)
            {
                for (auto k=1;k<=zDim;k++)
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
                for (auto k=1;k<=zDim;k++)
                    amplitude.Pack(i,j,k,ans[k-1]);

                zGlobalIntegrator->SetMatrix(idx,Z1,Z2,Z3);
                zGlobalIntegrator->SetData(idx,&amplitude(1,i,j,0,0),amplitude.Stride(3),amplitude.Stride(4));
                idx++;
            }

        zGlobalIntegrator->Parallelize();
        amplitude.CopyFromNeighbors();
    //	amplitude.ApplyBoundaryCondition();
        for (auto j=1;j<=yDim;j++)
            for (auto i=1;i<=xDim;i++)
            {
                if (task->n0[3]==MPI_PROC_NULL)
                    amplitude.Pack(i,j,0,amplitude(i,j,1)*eta_p/eta_m - amplitude(i,j,0)/eta_m);
                if (task->n1[3]==MPI_PROC_NULL)
                    amplitude.Pack(i,j,zDim+1,amplitude(i,j,zDim)*eta_p/eta_m + amplitude(i,j,zDim+1)/eta_m);
            }
    }
};
