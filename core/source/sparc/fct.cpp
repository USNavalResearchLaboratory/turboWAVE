module;

#include "tw_includes.h"

export module fct;
import fields;

export struct FCT_Engine
{
	// cells is the number of grid cells, not counting the 2 ghost cells
	// The V and A arrays have to be set up to contain cell volumes and low-side cell wall areas
	tw::Int cells;
	std::valarray<tw::Float> V, A, scratch;

	FCT_Engine(tw::Int ax, const MetricSpace& m);
	void Reset(const tw::strip& s, const MetricSpace& m, ScalarField* fluxMask);
	void Transport(std::valarray<tw::Float>& vel, std::valarray<tw::Float>& rho, std::valarray<tw::Float>& rho1, std::valarray<tw::Float>& diff, std::valarray<tw::Float>& flux, tw::Float dt);
	void Diffuse(std::valarray<tw::Float>& vel, std::valarray<tw::Float>& rho, std::valarray<tw::Float>& diff, tw::Float dt);
	void Limiter(tw::Float& adiff, const tw::Float& maxLow, const tw::Float& maxHigh);
	void Clip(std::valarray<tw::Float>& rho, std::valarray<tw::Float>& adiff, tw::Float rho00);
	void AntiDiffuse(std::valarray<tw::Float>& rho, std::valarray<tw::Float>& adiff, std::valarray<tw::Float>& flux);
};

export struct FCT_Driver
{
	Element en; // indices of density components
	tw::Int vi; // index to velocity component
	Field* diff, * net_flux;

	// externally owned data
	Field* rho;
	Field* rho1;
	Field* vel;
	ScalarField* fluxMask;
	MetricSpace* ms;

	FCT_Driver(Field* rho, Field* rho1, Field* vel, ScalarField* fluxMask, MetricSpace* ms);
	~FCT_Driver();
	void SetDensityElements(const Element& e) { en = e; }
	void SetVelocityElement(tw::Int v) { vi = v; }
	void Convect(const tw::grid::axis& axis, tw::bc::fld low, tw::bc::fld high, tw::Float dt);
	void GetTrueFlux(Field& flux, const Element& dst, const Element& src)
	{
		CopyFieldData(flux, dst, *net_flux, src);
	}
};


////////////////////
//                //
//   FCT ENGINE   //
//                //
////////////////////


FCT_Engine::FCT_Engine(tw::Int ax,const MetricSpace& m)
{
	cells = m.Dim(ax);
	V.resize(cells+2);
	A.resize(cells+2);
	scratch.resize(cells+2);
}

void FCT_Engine::Reset(const tw::strip& s,const MetricSpace& m,ScalarField *fluxMask)
{
	tw::Int i,ax=s.Axis();
	cells = m.Dim(ax);
	// The engine uses V[0] in the clipping stage
	// The engine never uses A[0]
	for (i=0;i<=cells+1;i++)
		V[i] = m.dS(s,i,0);
	for (i=1;i<=cells+1;i++)
		A[i] = m.dS(s,i,ax);
	if (fluxMask!=NULL)
		for (i=1;i<=cells+1;i++)
			A[i] *= 1.0 - tw::Float( (*fluxMask)(s,i-1) + (*fluxMask)(s,i) == 1.0 );
}

void FCT_Engine::Transport(std::valarray<tw::Float>& vel,
	std::valarray<tw::Float>& rho,
	std::valarray<tw::Float>& rho1,
	std::valarray<tw::Float>& diff,
	std::valarray<tw::Float>& flux,
	tw::Float dt)
{
	// First step in FCT algorithm
	// vel, rho, rho1, flux, and diff are indexed from 0..cells+1
	// vel contains velocity in the direction of transport
	// vel is known at t = 0 for a 1rst order push, at t = 1/2 for a second order push
	// diff and flux are given at the low-side cell walls, all other quantities at the cell centers
	// rho is density at t = 0
	// rho1 is density at t = 0 for a 1rst order push (rho1 = rho), or at t = 1/2 for a 2nd order push
	// diff is an output giving the diffused "mass" into the low side of each cell
	// flux is an output giving the transported and diffused "mass" into the low side of each cell
	// after calling this, the values of the low and high ghost cells for rho must be supplied

	tw::Int i;

	// compute diffused mass
	#pragma omp simd
	for (i=1;i<=cells+1;i++)
		diff[i] = -0.25*A[i]*dt*std::fabs(vel[i]+vel[i-1])*(rho[i] - rho[i-1]);

	// compute transported mass
	#pragma omp simd
	for (i=1;i<=cells+1;i++)
		flux[i] = A[i]*dt*0.25*(vel[i]+vel[i-1])*(rho1[i]+rho1[i-1]);

	// compute transported density
	#pragma omp simd
	for (i=1;i<=cells;i++)
		rho[i] += (flux[i] - flux[i+1])/V[i];

	// compute transported and diffused mass
	#pragma omp simd
	for (i=1;i<=cells+1;i++)
		flux[i] += diff[i];
}

void FCT_Engine::Diffuse(std::valarray<tw::Float>& vel,
	std::valarray<tw::Float>& rho,
	std::valarray<tw::Float>& diff,
	tw::Float dt)
{
	// vel is the same as for FCT_Engine::Transport
	// rho is the output from FCT_Engine::Transport, with appropriate ghost cell values
	// diff is the array output from FCT_Engine::Transport, with appropriate ghost cell values
	// on output, diff is replaced by the anti-diffused mass into the low-side of each cell
	// after calling this, the values of the low/high ghost cells for rho must be supplied
	// HOWEVER, the GLOBAL ghost cells for rho should not be updated (leave them with the undiffused values)

	tw::Int i;
	scratch = rho;

	// compute transported and diffused density
	#pragma omp simd
	for (i=1;i<=cells;i++)
		rho[i] += (diff[i] - diff[i+1])/V[i];

	// compute anti-diffused mass (re-using diffused mass array)
	#pragma omp simd
	for (i=1;i<=cells+1;i++)
		diff[i] = 0.25*A[i]*dt*std::fabs(vel[i]+vel[i-1])*(scratch[i] - scratch[i-1]);
}

void FCT_Engine::Limiter(tw::Float& adiff,const tw::Float& maxLow,const tw::Float& maxHigh)
{
	tw::Float pos_channel=adiff,neg_channel=adiff,test;

	tw::Float maxHigh_pos = tw::Float(maxHigh>0.0);
	tw::Float maxLow_pos = tw::Float(maxLow>0.0);

	test = tw::Float(pos_channel>maxHigh);
	pos_channel = maxHigh_pos*(test*maxHigh + (1.0-test)*pos_channel);
	test = tw::Float(pos_channel>maxLow);
	pos_channel = maxLow_pos*(test*maxLow + (1.0-test)*pos_channel);

	test = tw::Float(neg_channel<maxHigh);
	neg_channel = (1.0-maxHigh_pos)*(test*maxHigh + (1.0-test)*neg_channel);
	test = tw::Float(neg_channel<maxLow);
	neg_channel = (1.0-maxLow_pos)*(test*maxLow + (1.0-test)*neg_channel);

	adiff = pos_channel*tw::Float(adiff>0.0) + neg_channel*tw::Float(adiff<=0.0);
}

void FCT_Engine::Clip(std::valarray<tw::Float>& rho,std::valarray<tw::Float>& adiff,tw::Float rho00)
{
	// Limit the antidiffusive fluxes given in adiff to ensure stability, positivity, etc.
	// rho00 gives the value in the cell below the low-side ghost cell
	// rho and adiff are the outputs from FCT_Engine::Diffuse
	// after calling, high ghost cell for adiff must be supplied

	tw::Int i;
	tw::Float maxAntiDiffusionLow,maxAntiDiffusionHigh;

	maxAntiDiffusionHigh = (rho[2]-rho[1])*V[1];
	maxAntiDiffusionLow = (rho[0] - rho00)*V[0];
	Limiter(adiff[1],maxAntiDiffusionLow,maxAntiDiffusionHigh);
	#pragma omp simd
	for (i=2;i<=cells;i++)
	{
		maxAntiDiffusionHigh = (rho[i+1]-rho[i])*V[i];
		maxAntiDiffusionLow = (rho[i-1]-rho[i-2])*V[i-1];
		Limiter(adiff[i],maxAntiDiffusionLow,maxAntiDiffusionHigh);
	}
}

void FCT_Engine::AntiDiffuse(std::valarray<tw::Float>& rho,std::valarray<tw::Float>& adiff,std::valarray<tw::Float>& flux)
{
	// adiff are the clipped fluxes from FCT_Engine::Clip
	// after calling, low and high ghost cells for rho must be supplied
	// flux is an ouptut giving the true flux through the cell walls
	// after calling, the low ghost cell for flux must be supplied

	tw::Int i;
	const tw::Float safety_factor = 0.999999;
	#pragma omp simd
	for (i=1;i<=cells;i++)
		rho[i] += safety_factor*(adiff[i] - adiff[i+1])/V[i];
	#pragma omp simd
	for (i=1;i<=cells+1;i++)
		flux[i] += safety_factor*adiff[i];
}



////////////////////
//                //
//   FCT DRIVER   //
//                //
////////////////////


FCT_Driver::FCT_Driver(Field *rho,Field *rho1,Field *vel,ScalarField *fluxMask,MetricSpace *ms)
{
	en = All(*rho);
	vi = 0;
	this->rho = rho;
	this->rho1 = rho1;
	this->vel = vel;
	this->fluxMask = fluxMask;
	this->ms = ms;
	diff = NULL;
	net_flux = NULL;
}

FCT_Driver::~FCT_Driver()
{
	if (diff!=NULL)
		delete diff;
	if (net_flux!=NULL)
		delete net_flux;
}

void FCT_Driver::Convect(const tw::grid::axis& axis,tw::bc::fld low,tw::bc::fld high,tw::Float dt)
{
	const tw::Int ax = tw::grid::naxis(axis);
	const tw::Int N = ms->Dim(ax);

	if (diff!=NULL)
		delete diff;
	if (net_flux!=NULL)
		delete net_flux;
	diff = new Field;
	net_flux = new Field;
	diff->Initialize(en.Components(),*rho,rho->task);
	net_flux->Initialize(en.Components(),*rho,rho->task);
	diff->SetBoundaryConditions(axis,low,high);
	net_flux->SetBoundaryConditions(axis,low,high);

	#pragma omp parallel
	{
		tw::Int c;
		std::valarray<tw::Float> va_rho(N+2),va_rho1(N+2),va_vel(N+2),va_diff(N+2),va_flux(N+2);
		tw::Float va_rho00;
		FCT_Engine engine(ax,*ms);

		// TRANSPORT

		for (c=en.low;c<=en.high;c++)
		{
			for (auto s : StripRange(*ms,ax,strongbool::yes))
			{
				engine.Reset(s,*ms,fluxMask);
				rho->GetStrip(va_rho,s,c);
				rho1->GetStrip(va_rho1,s,c);
				vel->GetStrip(va_vel,s,vi);
				engine.Transport(va_vel,va_rho,va_rho1,va_diff,va_flux,dt);
				rho->SetStrip(va_rho,s,c);
				diff->SetStrip(va_diff,s,c-en.low);
				net_flux->SetStrip(va_flux,s,c-en.low);
			}
		}
		#pragma omp barrier
		#pragma omp single
		{
			rho->DownwardCopy(axis,en,1);
			rho->UpwardCopy(axis,en,1);
			rho->ApplyBoundaryCondition(en);
		}

		// DIFFUSE

		for (c=en.low;c<=en.high;c++)
		{
			for (auto s : StripRange(*ms,ax,strongbool::yes))
			{
				engine.Reset(s,*ms,fluxMask);
				rho->GetStrip(va_rho,s,c);
				diff->GetStrip(va_diff,s,c-en.low);
				vel->GetStrip(va_vel,s,vi);
				engine.Diffuse(va_vel,va_rho,va_diff,dt);
				rho->SetStrip(va_rho,s,c);
				diff->SetStrip(va_diff,s,c-en.low);
			}
		}
		#pragma omp barrier
		#pragma omp single
		{
			rho->DownwardCopy(axis,en,1);
			rho->UpwardCopy(axis,en,2); // need 2 low-side ghost cells for clipping operation
			// leave ghost cell with transported un-diffused value (no call to rho->ApplyBoundaryCondition)
		}

		// CLIP

		for (c=en.low;c<=en.high;c++)
		{
			for (auto s : StripRange(*ms,ax,strongbool::yes))
			{
				engine.Reset(s,*ms,fluxMask);
				rho->GetStrip(va_rho,s,c);
				diff->GetStrip(va_diff,s,c-en.low);
				va_rho00 = (*rho)(s,-1,c);
				engine.Clip(va_rho,va_diff,va_rho00);
				diff->SetStrip(va_diff,s,c-en.low);
			}
		}
		#pragma omp barrier
		#pragma omp single
		{
			diff->DownwardCopy(axis,1);
			diff->ApplyBoundaryCondition();
		}

		// ANTI-DIFFUSE

		for (c=en.low;c<=en.high;c++)
		{
			for (auto s : StripRange(*ms,ax,strongbool::yes))
			{
				engine.Reset(s,*ms,fluxMask);
				rho->GetStrip(va_rho,s,c);
				diff->GetStrip(va_diff,s,c-en.low);
				net_flux->GetStrip(va_flux,s,c-en.low);
				engine.AntiDiffuse(va_rho,va_diff,va_flux);
				rho->SetStrip(va_rho,s,c);
				net_flux->SetStrip(va_flux,s,c-en.low);
			}
		}
		#pragma omp barrier
		#pragma omp single
		{
			rho->DownwardCopy(axis,en,1);
			rho->UpwardCopy(axis,en,1);
			rho->ApplyBoundaryCondition(en);
			net_flux->UpwardCopy(axis,1);
			net_flux->ApplyBoundaryCondition();
		}
	}
}
