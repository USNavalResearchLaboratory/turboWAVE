module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_test.h"
#include "tw_logger.h"

export module rel_cold_fluid;
import input;
import driver;
import fct;
import fields;
import diagnostics;
import physics;
import injection;
import chemistry;
import parabolic;
import elliptic;
import logger;

using namespace tw::bc;

export struct Fluid:Driver
{
	tw::Float charge,mass,thermalMomentum,enCrossSection,initialIonizationFraction;
	Field state0,state1; // density,p1,p2,p3 (unlike SPARC state, p is not a momentum density)
	ScalarField fixed,gas;
	bool coulombCollisions;

	// temporaries used in Update
	Field vel; // gammaAvg,v1,v2,v3

	std::shared_ptr<Ionizer> ionizer;
	tw::Profiles profiles;
	tw::Waves waves;

	Field *EM,*J4;
	Field *laser;
	ComplexField *chi;
	tw::Float *carrierFrequency;

	Fluid(const std::string& name,MetricSpace *ms,Task *tsk);
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void VerifyInput();
	virtual void Initialize();
	virtual void Update();
	virtual void MoveWindow();
	virtual void AddDensity(tw::Float densityToAdd,tw::Int i,tw::Int j,tw::Int k);

	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void Report(Diagnostic&);
	virtual void RegisterTests() {
		REGISTER(Fluid,AdvectionTest);
		REGISTER(Fluid,ConservationTest);
	}
	void AdvectionTest();
	void ConservationTest();
};

Fluid::Fluid(const std::string& name,MetricSpace *ms,Task *tsk):Driver(name,ms,tsk)
{
	if (native.unit_system!=tw::units::plasma)
		throw tw::FatalError("Fluid module requires <native units = plasma>");

	charge = -1.0;
	mass = 1.0;
	thermalMomentum = 0.0;
	enCrossSection = 0.0;
	initialIonizationFraction = 1.0;
	coulombCollisions = false;

	state0.Initialize(4,*space,task);
	state1.Initialize(4,*space,task);
	fixed.Initialize(*space,task);
	gas.Initialize(*space,task);

	vel.Initialize(4,*space,task);

	EM = J4 = NULL;
	laser = NULL;
	chi = NULL;
	carrierFrequency = NULL;

	directives.Add("charge",new tw::input::Float(&charge),false);
	directives.Add("mass",new tw::input::Float(&mass),false);
	directives.Add("neutral cross section",new tw::input::Float(&enCrossSection),false);
	directives.Add("initial ionization fraction",new tw::input::Float(&initialIonizationFraction),false);
	directives.Add("coulomb collisions",new tw::input::Bool(&coulombCollisions),false);
}

bool Fluid::InspectResource(void* resource,const std::string& description)
{
	if (description=="electromagnetic:F")
	{
		EM = (Field*)resource;
		return true;
	}

	if (description=="electromagnetic:sources")
	{
		J4 = (Field*)resource;
		return true;
	}

	if (description=="laser:F")
	{
		laser = (Field*)resource;
		return true;
	}

	if (description=="laser:chi")
	{
		chi = (ComplexField*)resource;
		return true;
	}

	if (description=="laser:carrierFrequency")
	{
		carrierFrequency = (tw::Float*)resource;
		return true;
	}

	return false;
}

void Fluid::VerifyInput()
{
	Driver::VerifyInput();
	for (auto tool : tools) {
		if (std::dynamic_pointer_cast<Ionizer>(tool)) {
			ionizer = std::dynamic_pointer_cast<Ionizer>(tool);
		} else if (std::dynamic_pointer_cast<Wave>(tool)) {
			waves.push_back(std::dynamic_pointer_cast<Wave>(tool));
		} else if (std::dynamic_pointer_cast<Profile>(tool)) {
			profiles.push_back(std::dynamic_pointer_cast<Profile>(tool));
		}
	}
}

void Fluid::Initialize()
{
	Driver::Initialize();

	tw::bc::fld xNormal,xParallel,other;
	xNormal = space->bc0[1]==par::axisymmetric || space->bc0[1]==par::reflecting ? fld::dirichletWall : fld::neumannWall;
	xParallel = fld::neumannWall;
	other = fld::neumannWall;

	if (space->IsStdMovingWindow())
	{
		state0.SetBoundaryConditions(All(state0),tw::grid::x,xParallel,other);
		state0.SetBoundaryConditions(All(state0),tw::grid::y,other,other);
		state0.SetBoundaryConditions(All(state0),tw::grid::z,other,fld::none);
		state0.SetBoundaryConditions(Rng(1),tw::grid::x,xNormal,other);

		state1.SetBoundaryConditions(All(state1),tw::grid::x,xParallel,other);
		state1.SetBoundaryConditions(All(state1),tw::grid::y,other,other);
		state1.SetBoundaryConditions(All(state1),tw::grid::z,other,fld::none);
		state1.SetBoundaryConditions(Rng(1),tw::grid::x,xNormal,other);

		gas.SetBoundaryConditions(tw::grid::x,xParallel,other);
		gas.SetBoundaryConditions(tw::grid::y,other,other);
		gas.SetBoundaryConditions(tw::grid::z,other,fld::none);

		vel.SetBoundaryConditions(All(vel),tw::grid::x,xParallel,other);
		vel.SetBoundaryConditions(All(vel),tw::grid::y,other,other);
		vel.SetBoundaryConditions(All(vel),tw::grid::z,other,other);
	}
	else
	{
		state0.SetBoundaryConditions(All(state0),tw::grid::x,fld::neumannWall,fld::neumannWall);
		state0.SetBoundaryConditions(All(state0),tw::grid::y,fld::neumannWall,fld::neumannWall);
		state0.SetBoundaryConditions(All(state0),tw::grid::z,fld::neumannWall,fld::neumannWall);
		state0.SetBoundaryConditions(Rng(1),tw::grid::x,xNormal,fld::neumannWall);
		state0.SetBoundaryConditions(Rng(3),tw::grid::z,fld::dirichletWall,fld::dirichletWall);

		state1.SetBoundaryConditions(All(state1),tw::grid::x,fld::neumannWall,fld::neumannWall);
		state1.SetBoundaryConditions(All(state1),tw::grid::y,fld::neumannWall,fld::neumannWall);
		state1.SetBoundaryConditions(All(state1),tw::grid::z,fld::neumannWall,fld::neumannWall);
		state1.SetBoundaryConditions(Rng(1),tw::grid::x,xNormal,fld::neumannWall);
		state1.SetBoundaryConditions(Rng(3),tw::grid::z,fld::dirichletWall,fld::dirichletWall);

		gas.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
		gas.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
		gas.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);

		vel.SetBoundaryConditions(All(vel),tw::grid::x,fld::neumannWall,fld::neumannWall);
		vel.SetBoundaryConditions(All(vel),tw::grid::y,fld::neumannWall,fld::neumannWall);
		vel.SetBoundaryConditions(All(vel),tw::grid::z,fld::neumannWall,fld::neumannWall);
	}

	#pragma omp parallel
	{
		tw::vec3 A0,A1;
		tw::Float density;
		const tw::Float dth = 0.5*dx(0);

		for (auto cell : EntireCellRange(*this,1))
		{
			auto pos4 = space->Pos4(cell);
			auto pos = pos4.spatial();

			for (auto prof : profiles) {
				for (tw::Int c=1;c<=3;c++) {
					state0(cell,c) = prof->DriftMomentum(1.0)[c-1];
					state1(cell,c) = prof->DriftMomentum(1.0)[c-1];
				}
				density = prof->GetValue(pos4,*space);
				gas(cell) += (1.0 - initialIonizationFraction)*density;
				state0(cell,0) += initialIonizationFraction*density; // ionization fraction should not be zero
				state1(cell,0) += initialIonizationFraction*density;
				if (space->neutralize) {
					fixed(cell) += initialIonizationFraction*density;
				}
			}

			if (carrierFrequency==NULL) {
				A0 = A1 = tw::vec3(0,0,0);
				for (auto w : waves) {
					A0 += w->VectorPotential(-dth,pos);
					A1 += w->VectorPotential(0.0,pos);
				}
				state0(cell,1) += A0.x;
				state0(cell,2) += A0.y;
				state0(cell,3) += 0.5*(A0.x*A0.x + A0.y*A0.y);
				state1(cell,1) += A1.x;
				state1(cell,2) += A1.y;
				state1(cell,3) += 0.5*(A1.x*A1.x + A1.y*A1.y);
			} else {
				for (auto w : waves) {
					state0(cell,3) += 0.25*norm(w->VectorPotentialEnvelope(-dth,pos,*carrierFrequency));
					state1(cell,3) += 0.25*norm(w->VectorPotentialEnvelope(0.0,pos,*carrierFrequency));
				}
			}
		}
	}

	// temperature is set to last profile's temperature
	if (profiles.size()==0)
		throw tw::FatalError("Fluid module needs a profile.");
	thermalMomentum = profiles.back()->thermalMomentum.x;
	if (profiles.back()->temperature!=0.0)
		thermalMomentum = std::sqrt(profiles.back()->temperature*mass); // appropriate for std::exp(-v^2/(2*vth^2)) convention
	if (thermalMomentum==0.0)
		throw tw::FatalError("Fluid module requires temperature specification.");
	if (initialIonizationFraction<=0.0 || initialIonizationFraction>1.0)
		throw tw::FatalError("Fluid module reports initial ionization fraction out of range.");

	state0.CopyFromNeighbors(All(state0));
	state0.ApplyBoundaryCondition(All(state0));
	state1.CopyFromNeighbors(All(state1));
	state1.ApplyBoundaryCondition(All(state1));
	gas.CopyFromNeighbors();
	gas.ApplyBoundaryCondition();
}

void Fluid::MoveWindow()
{
	const tw::Float dth = 0.5*dx(0);
	Driver::MoveWindow();
	#pragma omp parallel
	{
		for (auto s : StripRange(*this,3,0,1,strongbool::yes))
		{
			tw::Int k = Dim(s.Axis())+1;
			tw::vec3 A0,A1;
			tw::Float incomingGas,incomingPlasma[4];
			auto pos4 = space->Pos4(s,k);
			auto pos = pos4.spatial();
			incomingGas = incomingPlasma[0] = incomingPlasma[1] = incomingPlasma[2] = incomingPlasma[3] = 0.0;
			for (auto profile : profiles) {
				if (ionizer==NULL) {
					incomingPlasma[0] += profile->GetValue(pos4,*space);
				} else {
					incomingGas += profile->GetValue(pos4,*space);
					incomingPlasma[0] += 1e-6*profile->GetValue(pos4,*space); // add a little plasma
				}
			}
			state0.Shift(All(state0),s,-1,incomingPlasma);
			state1.Shift(All(state1),s,-1,incomingPlasma);
			fixed.Shift(All(fixed),s,-1,0.0);
			gas.Shift(All(gas),s,-1,incomingGas);
			if (space->neutralize)
				fixed(s,k) += incomingPlasma[0];

			A0 = A1 = tw::vec3(0,0,0);
			for (auto wave : waves) {
				A0 += wave->VectorPotential(space->WindowPos(0)-dth,pos);
				A1 += wave->VectorPotential(space->WindowPos(0),pos);
			}
			state0(s,k,1) += A0.x;
			state0(s,k,2) += A0.y;
			state0(s,k,3) += 0.5*(A0.x*A0.x + A0.y*A0.y);
			state1(s,k,1) += A1.x;
			state1(s,k,2) += A1.y;
			state1(s,k,3) += 0.5*(A1.x*A1.x + A1.y*A1.y);
		}
	}

	state0.DownwardCopy(All(state0),tw::grid::z,1);
	state1.DownwardCopy(All(state1),tw::grid::z,1);
	fixed.DownwardCopy(tw::grid::z,1);
	gas.DownwardCopy(tw::grid::z,1);
}

void Fluid::AddDensity(tw::Float densToAdd,tw::Int i,tw::Int j,tw::Int k)
{
	if (densToAdd>0.0)
	{
		tw::Float scaleFactor = state1(1,i,j,k,0) / (state1(1,i,j,k,0) + densToAdd);
		state1(1,i,j,k,0) += densToAdd;
		state1(1,i,j,k,1) *= scaleFactor;
		state1(1,i,j,k,2) *= scaleFactor;
		state1(1,i,j,k,3) *= scaleFactor;
		if (space->neutralize)
			fixed(i,j,k) += densToAdd;
	}
}

void Fluid::Update()
{
	logger::TRACE("cold fluid update");

	tw::bc::fld bc;
	tw::Float field,ionizedDensity;

	const tw::Float dt = dx(0);
	const tw::Float dti = dk(0);
	const tw::Float dth = 0.5*dx(0);
	const tw::Float q0 = charge;
	const tw::Float m0 = mass;
	const tw::Float m0i = 1.0/m0;
	const tw::Float kT = sqr(thermalMomentum)/m0; // appropriate for std::exp(-v^2/(2*vth^2)) convention
	const tw::Float coulomb = coulombCollisions ? 1.0 : 0.0;

	// Time centering information upon entry:
	// qty		time		qty		time
	// ---		----		---		----
	// n0		-1/2		n1		0
	// 						p1		-1/2
	// a0		-1/2		a1		1/2
	// aa0		0			aa1		1/2
	// E		0			B		1/2

	// Time centering information on output:
	// qty		time		qty					time
	// ---		----		---					----
	// rho		1			chi					1/2
	// J		1/2

	// Advance the momentum to t=1/2

	#pragma omp parallel
	{
		// Apply half of E-field impulse and put relativistic mass into vel(cell,0)
		for (auto v : VectorStripRange<3>(*this,0,1,false))
		{
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
			{
				state1(v,k,1) += dth*q0*(*EM).sfwd(v,k,0,1);
				state1(v,k,2) += dth*q0*(*EM).sfwd(v,k,1,2);
				state1(v,k,3) += dth*q0*(*EM).sfwd(v,k,2,3);
				vel(v,k,0) = m0*m0 + sqr(state1(v,k,1)) + sqr(state1(v,k,2)) + sqr(state1(v,k,3));
			}
			if (laser)
			{
				#pragma omp simd
				for (tw::Int k=1;k<=dim[3];k++)
					vel(v,k,0) += 0.5*q0*q0*(*laser)(v,k,6);
			}
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
				vel(v,k,0) = std::sqrt(vel(v,k,0));
		}
	}

	vel.CopyFromNeighbors(Rng(0));
	vel.ApplyBoundaryCondition(Rng(0));
	#pragma omp parallel
	{
		tw::Float kT_eff,temp;
		std::valarray<tw::Float> nuColl(dim[3]+1);
		const tw::Float nconv = 1.0*tw::dims::density >> native >> cgs;
		const tw::Float Tconv = 1.0*tw::dims::temperature >> native >> cgs;
		const tw::Float fconv = 1.0*tw::dims::frequency >> cgs >> native;
		for (auto v : VectorStripRange<3>(*this,0,1,false))
		{
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
			{
				kT_eff = kT + (vel(v,k,0) - m0);
				nuColl[k] = gas(v,k) * enCrossSection * std::sqrt(kT_eff*m0i);
				temp = coulomb * 1e-5 * fixed(v,k) * nconv;
				temp *= std::pow(kT_eff*Tconv,-1.5);
				nuColl[k] += temp*fconv;
			}
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
			{
				state1(v,k,1) -= dt*vel.d1(v,k,0,1);
				state1(v,k,2) -= dt*vel.d1(v,k,0,2);
				state1(v,k,3) -= dt*vel.d1(v,k,0,3);
				state1(v,k,1) += dth*q0*(*EM).sfwd(v,k,0,1);
				state1(v,k,2) += dth*q0*(*EM).sfwd(v,k,1,2);
				state1(v,k,3) += dth*q0*(*EM).sfwd(v,k,2,3);
				state1(v,k,1) /= 1.0 + nuColl[k]*dt;
				state1(v,k,2) /= 1.0 + nuColl[k]*dt;
				state1(v,k,3) /= 1.0 + nuColl[k]*dt;
			}
		}
	}
	state1.CopyFromNeighbors(All(state1));
	state1.ApplyBoundaryCondition(All(state1));

	// Compute 1/mass and velocity at t = 1/2

	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,0,1,true))
		{
			#pragma omp simd
			for (tw::Int k=lfg[3];k<=ufg[3];k++)
				vel(v,k,0) = m0*m0 + sqr(state1(v,k,1)) + sqr(state1(v,k,2)) + sqr(state1(v,k,3));
			if (laser)
			{
				#pragma omp simd
				for (tw::Int k=lfg[3];k<=ufg[3];k++)
					vel(v,k,0) += 0.5*q0*q0*(*laser)(v,k,7);
			}
			#pragma omp simd
			for (tw::Int k=lfg[3];k<=ufg[3];k++)
			{
				vel(v,k,0) = 1.0/std::sqrt(vel(v,k,0));
				vel(v,k,1) = state1(v,k,1)*vel(v,k,0);
				vel(v,k,2) = state1(v,k,2)*vel(v,k,0);
				vel(v,k,3) = state1(v,k,3)*vel(v,k,0);
			}
		}
	}

	// FCT Update of density - trial step

	state0 = state1; // dens0(t=0) , dens1(t=0)
	FCT_Driver convector(&state0,&state1,&vel,NULL,space);
	convector.SetDensityElements(Rng(0));
	if (dim[1]>1)
	{
		convector.SetVelocityElement(1);
		convector.Convect(tw::grid::x, tw::bc::fld::dirichletCell, tw::bc::fld::dirichletCell, dth);
	}
	if (dim[2]>1)
	{
		convector.SetVelocityElement(2);
		convector.Convect(tw::grid::y, tw::bc::fld::dirichletCell, tw::bc::fld::dirichletCell, dth);
	}
	if (dim[3]>1)
	{
		bc = space->IsStdMovingWindow() ? tw::bc::fld::neumannWall : tw::bc::fld::dirichletCell;
		convector.SetVelocityElement(3);
		convector.Convect(tw::grid::z,bc,bc,dth);
	}

	// FCT Update of density - full step

	Swap(state0,state1); // dens0(t=0) , dens1(t=1/2)
	if (dim[1]>1)
	{
		convector.SetVelocityElement(1);
		convector.Convect(tw::grid::x, tw::bc::fld::dirichletCell, tw::bc::fld::dirichletCell, dt);
		convector.GetTrueFlux(vel,Rng(1),Rng(0));
	}
	if (dim[2]>1)
	{
		convector.SetVelocityElement(2);
		convector.Convect(tw::grid::y, tw::bc::fld::dirichletCell, tw::bc::fld::dirichletCell, dt);
		convector.GetTrueFlux(vel,Rng(2),Rng(0));
	}
	if (dim[3]>1)
	{
		bc = space->IsStdMovingWindow() ? tw::bc::fld::neumannWall : tw::bc::fld::dirichletCell;
		convector.SetVelocityElement(3);
		convector.Convect(tw::grid::z,bc,bc,dt);
		convector.GetTrueFlux(vel,Rng(3),Rng(0));
	}
	Swap(state0,state1); // dens0(t=1/2) , dens1(t=1)

	// Ionization

	if (ionizer!=NULL)
	{
		for (tw::Int i=lfg[1];i<=ufg[1];i++)
			for (tw::Int j=lfg[2];j<=ufg[2];j++)
				for (tw::Int k=lfg[3];k<=ufg[3];k++)
				{
					field = std::sqrt(sqr((*EM)(1,i,j,k,0))+sqr((*EM)(1,i,j,k,1))+sqr((*EM)(1,i,j,k,2)));
					ionizedDensity = gas(i,j,k)*dt*ionizer->InstantRate(1e-6,field);
					if (ionizedDensity > gas(i,j,k))
						ionizedDensity = gas(i,j,k);
					gas(i,j,k) = std::fabs(gas(i,j,k) - ionizedDensity);
					AddDensity(ionizedDensity,i,j,k);
				}
	}

	// Deposit current at t = 1/2 and charge at t = 1
	// Get conservative current based on actual FCT fluxes which are stored
	// in vel(cell,i) ; these are currents integrated over a cell wall.

	if (J4)
	{
		// This is a deposition routine.
		// Writing to ghost cells would double count the sources.
		#pragma omp parallel
		{
			// Need special treatment for current in ignorable directions.
			// Here cell walls and centers need not be distinguished (dim=1).
			if (dim[1]==1)
			{
				for (auto v : VectorStripRange<3>(*this,0,1,false))
					#pragma omp simd
					for (tw::Int k=1;k<=dim[3];k++)
						vel(v,k,1) *= space->dS(v,k,1) * dt * state0(v,k,0);
			}
			if (dim[2]==1)
			{
				for (auto v : VectorStripRange<3>(*this,0,1,false))
					#pragma omp simd
					for (tw::Int k=1;k<=dim[3];k++)
						vel(v,k,2) *= space->dS(v,k,2) * dt * state0(v,k,0);
			}
			if (dim[3]==1)
			{
				for (auto v : VectorStripRange<3>(*this,0,1,false))
					#pragma omp simd
					for (tw::Int k=1;k<=dim[3];k++)
						vel(v,k,3) *= space->dS(v,k,3) * dt * state0(v,k,0);
			}
		}
		#pragma omp parallel
		{
			for (auto v : VectorStripRange<3>(*this,0,1,false))
			{
				#pragma omp simd
				for (tw::Int k=1;k<=dim[3];k++)
				{
					(*J4)(v,k,0) += space->dS(v,k,0) * q0 * (state1(v,k,0) - fixed(v,k));
					(*J4)(v,k,1) += q0 * vel(v,k,1) * dti;
					(*J4)(v,k,2) += q0 * vel(v,k,2) * dti;
					(*J4)(v,k,3) += q0 * vel(v,k,3) * dti;
				}
			}
		}
	}

	if (chi)
	{
		#pragma omp parallel
		{
			for (auto cell : InteriorCellRange(*this,1)) {
				(*chi)(cell,0) -= space->dS(cell,0)*q0*q0*state0(cell,0)/(m0*vel(cell,0));
			}
		}
	}
}

void Fluid::ReadCheckpoint(std::ifstream& inFile)
{
	Driver::ReadCheckpoint(inFile);
	state0.ReadCheckpoint(inFile);
	state1.ReadCheckpoint(inFile);
	gas.ReadCheckpoint(inFile);
	fixed.ReadCheckpoint(inFile);
}

void Fluid::WriteCheckpoint(std::ofstream& outFile)
{
	Driver::WriteCheckpoint(outFile);
	state0.WriteCheckpoint(outFile);
	state1.WriteCheckpoint(outFile);
	gas.WriteCheckpoint(outFile);
	fixed.WriteCheckpoint(outFile);
}

void Fluid::Report(Diagnostic& diagnostic)
{
	Driver::Report(diagnostic);
	diagnostic.ReportField(name+"_e",state1,1,0,tw::dims::density,"$n_e$");
	diagnostic.ReportField(name+"_n",gas,1,0,tw::dims::density,"$n_g$");
}

