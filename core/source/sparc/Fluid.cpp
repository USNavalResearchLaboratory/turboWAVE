module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module fluid;
import input;
import compute_tool;
import twmodule;
import fct;
import fields;
import diagnostics;
import physics;
import injection;
import chemistry;
import parabolic;
import elliptic;
import logger;
#include "tw_logger.h"

using namespace tw::bc;

namespace sparc
{
	enum laserModel { vacuum, isotropic };
	enum radiationModel { noRadiation, thin, thick };
	enum plasmaModel { neutral, quasineutral };
}

export struct Fluid:Module
{
	tw::Float charge,mass,thermalMomentum,enCrossSection,initialIonizationFraction;
	Field state0,state1; // density,p1,p2,p3 (unlike SPARC state, p is not a momentum density)
	ScalarField fixed,gas;
	bool coulombCollisions;

	// temporaries used in Update
	Field vel; // gammaAvg,v1,v2,v3

	Ionizer *ionizer;

	Field *EM,*J4;
	Field *laser;
	ComplexField *chi;
	tw::Float *carrierFrequency;

	Fluid(const std::string& name,Simulation* sim);
	~Fluid();
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void VerifyInput();
	virtual void Initialize();
	virtual void Update();
	virtual void MoveWindow();
	virtual void AddDensity(tw::Float densityToAdd,tw::Int i,tw::Int j,tw::Int k);

	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void Report(Diagnostic&);
	virtual bool Test(tw::Int& id);
	void AdvectionTest();
	void ConservationTest();
};

export struct EquilibriumGroup;

export struct Chemical:Module
{
	EquilibriumGroup *group; // explictly typed super
	EOSComponent *eosData;
	Ionizer *ionizer;
	UniformProfile *background;
	sparc::material mat;
	tw::Int indexInState;

	Chemical(const std::string& name,Simulation* sim);
	virtual ~Chemical();
	void SetupIndexing();
	virtual void VerifyInput();

	bool LoadFluid(Field& hydro);
	void LoadInternalEnergy(Field& hydro,Field& eos);
};

export struct EquilibriumGroup:Module
{
	std::vector<Chemical*> chemical; // explicitly typed submodule list
	EOSMixture *eosMixData;
	bool mobile;

	// The hydro set contains indices into the state vector for this group.
	// The mass density index corresponds to the first chemical in the group.
	sparc::hydro_set hidx;
	// The eos set contains indices into the eos vector for this group.
	sparc::eos_set eidx;
	// The material data is packed and processed for optimization
	sparc::material_set matset;
	// Following is used to limit motion by zeroing forces on this group
	tw::Float forceFilter;

	tw::Float DensitySum(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=hidx.first;s<hidx.first+hidx.num;s++)
			ans += f(cell,s);
		return ans;
	}
	tw::Float DensityWeightedSum(const Field& f,std::valarray<tw::Float>& qty,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=hidx.first;s<hidx.first+hidx.num;s++)
			ans += f(cell,s)*qty[s-hidx.first];
		return ans;
	}
	void LoadMassDensity(ScalarField& nm,const Field& f)
	{
		for (auto cell : EntireCellRange(*this))
			nm(cell) = DensityWeightedSum(f,matset.mass,cell);
	}
	void LoadMassDensityCv(ScalarField& nmcv,const Field& f)
	{
		for (auto cell : EntireCellRange(*this))
			nmcv(cell) = DensityWeightedSum(f,matset.cvm,cell);
	}
	tw::vec3 Velocity(const Field& f,const tw::cell& cell)
	{
		tw::Float nm = DensityWeightedSum(f,matset.mass,cell);
		return tw::vec3(f(cell,hidx.npx),f(cell,hidx.npy),f(cell,hidx.npz))/(tw::small_pos + nm);
	}
	void LoadVelocity(ScalarField& vel,const Field& f,tw::Int ax)
	{
		// assumes velocity components appear in order in state vector
		tw::Float nm;
		for (auto cell : EntireCellRange(*this))
		{
			nm = DensityWeightedSum(f,matset.mass,cell);
			vel(cell) = f(cell,hidx.npx+ax-1)/(tw::small_pos + nm);
		}
	}

	EquilibriumGroup(const std::string& name,Simulation* sim);
	virtual ~EquilibriumGroup(); // ASHER_MOD
	void SetupIndexing();
	virtual void VerifyInput();
	bool GenerateFluid(Field& hydro,Field& eos);
};

export namespace sparc
{
struct HydroManager:Module
{
	std::vector<EquilibriumGroup*> group; // explicitly type submodule list
	std::vector<Reaction*> reaction;
	std::vector<Excitation*> excitation;
	std::vector<Collision*> collision;
	Field state0,state1,creationRate,destructionRate,eos0,eos1;

	ParabolicSolver *parabolicSolver;
	EllipticSolver *ellipticSolver;
	IsotropicPropagator *laserPropagator;
	tw::vec3 dipoleCenter;
	tw::Float laserFrequency; // derived from pulse list during initialization
	tw::Float backgroundDensity,backgroundTemperature;

	sparc::radiationModel radModel;
	sparc::laserModel lasModel;
	sparc::plasmaModel plasModel;
	bool electrostaticHeating;
	Chemical *electrons; // pointer to electron species, if present
	tw::Int ie; // index to electron density in state vector
	ScalarField scratch,scratch2,fluxMask;
	ScalarField rho0,rho,phi,nu_e,me_eff,radiativeLosses,radiationIntensity;
	ComplexField laserAmplitude,refractiveIndex;

	// items needed for step size control and user feedback
	tw::Float epsilonFactor;
	std::stringstream statusMessage;

	// The distinction between creation and destruction is used for adaptive timestep adjustment.
	// In the case of constant time step, there is no difference between create(x) and destroy(-x).
	// For signed quantities like momentum, the distinction is never used at all.
	void CreateMass(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ creationRate(cell,h.ni) += val; }
	void DestroyMass(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ destructionRate(cell,h.ni) += val; }
	void CreateMomentum(const tw::cell& cell,const tw::Int& ax,const tw::Float& val,const sparc::hydro_set& h)
		{ creationRate(cell,h.npx+ax-1) += val; }
	void DestroyMomentum(const tw::cell& cell,const tw::Int& ax,const tw::Float& val,const sparc::hydro_set& h)
		{ destructionRate(cell,h.npx+ax-1) += val; }
	void CreateTotalEnergy(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ creationRate(cell,h.u) += val; }
	void DestroyTotalEnergy(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ destructionRate(cell,h.u) += val; }
	void CreateVibrations(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ creationRate(cell,h.x) += val; }
	void DestroyVibrations(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
		{ destructionRate(cell,h.x) += val; }
	void CreateTotalAndVibrational(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
	{
		creationRate(cell,h.u) += val;
		creationRate(cell,h.x) += val;
	}
	void DestroyTotalAndVibrational(const tw::cell& cell,const tw::Float& val,const sparc::hydro_set& h)
	{
		destructionRate(cell,h.u) += val;
		destructionRate(cell,h.x) += val;
	}

	HydroManager(const std::string& name,Simulation* sim);
	virtual ~HydroManager();
	void SetupIndexing();
	virtual void Initialize();
	virtual void Reset();
	void LoadCollisionRate(Collision *coll,ScalarField& R);
	void ComputeElectronCollisionFrequency();
	void ComputeCollisionalSources();
	void ComputeRadiativeSources();
	void ComputeHydroSources();
	tw::vec3 ComputeForceOnBody(tw::Int i,tw::Int j,tw::Int k);
	void ComputeSources();
	tw::Float EstimateTimeStep();

	void HydroAdvance(const tw::grid::axis& axis,tw::Float dt);
	void LaserAdvance(tw::Float dt);
	void ChemAdvance(tw::Float dt);
	void DiffusionAdvance(tw::Float dt);
	void FieldAdvance(tw::Float dt);
	void EOSAdvance(tw::Float dt);
	void FirstOrderAdvance(tw::Float dt,bool computeSources);
	virtual void Update();

	virtual void VerifyInput();
	virtual bool ReadQuasitoolBlock(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);

	virtual void Report(Diagnostic&);
	virtual void StatusMessage(std::ostream *dest);
};
}

/////////////////////////////
//                         //
// RELATIVISTIC COLD FLUID //
//                         //
/////////////////////////////

Fluid::Fluid(const std::string& name,Simulation* sim):Module(name,sim)
{
	if (native.native!=tw::units::plasma)
		throw tw::FatalError("Fluid module requires <native units = plasma>");

	charge = -1.0;
	mass = 1.0;
	thermalMomentum = 0.0;
	enCrossSection = 0.0;
	initialIonizationFraction = 1.0;
	coulombCollisions = false;

	state0.Initialize(4,*this,owner);
	state1.Initialize(4,*this,owner);
	fixed.Initialize(*this,owner);
	gas.Initialize(*this,owner);

	vel.Initialize(4,*this,owner);

	EM = J4 = NULL;
	laser = NULL;
	chi = NULL;
	carrierFrequency = NULL;
	ionizer = NULL;

	directives.Add("charge",new tw::input::Float(&charge),false);
	directives.Add("mass",new tw::input::Float(&mass),false);
	directives.Add("neutral cross section",new tw::input::Float(&enCrossSection),false);
	directives.Add("initial ionization fraction",new tw::input::Float(&initialIonizationFraction),false);
	directives.Add("coulomb collisions",new tw::input::Bool(&coulombCollisions),false);
}

Fluid::~Fluid()
{
	if (ionizer!=NULL)
		owner->RemoveTool(ionizer);
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
	Module::VerifyInput();
	for (auto tool : moduleTool)
	{
		ionizer = dynamic_cast<Ionizer*>(tool);
		if (ionizer!=NULL)
			break;
	}
}

void Fluid::Initialize()
{
	Module::Initialize();

	tw::bc::fld xNormal,xParallel,other;
	xNormal = owner->bc0[1]==par::axisymmetric || owner->bc0[1]==par::reflecting ? fld::dirichletWall : fld::neumannWall;
	xParallel = fld::neumannWall;
	other = fld::neumannWall;

	if (owner->movingWindow)
	{
		state0.SetBoundaryConditions(tw::grid::x,xParallel,other);
		state0.SetBoundaryConditions(tw::grid::y,other,other);
		state0.SetBoundaryConditions(tw::grid::z,other,fld::none);
		state0.SetBoundaryConditions(Element(1),tw::grid::x,xNormal,other);

		state1.SetBoundaryConditions(tw::grid::x,xParallel,other);
		state1.SetBoundaryConditions(tw::grid::y,other,other);
		state1.SetBoundaryConditions(tw::grid::z,other,fld::none);
		state1.SetBoundaryConditions(Element(1),tw::grid::x,xNormal,other);

		gas.SetBoundaryConditions(tw::grid::x,xParallel,other);
		gas.SetBoundaryConditions(tw::grid::y,other,other);
		gas.SetBoundaryConditions(tw::grid::z,other,fld::none);

		vel.SetBoundaryConditions(tw::grid::x,xParallel,other);
		vel.SetBoundaryConditions(tw::grid::y,other,other);
		vel.SetBoundaryConditions(tw::grid::z,other,other);
	}
	else
	{
		state0.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
		state0.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
		state0.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
		state0.SetBoundaryConditions(Element(1),tw::grid::x,xNormal,fld::neumannWall);
		state0.SetBoundaryConditions(Element(3),tw::grid::z,fld::dirichletWall,fld::dirichletWall);

		state1.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
		state1.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
		state1.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
		state1.SetBoundaryConditions(Element(1),tw::grid::x,xNormal,fld::neumannWall);
		state1.SetBoundaryConditions(Element(3),tw::grid::z,fld::dirichletWall,fld::dirichletWall);

		gas.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
		gas.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
		gas.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);

		vel.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
		vel.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
		vel.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	}

	#pragma omp parallel
	{
		tw::vec3 pos,A0,A1;
		tw::Float density;
		const tw::Float dth = 0.5*dx(0);

		for (auto cell : EntireCellRange(*this))
		{
			pos = owner->Pos(cell);

			for (auto prof : profile)
			{
				for (tw::Int c=1;c<=3;c++)
				{
					state0(cell,c) = prof->DriftMomentum(1.0)[c-1];
					state1(cell,c) = prof->DriftMomentum(1.0)[c-1];
				}
				density = prof->GetValue(pos,*owner);
				gas(cell) += (1.0 - initialIonizationFraction)*density;
				state0(cell,0) += initialIonizationFraction*density; // ionization fraction should not be zero
				state1(cell,0) += initialIonizationFraction*density;
				if (owner->neutralize)
					fixed(cell) += initialIonizationFraction*density;
			}

			if (carrierFrequency==NULL)
			{
				A0 = A1 = tw::vec3(0,0,0);
				for (auto w : wave)
				{
					A0 += w->VectorPotential(-dth,pos);
					A1 += w->VectorPotential(0.0,pos);
				}
				state0(cell,1) += A0.x;
				state0(cell,2) += A0.y;
				state0(cell,3) += 0.5*(A0.x*A0.x + A0.y*A0.y);
				state1(cell,1) += A1.x;
				state1(cell,2) += A1.y;
				state1(cell,3) += 0.5*(A1.x*A1.x + A1.y*A1.y);
			}
			else
			{
				for (auto w : wave)
				{
					state0(cell,3) += 0.25*norm(w->VectorPotentialEnvelope(-dth,pos,*carrierFrequency));
					state1(cell,3) += 0.25*norm(w->VectorPotentialEnvelope(0.0,pos,*carrierFrequency));
				}
			}
		}
	}

	// temperature is set to last profile's temperature
	if (profile.size()==0)
		throw tw::FatalError("Fluid module needs a profile.");
	thermalMomentum = profile.back()->thermalMomentum.x;
	if (profile.back()->temperature!=0.0)
		thermalMomentum = std::sqrt(profile.back()->temperature*mass); // appropriate for std::exp(-v^2/(2*vth^2)) convention
	if (thermalMomentum==0.0)
		throw tw::FatalError("Fluid module requires temperature specification.");
	if (initialIonizationFraction<=0.0 || initialIonizationFraction>1.0)
		throw tw::FatalError("Fluid module reports initial ionization fraction out of range.");

	state0.CopyFromNeighbors();
	state0.ApplyBoundaryCondition();
	state1.CopyFromNeighbors();
	state1.ApplyBoundaryCondition();
	gas.CopyFromNeighbors();
	gas.ApplyBoundaryCondition();
}

void Fluid::MoveWindow()
{
	const tw::Float dth = 0.5*dx(0);
	Module::MoveWindow();
	#pragma omp parallel
	{
		for (auto s : StripRange(*this,3,strongbool::yes))
		{
			tw::Int k = Dim(s.Axis())+1;
			tw::vec3 pos,A0,A1;
			tw::Float incomingGas,incomingPlasma[4];
			pos = owner->Pos(s,k);
			incomingGas = incomingPlasma[0] = incomingPlasma[1] = incomingPlasma[2] = incomingPlasma[3] = 0.0;
			for (tw::Int p=0;p<profile.size();p++)
			{
				if (ionizer==NULL)
					incomingPlasma[0] += profile[p]->GetValue(pos,*owner);
				else
				{
					incomingGas += profile[p]->GetValue(pos,*owner);
					incomingPlasma[0] += 1e-6*profile[p]->GetValue(pos,*owner); // add a little plasma
				}
			}
			state0.Shift(s,-1,incomingPlasma);
			state1.Shift(s,-1,incomingPlasma);
			fixed.Shift(s,-1,0.0);
			gas.Shift(s,-1,incomingGas);
			if (owner->neutralize)
				fixed(s,k) += incomingPlasma[0];

			A0 = A1 = tw::vec3(0,0,0);
			for (tw::Int w=0;w<wave.size();w++)
			{
				A0 += wave[w]->VectorPotential(owner->WindowPos(0)-dth,pos);
				A1 += wave[w]->VectorPotential(owner->WindowPos(0),pos);
			}
			state0(s,k,1) += A0.x;
			state0(s,k,2) += A0.y;
			state0(s,k,3) += 0.5*(A0.x*A0.x + A0.y*A0.y);
			state1(s,k,1) += A1.x;
			state1(s,k,2) += A1.y;
			state1(s,k,3) += 0.5*(A1.x*A1.x + A1.y*A1.y);
		}
	}

	state0.DownwardCopy(tw::grid::z,1);
	state1.DownwardCopy(tw::grid::z,1);
	fixed.DownwardCopy(tw::grid::z,1);
	gas.DownwardCopy(tw::grid::z,1);
}

void Fluid::AddDensity(tw::Float densToAdd,tw::Int i,tw::Int j,tw::Int k)
{
	if (densToAdd>0.0)
	{
		tw::Float scaleFactor = state1(i,j,k,0) / (state1(i,j,k,0) + densToAdd);
		state1(i,j,k,0) += densToAdd;
		state1(i,j,k,1) *= scaleFactor;
		state1(i,j,k,2) *= scaleFactor;
		state1(i,j,k,3) *= scaleFactor;
		if (owner->neutralize)
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
		for (auto v : VectorStripRange<3>(*this,false))
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

	vel.CopyFromNeighbors(Element(0));
	vel.ApplyBoundaryCondition(Element(0));
	#pragma omp parallel
	{
		tw::Float kT_eff,temp;
		std::valarray<tw::Float> nuColl(dim[3]+1);
		const tw::Float nconv = 1.0*tw::dims::density >> native >> cgs;
		const tw::Float Tconv = 1.0*tw::dims::temperature >> native >> cgs;
		const tw::Float fconv = 1.0*tw::dims::frequency >> cgs >> native;
		for (auto v : VectorStripRange<3>(*this,false))
		{
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
			{
				kT_eff = kT + (vel(v,k,0) - m0);
				nuColl[k] = gas(v,k) * enCrossSection * std::sqrt(kT_eff*m0i);
				temp = coulomb * 1e-5 * fixed(v,k) * nconv;
				temp *= pow(kT_eff*Tconv,-1.5);
				nuColl[k] += temp*fconv;
			}
			#pragma omp simd
			for (tw::Int k=1;k<=dim[3];k++)
			{
				state1(v,k,1) -= dt*vel(v,k,0,1);
				state1(v,k,2) -= dt*vel(v,k,0,2);
				state1(v,k,3) -= dt*vel(v,k,0,3);
				state1(v,k,1) += dth*q0*(*EM).sfwd(v,k,0,1);
				state1(v,k,2) += dth*q0*(*EM).sfwd(v,k,1,2);
				state1(v,k,3) += dth*q0*(*EM).sfwd(v,k,2,3);
				state1(v,k,1) /= 1.0 + nuColl[k]*dt;
				state1(v,k,2) /= 1.0 + nuColl[k]*dt;
				state1(v,k,3) /= 1.0 + nuColl[k]*dt;
			}
		}
	}
	state1.CopyFromNeighbors();
	state1.ApplyBoundaryCondition();

	// Compute 1/mass and velocity at t = 1/2

	#pragma omp parallel
	{
		for (auto v : VectorStripRange<3>(*this,true))
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
	FCT_Driver convector(&state0,&state1,&vel,NULL,owner);
	convector.SetDensityElements(Element(0));
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
		bc = owner->movingWindow ? tw::bc::fld::neumannWall : tw::bc::fld::dirichletCell;
		convector.SetVelocityElement(3);
		convector.Convect(tw::grid::z,bc,bc,dth);
	}

	// FCT Update of density - full step

	Swap(state0,state1); // dens0(t=0) , dens1(t=1/2)
	if (dim[1]>1)
	{
		convector.SetVelocityElement(1);
		convector.Convect(tw::grid::x, tw::bc::fld::dirichletCell, tw::bc::fld::dirichletCell, dt);
		convector.GetTrueFlux(vel,Element(1),Element(0));
	}
	if (dim[2]>1)
	{
		convector.SetVelocityElement(2);
		convector.Convect(tw::grid::y, tw::bc::fld::dirichletCell, tw::bc::fld::dirichletCell, dt);
		convector.GetTrueFlux(vel,Element(2),Element(0));
	}
	if (dim[3]>1)
	{
		bc = owner->movingWindow ? tw::bc::fld::neumannWall : tw::bc::fld::dirichletCell;
		convector.SetVelocityElement(3);
		convector.Convect(tw::grid::z,bc,bc,dt);
		convector.GetTrueFlux(vel,Element(3),Element(0));
	}
	Swap(state0,state1); // dens0(t=1/2) , dens1(t=1)

	// Ionization

	if (ionizer!=NULL)
	{
		for (tw::Int i=lfg[1];i<=ufg[1];i++)
			for (tw::Int j=lfg[2];j<=ufg[2];j++)
				for (tw::Int k=lfg[3];k<=ufg[3];k++)
				{
					field = std::sqrt(sqr((*EM)(i,j,k,0))+sqr((*EM)(i,j,k,1))+sqr((*EM)(i,j,k,2)));
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
				for (auto v : VectorStripRange<3>(*this,false))
					#pragma omp simd
					for (tw::Int k=1;k<=dim[3];k++)
						vel(v,k,1) *= owner->dS(v,k,1) * dt * state0(v,k,0);
			}
			if (dim[2]==1)
			{
				for (auto v : VectorStripRange<3>(*this,false))
					#pragma omp simd
					for (tw::Int k=1;k<=dim[3];k++)
						vel(v,k,2) *= owner->dS(v,k,2) * dt * state0(v,k,0);
			}
			if (dim[3]==1)
			{
				for (auto v : VectorStripRange<3>(*this,false))
					#pragma omp simd
					for (tw::Int k=1;k<=dim[3];k++)
						vel(v,k,3) *= owner->dS(v,k,3) * dt * state0(v,k,0);
			}
		}
		#pragma omp parallel
		{
			for (auto v : VectorStripRange<3>(*this,false))
			{
				#pragma omp simd
				for (tw::Int k=1;k<=dim[3];k++)
				{
					(*J4)(v,k,0) += owner->dS(v,k,0) * q0 * (state1(v,k,0) - fixed(v,k));
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
			for (auto cell : InteriorCellRange(*this))
				(*chi)(cell) -= owner->dS(cell,0)*q0*q0*state0(cell,0)/(m0*vel(cell,0));
		}
	}
}

void Fluid::ReadCheckpoint(std::ifstream& inFile)
{
	Module::ReadCheckpoint(inFile);
	state0.ReadCheckpoint(inFile);
	state1.ReadCheckpoint(inFile);
	gas.ReadCheckpoint(inFile);
	fixed.ReadCheckpoint(inFile);
}

void Fluid::WriteCheckpoint(std::ofstream& outFile)
{
	Module::WriteCheckpoint(outFile);
	state0.WriteCheckpoint(outFile);
	state1.WriteCheckpoint(outFile);
	gas.WriteCheckpoint(outFile);
	fixed.WriteCheckpoint(outFile);
}

void Fluid::Report(Diagnostic& diagnostic)
{
	Module::Report(diagnostic);
	diagnostic.ReportField(name+"_e",state1,0,tw::dims::density,"$n_e$");
	diagnostic.ReportField(name+"_n",gas,0,tw::dims::density,"$n_g$");
}

/////////////////////
//                 //
// CHEMICAL MODULE //
//                 //
/////////////////////


Chemical::Chemical(const std::string& name,Simulation* sim):Module(name,sim)
{
	if (native.native!=tw::units::plasma)
		throw tw::FatalError("Chemical module requires <native units = plasma>");

	mat.charge = -1.0;
	mat.mass = 1.0;
	mat.cvm = 1.5;
	mat.excitationEnergy = 0.0;
	mat.thermometricConductivity = 0.0;
	mat.kinematicViscosity = 0.0;
	mat.eps[0] = 1.0;
	mat.eps[1] = 0.0;
	eosData = NULL;
	ionizer = NULL;
	background = NULL;
	mat.AddDirectives(directives);
}

// ASHER_MOD -- Chemical now needs a destructor because of the eosData
Chemical::~Chemical()
{
	if (eosData!=NULL)
		owner->RemoveTool(eosData);
	if (ionizer!=NULL)
		owner->RemoveTool(ionizer);
}

void Chemical::VerifyInput()
{
	Module::VerifyInput();
	group = (EquilibriumGroup*)super;
	// DFG - Find an EOS on the list of tools associated with this module.
	// This list is populated automatically by the base Module class.
	for (auto tool : moduleTool)
	{
		eosData = dynamic_cast<EOSComponent*>(tool);
		if (eosData!=NULL)
			break;
	}
	for (auto tool : moduleTool)
	{
		ionizer = dynamic_cast<Ionizer*>(tool);
		if (ionizer!=NULL)
			break;
	}
	// If the EOS tool could not be found, create one automatically.
	// Another approach would be to throw an error, if we want to force the user to be explicit (this would break old input files).
	// Ionization tool is optional, so if one was not found do nothing.
	if (eosData==NULL)
	{
		if (mat.mass==1.0)
			eosData = (EOSComponent*)owner->CreateTool("default_hot_electrons",tw::tool_type::eosHotElectrons);
		else
			eosData = (EOSComponent*)owner->CreateTool("default_ideal_gas",tw::tool_type::eosIdealGas);
	}
	// Add a uniform profile for the automatic background fluid.
	sparc::HydroManager *hydro = dynamic_cast<sparc::HydroManager*>(super->super);
	if (hydro->backgroundDensity > 0.0)
	{
		background = (UniformProfile*)owner->CreateTool("auto_background",tw::tool_type::uniformProfile);
		profile.push_back(background);
	}
}

void Chemical::SetupIndexing()
{
	// Called by EquilibriumGroup::SetupIndexing()
	// DFG - the ionization object is now a ComputeTool
	if (ionizer!=NULL)
	{
		// Setup the indexing for photoionization here (ionization tool cannot do it)
		// Since Hydromanager starts this chain, we can count on all EquilibriumGroup indexing to be in place.
		Chemical *echem = (Chemical*)owner->GetModule(ionizer->electron_name);
		Chemical *ichem = (Chemical*)owner->GetModule(ionizer->ion_name);
		ionizer->hgas = group->hidx;
		ionizer->hgas.ni = indexInState;
		ionizer->he = ((EquilibriumGroup*)echem->super)->hidx;
		ionizer->he.ni = echem->indexInState;
		ionizer->hi = ((EquilibriumGroup*)ichem->super)->hidx;
		ionizer->hi.ni = ichem->indexInState;
	}
	// Have to send indexing data to EOS
	eosData->SetupIndexing(indexInState,group->hidx,group->eidx,mat);
	if (mat.mass==1.0)
	{
		sparc::HydroManager *master = (sparc::HydroManager*)(super->super);
		master->ie = indexInState;
		master->electrons = this;
		group->forceFilter = 0.0;
	}
}

bool Chemical::LoadFluid(Field& hydro)
{
	// Load the mass, momentum, kinetic energy, and explicitly loaded internal energy.
	// Implicit internal energy is handled in a subsequent sweep.
	// Explicitly tracked vibrational energy is added here.

	bool massLoaded = false;
	tw::Float add = 0.0;

	// DFG - indices into the state vector are encapsulated in hidx
	// the type of hidx is sparc::hydro_set defined in physics.h
	// (similarly, EOS state indices are in eidx)
	const tw::Int ns = indexInState;
	const tw::Int npx = group->hidx.npx;
	const tw::Int npy = group->hidx.npy;
	const tw::Int npz = group->hidx.npz;
	const tw::Int U = group->hidx.u;
	const tw::Int Xi = group->hidx.x;

	for (auto prof : profile)
	{
		if ( prof->TimeGate(owner->WindowPos(0),&add) )
		{
			const tw::vec3 p0 = prof->DriftMomentum(mat.mass);
			for (auto cell : EntireCellRange(*this))
			{
				const tw::Float dens = prof->GetValue(owner->Pos(cell),*owner);
				if (prof->whichQuantity==tw::profile::quantity::density && dens>0.0)
				{
					massLoaded = true;
					const tw::Float kT = prof->Temperature(mat.mass);
					const tw::Float kinetic = 0.5*Norm(dens*p0)/(tw::small_pos + mat.mass*dens);
					const tw::Float vibrational = dens*mat.excitationEnergy/(std::fabs(std::exp(mat.excitationEnergy/kT) - 1.0) + tw::small_pos);
					hydro(cell,ns) = add*hydro(cell,ns) + dens;
					hydro(cell,npx) = add*hydro(cell,npx) + dens*p0.x;
					hydro(cell,npy) = add*hydro(cell,npy) + dens*p0.y;
					hydro(cell,npz) = add*hydro(cell,npz) + dens*p0.z;
					hydro(cell,U) = add*hydro(cell,U) + kinetic + vibrational; // internal energy added in subsequent sweep
					hydro(cell,Xi) = add*hydro(cell,Xi) + vibrational;
				}
				if (prof->whichQuantity==tw::profile::quantity::energy)
				{
					hydro(cell,U) = add*hydro(cell,U) + dens;
				}
				if (prof->whichQuantity==tw::profile::quantity::px)
				{
					hydro(cell,npx) = add*hydro(cell,npx) + dens;
				}
				if (prof->whichQuantity==tw::profile::quantity::py)
				{
					hydro(cell,npy) = add*hydro(cell,npy) + dens;
				}
				if (prof->whichQuantity==tw::profile::quantity::pz)
				{
					hydro(cell,npz) = add*hydro(cell,npz) + dens;
				}
			}
		}
	}
	hydro.ApplyBoundaryCondition(Element(ns,Xi));
	return massLoaded;
}

void Chemical::LoadInternalEnergy(Field& hydro,Field& eos)
{
	// Load internal energy, attempting to respect profile targets.
	// In this process the eos array is overwritten with temporary data.
	// Therefore eos has to be rebuilt after calling this routine.

	sparc::HydroManager *master = (sparc::HydroManager*)(super->super);

	for (auto prof : profile)
	{
		if (prof->whichQuantity==tw::profile::quantity::density)
		{
			// Get the target temperature for this profile
			const tw::Float kT = prof->Temperature(mat.mass);
			for (auto cell : EntireCellRange(*this))
			{
				// Put the target mass density in the scratch array
				master->scratch(cell) = mat.mass * prof->GetValue(owner->Pos(cell),*owner);
				// Put the target temperature into the eos array
				eos(cell,group->eidx.T) = kT;
				// Use scratch2 to store the reference temperature, currently hard coded to zero
				master->scratch2(cell) = 0.0;
			}
			// Set heat capacity based on the target mass density and the target temperature
			eosData->SetHeatCapacity(master->scratch,eos);
			// Add internal energy based on target mass density, reference temperature, heat capacity, and target temperature
			// The hydro array is already loaded with the totals for the mass and momentum densities.
			// The tool can access everything via its own indexing data.
			group->eosMixData->UpdateEnergy(master->scratch,master->scratch2,hydro,eos);
		}
	}
	hydro.ApplyBoundaryCondition(Element(group->hidx.u));
}


/////////////////////
//                 //
//   GROUP MODULE  //
//                 //
/////////////////////


EquilibriumGroup::EquilibriumGroup(const std::string& name,Simulation* sim):Module(name,sim)
{
	if (native.native!=tw::units::plasma)
		throw tw::FatalError("EquilibriumGroup module requires <native units = plasma>");

	mobile = true;
	forceFilter = 1.0;

	// DFG - Start with NULL tool.
	// User can select one by name, or let it be created automatically in VerifyInput()
	eosMixData = NULL;

	directives.Add("mobile",new tw::input::Bool(&mobile));
}

// ASHER_MOD -- EquilibriumGroup now needs a destructor because of eosMixData
EquilibriumGroup::~EquilibriumGroup()
{
	if (eosMixData!=NULL)
		owner->RemoveTool(eosMixData);
}

void EquilibriumGroup::SetupIndexing()
{
	// Called by HydroManager::SetupIndexing()
	matset.Allocate(chemical.size());
	for (tw::Int i=0;i<chemical.size();i++)
		matset.AddMaterial(chemical[i]->mat,i);
	eosMixData->SetupIndexing(hidx,eidx,matset);
	// First need to set indexInState for all chemicals
	for (tw::Int i=0;i<chemical.size();i++)
		chemical[i]->indexInState = hidx.first + i;
	// Now we can setup indexing for all Chemical modules
	for (auto chem : chemical)
		chem->SetupIndexing();
}

void EquilibriumGroup::VerifyInput()
{
	Module::VerifyInput();
	// Extract Chemical modules from the submodule list
	for (auto sub : submodule)
	{
		Chemical *chem = dynamic_cast<Chemical*>(sub);
		if (chem!=NULL)
			chemical.push_back(chem);
	}
	// DFG - Find an EOS mixture on the list of tools associated with this module.
	// This list is populated automatically by the base Module class.
	for (auto tool : moduleTool)
	{
		eosMixData = dynamic_cast<EOSMixture*>(tool);
		if (eosMixData!=NULL)
			break;
	}
	// If the tool could not be found, create one automatically.
	// Another approach would be to throw an error, if we want to force the user to be explicit (this would break old input files).
	if (eosMixData==NULL)
		eosMixData = (EOSMixture *)owner->CreateTool("default_eos_mix",tw::tool_type::eosMixture);
}

bool EquilibriumGroup::GenerateFluid(Field& hydro,Field& eos)
{
	bool didGenerate = false;
	std::vector<bool> massLoaded(chemical.size());
	for (tw::Int i=0;i<chemical.size();i++)
		massLoaded[i] = chemical[i]->LoadFluid(hydro);
	for (tw::Int i=0;i<chemical.size();i++)
		if (massLoaded[i])
			chemical[i]->LoadInternalEnergy(hydro,eos);
	for (auto loaded : massLoaded)
		didGenerate |= loaded;
	return didGenerate;
}


///////////////////////////////////////
//                                   //
// SPARC Hydrodynamics Master Module //
//                                   //
///////////////////////////////////////


sparc::HydroManager::HydroManager(const std::string& name,Simulation* sim):Module(name,sim)
{
	if (native.native!=tw::units::plasma)
		throw tw::FatalError("HydroManager module requires <native units = plasma>");

	epsilonFactor = 1e-4;
	laserFrequency = 1.0;
	backgroundDensity = backgroundTemperature = 0.0;

	// some arrays must be initialized later since we don't know how many elements they have yet
	scratch.Initialize(*this,owner);
	scratch2.Initialize(*this,owner);
	fluxMask.Initialize(*this,owner);
	rho0.Initialize(*this,owner);
	rho.Initialize(*this,owner);
	phi.Initialize(*this,owner);
	nu_e.Initialize(*this,owner);
	me_eff.Initialize(*this,owner);
	laserAmplitude.Initialize(*this,owner);
	radiativeLosses.Initialize(*this,owner);
	radiationIntensity.Initialize(*this,owner);
	refractiveIndex.Initialize(*this,owner);

	// DFG - start all tools with null pointers
	parabolicSolver = NULL;
	ellipticSolver = NULL;
	laserPropagator = NULL;

	radModel = sparc::noRadiation;
	lasModel = sparc::vacuum;
	plasModel = sparc::neutral;
	electrons = NULL;
	electrostaticHeating = false;

	directives.Add("epsilon factor",new tw::input::Float(&epsilonFactor),false);
	std::map<std::string,sparc::radiationModel> rad = {{"none",sparc::noRadiation},{"thin",sparc::thin},{"thick",sparc::thick}};
	directives.Add("radiation model",new tw::input::Enums<sparc::radiationModel>(rad,&radModel),false);
	std::map<std::string,sparc::laserModel> las = {{"vacuum",sparc::vacuum},{"isotropic",sparc::isotropic}};
	directives.Add("laser model",new tw::input::Enums<sparc::laserModel>(las,&lasModel),false);
	std::map<std::string,sparc::plasmaModel> plas = {{"neutral",sparc::neutral},{"quasineutral",sparc::quasineutral}};
	directives.Add("plasma model",new tw::input::Enums<sparc::plasmaModel>(plas,&plasModel),false);
	directives.Add("dipole center",new tw::input::Vec3(&dipoleCenter),false);
	directives.Add("background density",new tw::input::Float(&backgroundDensity),false);
	directives.Add("background temperature",new tw::input::Float(&backgroundTemperature),false);
	directives.Add("electrostatic heating",new tw::input::Bool(&electrostaticHeating),false);
}

sparc::HydroManager::~HydroManager()
{
	if (parabolicSolver!=NULL)
		owner->RemoveTool(parabolicSolver);
	if (ellipticSolver!=NULL)
		owner->RemoveTool(ellipticSolver);
	if (laserPropagator!=NULL)
		owner->RemoveTool(laserPropagator);
	for (auto quasitool : reaction)
		delete quasitool;
	for (auto quasitool : excitation)
		delete quasitool;
	for (auto quasitool : collision)
		delete quasitool;
}

void sparc::HydroManager::SetupIndexing()
{
	tw::Int r;

	// EOS indexing
	// DFG - note how setting up the indexing has simplified.
	r = 0; // running index
	for (auto grp : group)
		r = grp->eidx.Load(r);

	// Hydro indexing
	r = 0; // running index
	for (auto grp : group)
		r = grp->hidx.Load(r,grp->chemical.size());

	// Pass indexing down to lower level objects.
	for (auto grp : group)
		grp->SetupIndexing();

	// Under new system, collisions and reactions also have to be indexed.
	// First define some lambdas to help, then index.

	auto GetHydroSet = [&] (const std::string& module_name)
	{
		sparc::hydro_set ans;
		Chemical *chem = (Chemical*)owner->GetModule(module_name);
		ans = chem->group->hidx;
		ans.ni = chem->indexInState;
		return ans;
	};

	auto GetEOSSet = [&] (const std::string& module_name)
	{
		Chemical *chem = (Chemical*)owner->GetModule(module_name);
		return chem->group->eidx;
	};

	auto GetMaterial = [&] (const std::string& module_name)
	{
		Chemical *chem = (Chemical*)owner->GetModule(module_name);
		return chem->mat;
	};

	for (auto rxn : reaction)
	{
		rxn->catalyst = GetEOSSet(rxn->catalyst_name);
		for (auto sub : rxn->sub)
		{
			for (auto str : sub->reactant_names)
			{
				sub->reactants.push_back(GetHydroSet(str));
				sub->mat_r.push_back(GetMaterial(str));
			}
			for (auto str : sub->product_names)
			{
				sub->products.push_back(GetHydroSet(str));
				sub->mat_p.push_back(GetMaterial(str));
			}
		}
	}

	for (auto coll : collision)
	{
		coll->h1 = GetHydroSet(coll->name1);
		coll->e1 = GetEOSSet(coll->name1);
		coll->h2 = GetHydroSet(coll->name2);
		coll->e2 = GetEOSSet(coll->name2);
		coll->m1 = GetMaterial(coll->name1);
		coll->m2 = GetMaterial(coll->name2);
	}

	for (auto x : excitation)
	{
		x->h1 = GetHydroSet(x->name1);
		x->e1 = GetEOSSet(x->name1);
		x->h2 = GetHydroSet(x->name2);
		x->e2 = GetEOSSet(x->name2);
		x->m1 = GetMaterial(x->name1);
		x->m2 = GetMaterial(x->name2);
	}
}

void sparc::HydroManager::VerifyInput()
{
	Module::VerifyInput();

	if (backgroundDensity!=0.0)
	{
		if (backgroundDensity<=0.0)
			throw tw::FatalError("Hydro module background density must be positive.");
		if (backgroundTemperature<=0.0)
			throw tw::FatalError("Hydro module background temperature must be positive.");
	}

	// Search submodule list for EquilibriumGroup modules

	for (auto sub : submodule)
	{
		EquilibriumGroup *grp = dynamic_cast<EquilibriumGroup*>(sub);
		if (grp!=NULL)
			group.push_back(grp);
	}

	// DFG - here is an example of where we have more than 1 kind of tool to consider.
	// Simple appoach is just repeat the structure for each one.

	for (auto tool : moduleTool)
	{
		ellipticSolver = dynamic_cast<EllipticSolver*>(tool);
		if (ellipticSolver!=NULL)
			break;
	}
	for (auto tool : moduleTool)
	{
		parabolicSolver = dynamic_cast<ParabolicSolver*>(tool);
		if (parabolicSolver!=NULL)
			break;
	}
	for (auto tool : moduleTool)
	{
		laserPropagator = dynamic_cast<IsotropicPropagator*>(tool);
		if (laserPropagator!=NULL)
			break;
	}

	if (parabolicSolver==NULL)
		parabolicSolver = (ParabolicSolver*)owner->CreateTool("default_parabolic_solver",tw::tool_type::generalParabolicPropagator);

	if (ellipticSolver==NULL)
	{
		if (owner->SpatialDims()==1)
			ellipticSolver = (EllipticSolver*)owner->CreateTool("default_elliptic_solver",tw::tool_type::ellipticSolver1D);
		else
			ellipticSolver = (EllipticSolver*)owner->CreateTool("default_elliptic_solver",tw::tool_type::iterativePoissonSolver);
		// Default tool needs some default boundaries; if a tool was attached the user will have specified them.
		tw::bc::fld x0 = owner->bc0[1]==par::axisymmetric ? fld::neumannWall : fld::dirichletCell;
		tw::bc::fld z0 = fld::dirichletCell;
		tw::bc::fld z1 = fld::neumannWall;
		ellipticSolver->SetBoundaryConditions(x0,z0,z0,z0,z0,z1);
	}

	if (laserPropagator==NULL)
		laserPropagator = (IsotropicPropagator*)owner->CreateTool("default_laser_propagator",tw::tool_type::isotropicPropagator);
}

void sparc::HydroManager::Initialize()
{
	Module::Initialize();

	// Initialize laser frequency and refractive index
	if (wave.size())
	{
		laserFrequency = 0.0;
		for (auto pulse : wave)
			laserFrequency += pulse->w;
		laserFrequency /= tw::Float(wave.size());
	}
	refractiveIndex = tw::Complex(1.0,0.0);

	ellipticSolver->SetFieldsBoundaryConditions(phi,Element(0));
	rho.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	rho.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	rho.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	scratch.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	scratch.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	scratch.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);

	// Find number of state variables characterizing the whole system
	// DFG - now we have a static variable to help keep count of these more reliably
	// For more see sparc::hydro_set and sparc::eos_set in physics.h
	tw::Int hydroElements = 0;
	for (auto grp : group)
		hydroElements += grp->chemical.size() + sparc::hydro_set::count;
	tw::Int eosElements = group.size()*sparc::eos_set::count;
	state0.Initialize(hydroElements,*this,owner);
	state1.Initialize(hydroElements,*this,owner);
	creationRate.Initialize(hydroElements,*this,owner);
	destructionRate.Initialize(hydroElements,*this,owner);
	eos0.Initialize(eosElements,*this,owner);
	eos1.Initialize(eosElements,*this,owner);

	// variables defining normal component boundary conditions
	tw::bc::fld bc0[4],bc1[4];
	for (tw::Int ax=1;ax<=3;ax++)
	{
		bc0[ax] = (owner->bc0[ax] == par::reflecting || owner->bc0[ax] == par::axisymmetric) ? fld::dirichletWall : fld::neumannWall;
		bc1[ax] = (owner->bc1[ax] == par::reflecting || owner->bc1[ax] == par::axisymmetric) ? fld::dirichletWall : fld::neumannWall;
	}

	// setup the flux mask used to make conductors impermeable
	// also used for reflecting boundary conditions at simulation walls
	fluxMask = 1.0;
	for (auto cell : EntireCellRange(*this))
		for (auto c : conductor)
			if (c->theRgn->Inside(owner->Pos(cell),*owner))
				fluxMask(cell) = 0.0;
	for (tw::Int ax=1;ax<=3;ax++)
		for (auto strip : StripRange(*this,ax,strongbool::yes))
		{
			if (bc0[ax]==fld::dirichletWall && owner->n0[ax]==MPI_PROC_NULL && dim[ax]>1)
				fluxMask(strip,lfg[ax]) = fluxMask(strip,lng[ax]) = 0.0;
			if (bc1[ax]==fld::dirichletWall && owner->n1[ax]==MPI_PROC_NULL && dim[ax]>1)
				fluxMask(strip,ufg[ax]) = fluxMask(strip,ung[ax]) = 0.0;
		}

	// Default boundary conditions, refine after setting up indexing.
	state0.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	state0.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	state0.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	state1.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	state1.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	state1.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	eos0.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	eos0.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	eos0.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	eos1.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	eos1.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	eos1.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	creationRate.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	creationRate.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	creationRate.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	destructionRate.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	destructionRate.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	destructionRate.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	nu_e.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	nu_e.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	nu_e.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);

	// Prepare for loading fluid: indexing, background profiles, and refined boundary conditions

	// This function indexes the entire hydro tree
	SetupIndexing();

	// Setup the automatic background profiles
	if (backgroundDensity!=0.0)
	{
		tw::Int Nminus=0,Nplus=0;
		for (auto grp : group)
			for (auto chem : grp->chemical)
			{
				if (chem->mat.charge < tw::small_neg) Nminus++;
				if (chem->mat.charge > tw::small_pos) Nplus++;
			}
		for (auto grp : group)
			for (auto chem : grp->chemical)
			{
				chem->background->density = backgroundDensity;
				chem->background->temperature = backgroundTemperature;
				const tw::Float Q = chem->mat.charge;
				if (Q < tw::small_neg) chem->background->density /= -Q*Nminus;
				if (Q > tw::small_pos) chem->background->density /= Q*Nplus;
			}
	}

	// Refine boundary conditions for particular quantities
	for (auto grp : group)
	{
		Element e;
		// Temperature
		e = Element(grp->eidx.T);
		parabolicSolver->SetFieldsBoundaryConditions(eos0,e);
		parabolicSolver->SetFieldsBoundaryConditions(eos1,e);
		// X-Component
		e = Element(grp->hidx.npx);
		state0.SetBoundaryConditions(e,tw::grid::x,bc0[1],bc1[1]);
		state1.SetBoundaryConditions(e,tw::grid::x,bc0[1],bc1[1]);
		creationRate.SetBoundaryConditions(e,tw::grid::x,bc0[1],bc1[1]);
		destructionRate.SetBoundaryConditions(e,tw::grid::x,bc0[1],bc1[1]);
		// Y-Component
		e = Element(grp->hidx.npy);
		state0.SetBoundaryConditions(e,tw::grid::y,bc0[2],bc1[2]);
		state1.SetBoundaryConditions(e,tw::grid::y,bc0[2],bc1[2]);
		creationRate.SetBoundaryConditions(e,tw::grid::y,bc0[2],bc1[2]);
		destructionRate.SetBoundaryConditions(e,tw::grid::y,bc0[2],bc1[2]);
		// Z-Component
		e = Element(grp->hidx.npz);
		state0.SetBoundaryConditions(e,tw::grid::z,bc0[3],bc1[3]);
		state1.SetBoundaryConditions(e,tw::grid::z,bc0[3],bc1[3]);
		creationRate.SetBoundaryConditions(e,tw::grid::z,bc0[3],bc1[3]);
		destructionRate.SetBoundaryConditions(e,tw::grid::z,bc0[3],bc1[3]);
	}

	// Add fluid from all profiles
	for (auto grp : group)
	{
		grp->forceFilter = grp->mobile ? 1.0 : 0.0;
		grp->GenerateFluid(state1,eos1);
	}

	me_eff = 1.0; // effective mass reserved for future use
	nu_e = 1.0; // put arbitrary value so initial transport coefficients don't blow up

	// Initialize EOS quantities
	EOSAdvance(0.0); // gets eos0 using state1 only
	eos1 = eos0;
	ComputeElectronCollisionFrequency();
}

void sparc::HydroManager::Reset()
{
	state0 = state1;
	eos0 = eos1;
	rho0 = rho;

	if (!owner->IsFirstStep())
	{
		// Generate new fluid due to sources
		// N.b. EquilibriumGroup::GenerateFluid is destructive to the eos field that is passed in.
		bool didGenerate = false;
		for (auto grp : group)
			didGenerate |= grp->GenerateFluid(state1,eos1);
		if (didGenerate)
		{
			EOSAdvance(dx(0)); // gets eos0 using state0 and state1
			eos1 = eos0;
			state0 = state1;
		}
	}
}

void sparc::HydroManager::LoadCollisionRate(Collision *coll,ScalarField& R)
{
	// Load collision rate into R (units of volume/time).  Can be called from OpenMP parallel section.
	// Deriving a particular frequency depends on the process.
	// Energy and momentum transfer have a factor of 3 difference, and different mass dependence.
	// If this coefficient is multiplied by the field density, we get what is typically called "collision frequency".

	// We will work in CGS units for these computations

	const tw::Float m1 = coll->m1.mass * tw::dims::mass >> native >> cgs;
	const tw::Float m2 = coll->m2.mass * tw::dims::mass >> native >> cgs;
	const tw::Float m12 = m1*m2/(m1+m2);
	const tw::Float q1 = coll->m1.charge * tw::dims::charge >> native >> cgs;
	const tw::Float q2 = coll->m2.charge * tw::dims::charge >> native >> cgs;
	const tw::Float sigma = coll->crossSection * tw::dims::cross_section >> native >> cgs;
	const tw::Float EFermi = coll->T_ref * tw::dims::energy >> native >> cgs;
	const tw::Float nref = coll->n_ref * tw::dims::density >> native >> cgs;
	// Get multiplicative factors so we don't have dimensional numbers in the inner loops
	const tw::Float nconv = 1.0 * tw::dims::density >> native >> cgs;
	const tw::Float econv = 1.0 * tw::dims::energy >> native >> cgs;
	const tw::Float rconv = 1.0 * tw::dims::rate_coefficient_2 >> cgs >> native;
	// CGS quantities that depend on the cell
	tw::Float N1,N2,T1,T2,v12,Ti,phonon,coulomb;

	if (coll->type==sparc::hard_sphere)
		for (auto cell : InteriorCellRange(*this))
		{
			T1 = eos1(cell,coll->e1.T) * econv;
			T2 = eos1(cell,coll->e2.T) * econv;
			v12 = std::sqrt(8.0*(T1/m1 + T2/m2)/pi);
			R(cell) = rconv * (4.0/3.0) * v12 * sigma;
		}

	if (coll->type==sparc::coulomb)
		for (auto cell : InteriorCellRange(*this))
		{
			N1 = state1(cell,coll->h1.ni) * nconv;
			N2 = state1(cell,coll->h2.ni) * nconv;
			T1 = eos1(cell,coll->e1.T) * econv;
			T2 = eos1(cell,coll->e2.T) * econv;
			v12 = std::sqrt(8.0*(T1/m1 + T2/m2)/pi);
			R(cell) = rconv * (4.0/3.0) * v12 * sparc::CoulombCrossSectionCGS(q1,q2,m12,v12,N1,N2,T1,T2);
		}

	if (coll->type==sparc::metallic)
		// Electron-phonon rate is the collision frequency divided by a reference density
		// The collision frequency is defined by the momentum loss rate, e.g., dp/dt = -nu*p
		for (auto cell : InteriorCellRange(*this))
		{
			N1 = state1(cell,coll->h1.ni) * nconv;
			N2 = state1(cell,coll->h2.ni) * nconv;
			T1 = eos1(cell,coll->e1.T) * econv;
			T2 = eos1(cell,coll->e2.T) * econv;
			v12 = std::sqrt(8.0*(T1/m1 + T2/m2)/pi);
			Ti = m1==1.0 ? T2 : T1;
			phonon = sparc::ElectronPhononFrequencyCGS(Ti,EFermi,coll->ks)/nref;
			coulomb = (4.0/3.0) * v12 * sparc::CoulombCrossSectionCGS(q1,q2,m12,v12,N1,N2,T1,T2);
			// This expression causes the rate to go over to Spitzer when phonon rate >> coulomb rate.
			// This is also where temperature is large, as required.
			R(cell) = rconv*coulomb*phonon / (coulomb + phonon);
		}
}

void sparc::HydroManager::ComputeElectronCollisionFrequency()
{
	#pragma omp parallel
	{
		for (auto cell : InteriorCellRange(*this))
			nu_e(cell) = 0.0;
		for (auto coll : collision)
		{
			LoadCollisionRate(coll,scratch);
			if (coll->h1.ni==ie)
				for (auto cell : InteriorCellRange(*this))
					nu_e(cell) += state1(cell,coll->h2.ni) * scratch(cell);
			if (coll->h2.ni==ie)
				for (auto cell : InteriorCellRange(*this))
					nu_e(cell) += state1(cell,coll->h1.ni) * scratch(cell);
		}
	}
	nu_e.CopyFromNeighbors();
	nu_e.ApplyBoundaryCondition();
}

void sparc::HydroManager::ComputeCollisionalSources()
{
	#pragma omp parallel
	{
		tw::Float rateNow;

		// REACTIONS

		for (auto rx : reaction)
			for (auto cell : InteriorCellRange(*this))
			{
				rateNow = rx->PrimitiveRate(eos1(cell,rx->catalyst.T));
				for (auto s : rx->sub)
					for (auto r : s->reactants)
						rateNow *= state1(cell,r.ni);

				// Update creation and destruction arrays
				if (rateNow>0.0)
				{
					for (auto s : rx->sub)
					{
						tw::Float powerDensity = 0.0;
						tw::Float vibrationalPowerDensity = 0.0;
						tw::vec3 forceDensity = 0.0;

						// Direct loss of reactant mass, momentum, and energy
						for (auto r : s->reactants)
						{
							// ASSUMES NO MIXING OF VIBRATING AND NON-VIBRATING CHEMICALS IN EQUILIBRIUM GROUPS
							// (do not confuse the "group" of reactants with the EquilibriumGroup of a particular reactant)
							const tw::Float V = 1.0/(tw::small_pos + r.DensitySum(state1,cell)); // specific volume of reactant's EquilibriumGroup
							const tw::Float nFx = V*rateNow*state1(cell,r.npx);
							const tw::Float nFy = V*rateNow*state1(cell,r.npy);
							const tw::Float nFz = V*rateNow*state1(cell,r.npz);
							const tw::Float nP = V*rateNow*state1(cell,r.u);
							const tw::Float nPv = V*rateNow*state1(cell,r.x);

							DestroyMass(cell,rateNow,r);
							DestroyMomentum(cell,1,nFx,r);
							DestroyMomentum(cell,2,nFy,r);
							DestroyMomentum(cell,3,nFz,r);
							DestroyTotalEnergy(cell,nP,r);
							DestroyVibrations(cell,nPv,r);

							forceDensity += tw::vec3(nFx,nFy,nFz);
							powerDensity += nP;
							vibrationalPowerDensity += nPv;
						}
						// parcel out conserved quantities weighted by stoichiometric coefficients
						// ALL VIBRATIONAL ENERGY IS CONVERTED TO TRANSLATIONAL
						const tw::Float weight = 1.0/tw::Float(s->products.size());
						for (auto p : s->products)
						{
							CreateMass(cell,rateNow,p);
							CreateTotalEnergy(cell,weight*powerDensity,p);
							CreateMomentum(cell,1,weight*forceDensity.x,p);
							CreateMomentum(cell,2,weight*forceDensity.y,p);
							CreateMomentum(cell,3,weight*forceDensity.z,p);
						}
						// Heat of reaction
						// Affects the group of the last chemical in the reactant or product list
						// Formerly it was all partitioned into product groups as a signed creation term
						powerDensity = rateNow*(s->heat + s->vheat);
						vibrationalPowerDensity = rateNow*(s->vheat);
						if (powerDensity<0.0)
							DestroyTotalEnergy(cell,-powerDensity,s->reactants.back());
						else
							CreateTotalEnergy(cell,powerDensity,s->products.back());
						if (vibrationalPowerDensity<0.0)
							DestroyVibrations(cell,-vibrationalPowerDensity,s->reactants.back());
						else
							CreateVibrations(cell,vibrationalPowerDensity,s->products.back());
					}
				}
			}

		// COLLISIONS (MOMENTUM TRANSFER, THERMAL ENERGY TRANSFER)

		for (auto coll : collision)
		{
			LoadCollisionRate(coll,scratch);
			for (auto cell : InteriorCellRange(*this))
			{
				const tw::Float R = scratch(cell);
				const tw::Float N1 = state1(cell,coll->h1.ni);
				const tw::Float N2 = state1(cell,coll->h2.ni);
				const tw::Float T1 = eos1(cell,coll->e1.T);
				const tw::Float T2 = eos1(cell,coll->e2.T);
				const tw::Float m1 = coll->m1.mass;
				const tw::Float m2 = coll->m2.mass;
				const tw::Float m12 = m1*m2/(m1+m2);

				for (tw::Int ax=1;ax<=3;ax++)
				{
					const tw::Float n1v1 = state1(cell,coll->h1.npx+ax-1)/m1;
					const tw::Float n2v2 = state1(cell,coll->h2.npx+ax-1)/m2;
					CreateMomentum( cell , ax , m12*R*(N2*n1v1-N1*n2v2) , coll->h2);
					DestroyMomentum( cell , ax , m12*R*(N2*n1v1-N1*n2v2) , coll->h1);
				}
				if (T1>T2)
				{
					DestroyTotalEnergy( cell , 3.0*m12*R*N1*N2 * (T1 - T2) / (m1 + m2) , coll->h1);
					CreateTotalEnergy( cell ,  3.0*m12*R*N1*N2 * (T1 - T2) / (m1 + m2) , coll->h2);
				}
				else
				{
					CreateTotalEnergy( cell ,  3.0*m12*R*N1*N2 * (T2 - T1) / (m1 + m2) , coll->h1);
					DestroyTotalEnergy( cell , 3.0*m12*R*N1*N2 * (T2 - T1) / (m1 + m2) , coll->h2);
				}
			}
		}

		// VIBRATIONS

		for (auto x : excitation)
			for (auto cell : InteriorCellRange(*this))
			{
				const tw::Float Te = eos1(cell,x->e1.T);
				const tw::Float Tv = eos1(cell,x->e2.Tv);
				const tw::Float energy = x->m2.excitationEnergy;
				const tw::Float level = x->level;
				const tw::Float Xv = x->PrimitiveRate(std::fabs(Te));
				const tw::Int i1 = x->h1.ni;
				const tw::Int i2 = x->h2.ni;

				if (level>0)
				{
					tw::Float n0 = state1(cell,i2) * (1.0 - std::exp(-energy/Tv));
					rateNow = energy * level * Xv * state1(cell,i1) * n0 * (1.0 - std::exp(energy*level/Te - energy*level/Tv));
				}
				else
				{
					rateNow = energy * Xv * state1(cell,i1) * state1(cell,i2) * (1.0 - std::exp(energy/Te - energy/Tv));
				}

				CreateTotalAndVibrational(cell,rateNow,x->h2);
				DestroyTotalEnergy(cell,rateNow,x->h1);
			}
	}
}

void sparc::HydroManager::ComputeRadiativeSources()
{
	#pragma omp parallel
	{
		const tw::Float stef_boltz = 5.67e-8; // W/m^2/K^4
		const tw::vec3 R = owner->GlobalPhysicalSize().spatial();
		const tw::vec3 Rmks(R.x*tw::dims::length>>native>>mks,R.y*tw::dims::length>>native>>mks,R.z*tw::dims::length>>native>>mks);
		const tw::Float Vmks = owner->car*Rmks.x*Rmks.y*Rmks.z + owner->cyl*pi*sqr(Rmks.x)*Rmks.z + owner->sph*1.33*pi*cub(Rmks.x);
		const tw::Float Smks = owner->car*2*(Rmks.x*Rmks.y + Rmks.y*Rmks.z + Rmks.z*Rmks.x) + owner->cyl*(2*pi*sqr(Rmks.x)+2*pi*Rmks.x*Rmks.z) + owner->sph*4*pi*sqr(Rmks.x);

		// Ohmic heating due to static and laser fields
		if (electrons)
		{
			for (auto cell : InteriorCellRange(*this))
			{
				const tw::Float sig_AC = nu_e(cell)*state1(cell,ie)/(sqr(laserFrequency) + sqr(nu_e(cell)));
				const tw::Float sig_DC = state1(cell,ie)/nu_e(cell);
				const tw::Float E2_AC = 0.5*norm(laserAmplitude(cell));
				const tw::Float E2_DC = sqr((phi.fwd(cell,0,1) - phi.bak(cell,0,1))/owner->dL(cell,1)) +
					sqr((phi.fwd(cell,0,2) - phi.bak(cell,0,2))/owner->dL(cell,2)) +
					sqr((phi.fwd(cell,0,3) - phi.bak(cell,0,3))/owner->dL(cell,3));
				CreateTotalEnergy(cell,sig_AC*E2_AC,electrons->group->hidx);
				if (electrostaticHeating)
					CreateTotalEnergy(cell,sig_DC*E2_DC,electrons->group->hidx);
			}
		}

		// Photoionization
		for (auto cell : InteriorCellRange(*this))
			if (radiationIntensity(cell)>0.0)
				for (auto grp : group)
					for (auto chem : grp->chemical)
					{
						if (chem->ionizer!=NULL)
						{
							const hydro_set& r = chem->ionizer->hgas;
							const hydro_set& he = chem->ionizer->he;
							const hydro_set& hi = chem->ionizer->hi;

							const tw::Float Emag = std::sqrt(norm(laserAmplitude(cell)));
							const tw::Float photoRate = state1(cell,r.ni)*chem->ionizer->AverageRate(laserFrequency,Emag);
							const tw::Float V = 1.0/(tw::small_pos + r.DensitySum(state1,cell)); // specific volume of reactant's EquilibriumGroup
							const tw::Float nFx = V*photoRate*state1(cell,r.npx);
							const tw::Float nFy = V*photoRate*state1(cell,r.npy);
							const tw::Float nFz = V*photoRate*state1(cell,r.npz);
							const tw::Float nP = V*photoRate*state1(cell,r.u);
							const tw::Float nPv = V*photoRate*state1(cell,r.x);

							DestroyMass(cell,photoRate,r);
							DestroyMomentum(cell,1,nFx,r);
							DestroyMomentum(cell,2,nFy,r);
							DestroyMomentum(cell,3,nFz,r);
							DestroyTotalEnergy(cell,nP,r);
							DestroyVibrations(cell,nPv,r);
							CreateMass(cell,photoRate,hi);
							CreateMomentum(cell,1,nFx,hi);
							CreateMomentum(cell,2,nFy,hi);
							CreateMomentum(cell,3,nFz,hi);
							CreateTotalEnergy(cell,0.5*nP,hi); // half energy to ions
							CreateMass(cell,photoRate,he);
							CreateTotalEnergy(cell,0.5*nP,he); // half energy to electrons
						}
					}

		// Compute radiative losses
		// For optically thin, estimate mean free path from Zel'dovich table 5.2
		// Strictly these formulae only work in LTE, we hope effect is small when this is violated.
		// The effective LTE temperature is estimated from (total pressure / total density)
		if (radModel!=sparc::noRadiation)
		{
			for (auto cell : InteriorCellRange(*this))
			{
				tw::Float Ptot=tw::small_pos, ntot=tw::small_pos;
				for (auto grp : group)
				{
					Ptot += eos1(cell,grp->eidx.P);
					ntot += grp->DensitySum(state1,cell);
				}
				const tw::Float TKelvin = (Ptot/ntot)*tw::dims::temperature >> native >> mks;
				const tw::Float Lmks = tw::small_pos + 8.0e-14 * sqr(TKelvin); // mean free path in meters
				const tw::Float Imks = 4.0*stef_boltz*pow(TKelvin,4);
				if (radModel==sparc::thin)
					radiativeLosses(cell) = (Imks/Lmks)*tw::dims::power_density >> mks >> native;
				if (radModel==sparc::thick)
					radiativeLosses(cell) = (Imks*Smks/Vmks)*tw::dims::power_density >> mks >> native;
				for (auto grp : group)
				{
					const tw::Float lossNow = eos1(cell,grp->eidx.P)*radiativeLosses(cell)/Ptot;
					DestroyTotalEnergy(cell,lossNow,grp->hidx);
				}
			}
		}
	}
}

tw::vec3 sparc::HydroManager::ComputeForceOnBody(tw::Int i,tw::Int j,tw::Int k)
{
	tw::vec3 ans;
	tw::Int s;
	EquilibriumGroup* g;

	for (s=0;s<group.size();s++)
	{
		g = group[s];
		const tw::Float Ax0 = owner->dS(i,j,k,1);
		const tw::Float Ax1 = owner->dS(i+1,j,k,1);
		const tw::Float Ay0 = owner->dS(i,j,k,2);
		const tw::Float Ay1 = owner->dS(i,j+1,k,2);
		const tw::Float Az0 = owner->dS(i,j,k,3);
		const tw::Float Az1 = owner->dS(i,j,k+1,3);

		const tw::Float Pc = eos1(i,j,k,g->eidx.P);

		if (fluxMask(i-1,j,k)==0.0 && fluxMask(i,j,k)==1.0)
			ans.x -= Pc*Ax0;

		if (fluxMask(i,j,k)==1.0 && fluxMask(i+1,j,k)==0.0)
			ans.x += Pc*Ax1;

		if (fluxMask(i,j-1,k)==0.0 && fluxMask(i,j,k)==1.0)
			ans.y -= Pc*Ay0;

		if (fluxMask(i,j,k)==1.0 && fluxMask(i,j+1,k)==0.0)
			ans.y += Pc*Ay1;

		if (fluxMask(i,j,k-1)==0.0 && fluxMask(i,j,k)==1.0)
			ans.z -= Pc*Az0;

		if (fluxMask(i,j,k)==1.0 && fluxMask(i,j,k+1)==0.0)
			ans.z += Pc*Az1;
	}
	return ans;
}

void sparc::HydroManager::ComputeHydroSources()
{
	#pragma omp parallel
	{
		for (auto g : group)
		{
			if (g->chemical[0]==electrons || g->forceFilter==0.0)
			{
				// electrons are not advanced using hydro
				// charge and current are computed implicitly in field advance
				// momentum density is fixed to give the heavy particle velocity (see ApplyEOS)
				// the velocity is then used to convect the energy density
				// mass and momentum density are convected, but then reset to restore quasineutrality and heavy particle velocity

				// take away any collisional sources of momentum
				for (auto cell : InteriorCellRange(*this))
				{
					creationRate(cell,g->hidx.npx) = 0.0;
					creationRate(cell,g->hidx.npy) = 0.0;
					creationRate(cell,g->hidx.npz) = 0.0;
					destructionRate(cell,g->hidx.npx) = 0.0;
					destructionRate(cell,g->hidx.npy) = 0.0;
					destructionRate(cell,g->hidx.npz) = 0.0;
				}
			}
			else
			{
				for (tw::Int ax=1;ax<=3;ax++)
				{
					#pragma omp barrier
					g->LoadVelocity(scratch,state1,ax);
					#pragma omp barrier
					for (auto cell : InteriorCellRange(*this))
					{
						tw::Float dV,dS0,dS1,dl0,dl1,P0,P1,v0,v1;
						tw::Float forceDensity = 0.0;
						tw::Float powerDensity = 0.0;
						owner->GetCellMetrics(cell,ax,&dV,&dS0,&dS1,&dl0,&dl1);

						const tw::Float E1 = (phi.bak(cell,0,ax) - phi.fwd(cell,0,ax)) / (dl0 + dl1);
						const tw::Float nm = g->DensityWeightedSum(state1,g->matset.mass,cell);
						const tw::Float nq = g->DensityWeightedSum(state1,g->matset.charge,cell);
						const tw::Float Pc = eos1(cell,g->eidx.P);
						const tw::Float vc = scratch(cell);
						const tw::Float f0 = fluxMask.bak(cell,0,ax);
						const tw::Float fc = fluxMask(cell);
						const tw::Float f1 = fluxMask.fwd(cell,0,ax);

						P0 = eos1.bak(cell,g->eidx.P,ax);
						P1 = eos1.fwd(cell,g->eidx.P,ax);
						v0 = scratch.bak(cell,0,ax);
						v1 = scratch.fwd(cell,0,ax);

						P0 = (P0*f0 + Pc*fc) / (f0 + fc + tw::small_pos);
						P1 = (P1*f1 + Pc*fc) / (f1 + fc + tw::small_pos);
						v0 = 0.5*f0*fc*(v0 + vc);
						v1 = 0.5*f1*fc*(v1 + vc);

						// external forces

						forceDensity += nq * E1;
						powerDensity += nq * (vc * E1);

						// pressure force

						forceDensity += (dS0*P0 - dS1*P1)/dV;

						// work done by pressure

						powerDensity += (v0*P0*dS0 - v1*P1*dS1)/dV;

						// N3 force from walls

						if (f0==0.0 && fc==1.0 && vc<0.0)
							forceDensity += nm*vc*vc/dl0;

						if (fc==1.0 && f1==0.0 && vc>0.0)
							forceDensity -= nm*vc*vc/dl1;

						// put everything into the sources

						CreateMomentum(cell,ax,g->forceFilter*fc*forceDensity,g->hidx);
						if (powerDensity>0.0)
							CreateTotalEnergy(cell,fc*powerDensity,g->hidx);
						else
							DestroyTotalEnergy(cell,-fc*powerDensity,g->hidx);
					}
				} // end loop over axes

				// Undifferentiated tensor divergence terms
				//#pragma omp barrier
				if (owner->gridGeometry==tw::grid::cylindrical)
					for (auto cell : InteriorCellRange(*this))
					{
						const tw::Float nm = g->DensityWeightedSum(state1,g->matset.mass,cell);
						const tw::Float Pc = eos1(cell,g->eidx.P);
						const tw::vec3 vc = g->Velocity(state1,cell);
						const tw::vec3 pos = owner->Pos(cell);
						CreateMomentum(cell,1,fluxMask(cell)*(nm*sqr(vc.y) + Pc)/pos.x,g->hidx);
						DestroyMomentum(cell,2,fluxMask(cell)*nm*vc.x*vc.y/pos.x,g->hidx);
					}
				if (owner->gridGeometry==tw::grid::spherical)
					for (auto cell : InteriorCellRange(*this))
					{
						const tw::Float nm = g->DensityWeightedSum(state1,g->matset.mass,cell);
						const tw::Float Pc = eos1(cell,g->eidx.P);
						const tw::vec3 vc = g->Velocity(state1,cell);
						const tw::vec3 pos = owner->Pos(cell);
						const tw::Float tanz = std::tan(pos.z);
						CreateMomentum(cell,1,fluxMask(cell)*(nm*(sqr(vc.y) + sqr(vc.z)) + 2.0*Pc)/pos.x,g->hidx);
						DestroyMomentum(cell,2,fluxMask(cell)*(nm*(vc.x*vc.z - vc.y*vc.y/tanz) - Pc/tanz)/pos.x,g->hidx);
						DestroyMomentum(cell,3,fluxMask(cell)*nm*(vc.x*vc.y + vc.y*vc.z/tanz)/pos.x,g->hidx);
					}

			} // end else
		} // end loop over groups
	} // end parallel region
}

void sparc::HydroManager::ComputeSources()
{
	creationRate = 0.0;
	destructionRate = 0.0;

	ComputeCollisionalSources();
	ComputeRadiativeSources();
	ComputeHydroSources();

	creationRate.CopyFromNeighbors();
	creationRate.ApplyBoundaryCondition();

	destructionRate.CopyFromNeighbors();
	destructionRate.ApplyBoundaryCondition();
}

void sparc::HydroManager::LaserAdvance(tw::Float dt)
{
	// Currently SPARC ignores differences in frequency between injected pulses.
	// The average frequency of all the pulses is used throughout (see also HydroManager::Initialize).

	auto vacuumE = [&] (tw::Complex A,tw::Float k,tw::Float z)
	{
		// take enveloped vector potential and get space resolved electric field
		// only works in vacuum (e.g., at boundaries)
		return std::exp(ii*k*z)*ii*std::fabs(k)*A;
	};

	if (wave.size() && lasModel==sparc::vacuum) // use a prescribed field
	{
		#pragma omp parallel
		{
			for (auto cell : InteriorCellRange(*this))
			{
				tw::Complex fwd = 0.0;
				tw::Complex bak = 0.0;
				const tw::vec3 pos = owner->Pos(cell);
				for (auto pulse : wave)
				{
					if (pulse->direction.z > 0.0)
						fwd += pulse->VectorPotentialEnvelope(owner->WindowPos(0),pos,laserFrequency);
					else
						bak += pulse->VectorPotentialEnvelope(owner->WindowPos(0),pos,laserFrequency);
				}
				laserAmplitude(cell) = vacuumE(fwd,laserFrequency,pos.z) + vacuumE(bak,-laserFrequency,pos.z);
				radiationIntensity(cell) = 0.5*norm(laserAmplitude(cell));
			}
		}
	}

	if (electrons && wave.size() && lasModel==sparc::isotropic && dim[3]>1)
	{
		#pragma omp parallel
		{
			tw::Complex a0,a1; // amplitudes of incoming waves on left
			tw::Complex aN,aN1; // amplitudes of incoming waves on right
			for (auto strip : StripRange(*this,3,strongbool::no))
			{
				a0 = a1 = aN = aN1 = 0.0;
				const tw::Float z0 = owner->Pos(strip,0).z;
				const tw::Float z1 = owner->Pos(strip,1).z;
				const tw::Float zN = owner->Pos(strip,dim[3]).z;
				const tw::Float zN1 = owner->Pos(strip,dim[3]+1).z;
				for (auto pulse : wave)
				{
					if (pulse->direction.z > 0.0)
					{
						a0 += pulse->VectorPotentialEnvelope(owner->WindowPos(0),z0,laserFrequency);
						a1 += pulse->VectorPotentialEnvelope(owner->WindowPos(0),z1,laserFrequency);
					}
					else
					{
						aN += pulse->VectorPotentialEnvelope(owner->WindowPos(0),zN,laserFrequency);
						aN1 += pulse->VectorPotentialEnvelope(owner->WindowPos(0),zN1,laserFrequency);
					}
				}
				a0 = vacuumE(a0,laserFrequency,z0); a1 = vacuumE(a1,laserFrequency,z1);
				aN = vacuumE(aN,-laserFrequency,zN); aN1 = vacuumE(aN1,-laserFrequency,zN1);
				laserPropagator->SetupIncomingWaveLeft(strip,laserAmplitude,a0,a1,laserFrequency);
				laserPropagator->SetupIncomingWaveRight(strip,laserAmplitude,aN,aN1,laserFrequency);

				for (tw::Int k=1;k<=dim[3];k++)
				{
					refractiveIndex(strip,k) = 0.0;
					// Add dispersionless part of susceptibility
					for (auto grp : group)
						for (auto chem : grp->chemical)
						{
							tw::Complex susceptibility(chem->mat.eps[0] - 1.0,chem->mat.eps[1]);
							refractiveIndex(strip,k) += state1(strip,k,chem->indexInState) * susceptibility;
						}
					// Add plasma contribution to susceptibility
					refractiveIndex(strip,k) -= state1(strip,k,ie)/sqr(laserFrequency)*(one - ii*nu_e(strip,k)/laserFrequency)/(one + sqr(nu_e(strip,k)/laserFrequency));
					// Convert susceptibility to refractive index
					refractiveIndex(strip,k) = std::sqrt(one + refractiveIndex(strip,k));
				}
			}
		}

		laserPropagator->Advance(laserAmplitude,refractiveIndex,nu_e,laserFrequency,dt);

		#pragma omp parallel
		{
			for (auto cell : InteriorCellRange(*this))
				radiationIntensity(cell) = real(refractiveIndex(cell))*0.5*norm(laserAmplitude(cell));
		}
	}
}

tw::Float sparc::HydroManager::EstimateTimeStep()
{
	std::vector<tw::Float> dtMax(tw::GetOMPMaxThreads());
	tw::Float dtMaxAllThreads,dtMaxAllNodes;

	#pragma omp parallel
	{
		const tw::Int tid = tw::GetOMPThreadNum();
		const tw::Float sqrt_eps = std::sqrt(epsilonFactor);
		const tw::Float creationDominance = 100.0;
		tw::Float dts;

		if (owner->IsFirstStep())
			dtMax[tid] = owner->dx(0);
		else
			dtMax[tid] = owner->MaxSpacing(0);

		// Courant condition

		for (auto cell : InteriorCellRange(*this))
			for (tw::Int s=0;s<group.size();s++)
			{
				const tw::vec3 vel = group[s]->Velocity(state1,cell);
				if (Norm(vel)!=0.0)
				{
					const tw::Float dxi2 = sqr(0.5*owner->dL(cell,1));
					const tw::Float deta2 = sqr(0.5*owner->dL(cell,2));
					const tw::Float dzeta2 = sqr(0.5*owner->dL(cell,3));
					dts = sqr(vel.x)/dxi2 + sqr(vel.y)/deta2 + sqr(vel.z)/dzeta2;
					dts = 0.9/std::sqrt(dts);
					if (dts < dtMax[tid])
					{
						dtMax[tid] = dts;
						//statusMessage.str("");
						//statusMessage << "Limited by Courant condition : " << group[s]->name << " v = " << Magnitude(vel) << std::endl;
					}
				}
			}

		// Reaction rates, collision rates, other source terms

		auto AsymptoticStep = [&] (const tw::cell& cell,tw::Int c)
		{
			tw::Float dt_temp;
			if (creationRate(cell,c) > creationDominance*destructionRate(cell,c))
				dt_temp = sqrt_eps*std::fabs( (tw::small_pos+state1(cell,c)) / destructionRate(cell,c) );
			else
				dt_temp = sqrt_eps*std::fabs( (tw::small_pos+state1(cell,c)) / (creationRate(cell,c) - destructionRate(cell,c)) );
			if (dt_temp < dtMax[tid])
				return dt_temp;
			else
				return dtMax[tid];
		};

		for (auto g : group)
		{
			for (auto chem : g->chemical)
				for (auto cell : InteriorCellRange(*this))
					dtMax[tid] = AsymptoticStep(cell,chem->indexInState);

			for (auto cell : InteriorCellRange(*this))
				dtMax[tid] = AsymptoticStep(cell,g->hidx.u);

			for (auto cell : InteriorCellRange(*this))
				dtMax[tid] = AsymptoticStep(cell,g->hidx.x);
		}
	} // end parallel region

	// Choose the smallest maximum step from all the threads
	dtMaxAllThreads = *std::min_element(std::begin(dtMax),std::end(dtMax));
	// Choose the smallest maximum step from all the nodes
	dtMaxAllNodes = owner->strip[0].GetMin(dtMaxAllThreads);

	// Ensure that step size is sufficently small to resolve each laser pulse.
	tw::Float dtMaxAllPulses = owner->MaxSpacing(0);
	for (auto pulse : wave)
	{
		const PulseShape& shape = pulse->pulseShape;
		const tw::Float pulseDuration = shape.t4 - shape.t1;
		// For upcoming pulses...
		if ( owner->WindowPos(0) < shape.delay )
		{
			// Find closest pulse
			if ( dtMaxAllPulses > shape.delay - owner->WindowPos(0) )
				dtMaxAllPulses = shape.delay - owner->WindowPos(0);
		}
		// While inside a pulse
		if (owner->WindowPos(0)>=shape.delay && owner->WindowPos(0)<shape.delay+pulseDuration)
		{
			// Ensure Step resolves shortest pulse
			if ( dtMaxAllPulses > std::sqrt(epsilonFactor)*pulseDuration )
				dtMaxAllPulses = std::sqrt(epsilonFactor)*pulseDuration;
		}
		if ( dtMaxAllNodes > dtMaxAllPulses )
			dtMaxAllNodes = dtMaxAllPulses;
	}

	// If the critical time step is reached, switch to fixed time step and hope for the best!
	if (dtMaxAllNodes < owner->CriticalSpacing(0))
	{
		dtMaxAllNodes = owner->dx(0); // use the starting time step
		owner->adaptiveTimestep = false; // no more adaptive stepping after this
	}
	// Don't let it fall below the minimum time step
	return dtMaxAllNodes < owner->MinSpacing(0) ? owner->MinSpacing(0) : dtMaxAllNodes;
}

void sparc::HydroManager::DiffusionAdvance(tw::Float dt)
{
	for (auto g : group)
	{
		// HEAT CONDUCTION

		for (auto c : conductor)
			parabolicSolver->FixTemperature(eos1,Element(g->eidx.T),c->theRgn,c->Temperature(owner->WindowPos(0)));
		g->LoadMassDensity(scratch,state1);
		CopyFieldData(scratch2,Element(0),eos1,Element(g->eidx.T));
		parabolicSolver->Advance(eos1,g->eidx.T,fluxMask,&eos1,g->eidx.nmcv,&eos1,g->eidx.K,dt);
		g->eosMixData->UpdateEnergy(scratch,scratch2,state1,eos1);

		// VISCOSITY

		for (tw::Int ax=1;ax<=3;ax++)
		{
			g->LoadVelocity(scratch2,state1,ax);
			CopyBoundaryConditions(scratch2,0,state1,g->hidx.npx+ax-1);
			parabolicSolver->Advance(scratch2,0,fluxMask,&scratch,0,&eos1,g->eidx.visc,dt);
			#pragma omp parallel
			{
				for (auto cell : InteriorCellRange(*this))
					state1(cell,g->hidx.npx+ax-1) = scratch(cell)*scratch2(cell);
			}
		}
	}

	state1.CopyFromNeighbors();
	state1.ApplyBoundaryCondition();
	eos1.CopyFromNeighbors();
	eos1.ApplyBoundaryCondition();
}

void sparc::HydroManager::FieldAdvance(tw::Float dt)
{
	if (!electrons)
		return;

	const tw::Float q0 = electrons->mat.charge;
	const tw::Float m0 = electrons->mat.mass;
	const tw::Int Pe = electrons->group->eidx.P;

	// set up source and coefficient
	if (plasModel==sparc::quasineutral)
	{
		#pragma omp parallel
		{
			tw::Float D1,D2,P0,P1,mu0,mu1,dV,dS0,dS1,dl0,dl1;
			for (auto cell : InteriorCellRange(*this))
			{
				rho(cell) = rho0(cell);
				for (tw::Int ax=1;ax<=3;ax++)
				{
					owner->GetCellMetrics(cell,ax,&dV,&dS0,&dS1,&dl0,&dl1);
					P0 = eos1.bak(cell,Pe,ax);
					P1 = eos1.fwd(cell,Pe,ax);
					mu0 = 2.0*(q0/m0)/(nu_e(cell)+nu_e.bak(cell,0,ax));
					mu1 = 2.0*(q0/m0)/(nu_e(cell)+nu_e.fwd(cell,0,ax));
					// Laplacian coefficients weighted by mobility
					// ( for calculation of div(mu*grad(P)) )
					D1 = mu0*dS0/(dl0*dV);
					D2 = mu1*dS1/(dl1*dV);
					// re-use rho to hold rho_eff  = rho(t=0) + dt*div(mu*grad(P))
					rho(cell) += dt*(D1*P0 - (D1+D2)*eos1(cell,Pe) + D2*P1);
				}
				// coefficients multiplying grad(phi) go into scratch array
				scratch(cell) = 1.0 + dt*q0*state1(cell,ie)*(q0/m0)/nu_e(cell);
			}
		}
		scratch.CopyFromNeighbors();
		scratch.ApplyBoundaryCondition();
	}
	else
	{
		rho = 0.0;
		scratch = 1.0;
	}

	// Solve the elliptical equation --- div(scratch*grad(phi)) = -rho_eff
	// Even if plasma model is neutral, still do this to allow for external fields
	for (auto c : conductor)
		ellipticSolver->FixPotential(phi,c->theRgn,c->Voltage(owner->WindowPos(0)));
	ellipticSolver->SetCoefficients(&scratch);
	ellipticSolver->Solve(phi,rho,-1.0);

	// Compute the new charge density
	// TODO: #pragma omp parallel for collapse(3) schedule(static)
	for (tw::Int i=1;i<=dim[1];i++)
		for (tw::Int j=1;j<=dim[2];j++)
			for (tw::Int k=1;k<=dim[3];k++)
			{
				rho(i,j,k) = (phi(i,j,k) - phi(i-1,j,k))*owner->dS(i,j,k,1) / owner->dl(i,j,k,1);
				rho(i,j,k) += (phi(i,j,k) - phi(i+1,j,k))*owner->dS(i+1,j,k,1) / owner->dl(i+1,j,k,1);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j-1,k))*owner->dS(i,j,k,2) / owner->dl(i,j,k,2);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j+1,k))*owner->dS(i,j+1,k,2) / owner->dl(i,j+1,k,2);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j,k-1))*owner->dS(i,j,k,3) / owner->dl(i,j,k,3);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j,k+1))*owner->dS(i,j,k+1,3) / owner->dl(i,j,k+1,3);
				rho(i,j,k) /= owner->dS(i,j,k,0);
			}
	rho.CopyFromNeighbors();
	rho.ApplyBoundaryCondition();
}

void sparc::HydroManager::HydroAdvance(const tw::grid::axis& axis,tw::Float dt)
{
	// Convect the fluid

	tw::Int ax = tw::grid::naxis(axis);

	if (dim[ax] > 1)
	{
		FCT_Driver convector(&state0,&state1,&scratch,&fluxMask,owner);
		tw::bc::fld bc0 = (owner->bc0[ax] == par::reflecting || owner->bc0[ax] == par::axisymmetric) ? fld::dirichletCell : fld::neumannWall;
		tw::bc::fld bc1 = owner->bc1[ax] == par::reflecting ? fld::dirichletCell : fld::neumannWall;

		for (auto g : group)
		{
			convector.SetDensityElements(Element(g->hidx.first,g->hidx.last));
			convector.SetVelocityElement(0);
			#pragma omp parallel
			{
				g->LoadVelocity(scratch,state1,ax);
				for (auto cell : EntireCellRange(*this))
					scratch(cell) *= fluxMask(cell);
			}
			convector.Convect(axis,bc0,bc1,dt);
		}
		state0.ApplyBoundaryCondition();
	}
}

void sparc::HydroManager::ChemAdvance(tw::Float dt)
{
	creationRate *= dt;
	state0 += creationRate;

	destructionRate *= dt;
	state0 -= destructionRate;

	state0.CopyFromNeighbors();
	state0.ApplyBoundaryCondition();
}

void sparc::HydroManager::EOSAdvance(tw::Float dt)
{
	// Load (P,T,Tv,K,visc) into eos using (n,np,u,x) from state.
	// dt = 0.0 signals that we want to use a method that involves only 1 time level.
	// Otherwise we are allowed to use different time levels per the chain rule.
	// E.g., d(u/nm)/dT = cv --> d(u/nm)/dt = cv*dT/dt
	// Here we also impose quasineutrality and fix electron velocity
	// Finally, we throw an error if numerical failure is detected

	// Time centering information upon entry*:
	// Trial Step:    Full Step:
	// qty    time    qty    time
	// ---    ----    ---    ----
	// state0 0      state0  1/2
	// state1 1/2    state1  1
	// eos0   0      eos0    0
	// eos1   0      eos1    1/2
	// * notice eos1 and state0 are always known at the same time.
	//   and that state1 is always 1/2 step after state0.
	//   We are trying to load eos0 with the time level of state1.

	// Time centering information upon output*:
	// Trial Step:    Full Step:
	// qty    time    qty    time
	// ---    ----    ---    ----
	// state0 0      state0  1/2
	// state1 1/2    state1  1
	// eos0   1/2    eos0    1
	// eos1   0      eos1    1/2
	// * eos0 and eos1 are swapped after exiting function

	// Quasineutrality and velocity forcing
	#pragma omp parallel
	{
		tw::Float ionChargeDensity;
		tw::vec3 ionVelocity;

		for (auto cell : EntireCellRange(*this))
		{
			ionChargeDensity = 0.0;
			ionVelocity = 0.0;
			for (auto g : group)
			{
				if (g->chemical[0]!=electrons)
				{
					ionChargeDensity += g->DensityWeightedSum(state1,g->matset.charge,cell);
					ionVelocity += g->Velocity(state1,cell);
				}
			}
			if (electrons)
			{
				// Impose quasineutrality and assume electron velocity = average heavy particle velocity
				const tw::Float ne = -ionChargeDensity/electrons->mat.charge;
				state1(cell,ie) = ne;
				state1(cell,electrons->group->hidx.npx) = ne*electrons->mat.mass*ionVelocity.x/tw::Float(group.size()-1);
				state1(cell,electrons->group->hidx.npy) = ne*electrons->mat.mass*ionVelocity.y/tw::Float(group.size()-1);
				state1(cell,electrons->group->hidx.npz) = ne*electrons->mat.mass*ionVelocity.z/tw::Float(group.size()-1);
			}
		}
	}

	// DFG - factorized EOS loop, avoids nested tools.
	eos0 = 0.0; // safest to explicitly reset here
	for (auto g : group)
	{
		for (auto chem : g->chemical)
			chem->eosData->AddHeatCapacity(state1,eos0);
		// The following loads T into eos, IE into scratch, and nm into scratch2
		// N.b. eos0 is sought at the time of state1, eos1 is known at the time of state0.
		if (dt==0.0)
			g->eosMixData->ComputeTemperature(scratch,scratch2,state1,eos0);
		else
			g->eosMixData->ComputeTemperature(scratch,scratch2,state0,state1,eos1,eos0);
		for (auto chem : g->chemical)
			chem->eosData->AddPKV(scratch,scratch2,nu_e,state1,eos0);
	}

	// Check for numerical failure, defined by NaN in the hydro state vector
	tw::Int badCells = 0;
	for (auto cell : EntireCellRange(*this))
		for (tw::Int c=0;c<state1.Components();c++)
			badCells += std::isnan(state1(cell,c));
	if (badCells)
		throw tw::FatalError("Encountered NaN in hydrodynamic state");

	// DFG - EntireCellRange has been updated, no need for message passing.
	// state1.CopyFromNeighbors();
	// state1.ApplyBoundaryCondition();
	// eos0.CopyFromNeighbors();
	// eos0.ApplyBoundaryCondition();
}

void sparc::HydroManager::FirstOrderAdvance(tw::Float dt,bool computeSources)
{
	if (computeSources)
	{
		ComputeElectronCollisionFrequency();
		ComputeSources();
	}
	HydroAdvance(tw::grid::x,dt);
	HydroAdvance(tw::grid::y,dt);
	HydroAdvance(tw::grid::z,dt);
	ChemAdvance(dt);
	Swap(state0,state1);

	EOSAdvance(dt);
	Swap(eos0,eos1);

	DiffusionAdvance(dt);
	FieldAdvance(dt);
}

void sparc::HydroManager::Update()
{
	const tw::Float dt = dx(0);
	tw::Float dts;

	if (!owner->adaptiveTimestep)
	{
		FirstOrderAdvance(0.5*dt,true);
		FirstOrderAdvance(dt,true);
	}
	else
	{
		ComputeElectronCollisionFrequency();
		ComputeSources();
		dts = EstimateTimeStep();
		owner->UpdateTimestep(dts);
		FirstOrderAdvance(0.5*dt,false);
		FirstOrderAdvance(dt,true);
	}

	LaserAdvance(dt);
}

bool sparc::HydroManager::ReadQuasitoolBlock(const TSTreeCursor *curs0,const std::string& src)
{
	std::string key = tw::input::node_kind(curs0);
	auto curs = tw::input::Cursor(curs0);
	if (key=="reaction")
	{
		reaction.push_back(new Reaction);
		reaction.back()->ReadInputFile(curs.get(),src,native);
		return true;
	}
	if (key=="excitation")
	{
		excitation.push_back(new Excitation);
		excitation.back()->ReadInputFile(curs.get(),src,native);
		return true;
	}
	if (key=="collision")
	{
		collision.push_back(new Collision);
		collision.back()->ReadInputFile(curs.get(),src,native);
		return true;
	}
	return false;
}

void sparc::HydroManager::ReadCheckpoint(std::ifstream& inFile)
{
	Module::ReadCheckpoint(inFile);
	eos1.ReadCheckpoint(inFile);
	state1.ReadCheckpoint(inFile);
	rho.ReadCheckpoint(inFile);
	phi.ReadCheckpoint(inFile);
}

void sparc::HydroManager::WriteCheckpoint(std::ofstream& outFile)
{
	Module::WriteCheckpoint(outFile);
	eos1.WriteCheckpoint(outFile);
	state1.WriteCheckpoint(outFile);
	rho.WriteCheckpoint(outFile);
	phi.WriteCheckpoint(outFile);
}

void sparc::HydroManager::Report(Diagnostic& diagnostic)
{
	Module::Report(diagnostic);

	std::map<tw::Int,std::string> xyz = {{1,"x"},{2,"y"},{3,"z"}};

	auto WriteSubmoduleData = [&] (Field& theData,tw::Int comp,const std::string& filename,const tw::dims units=tw::dims::none,const std::string& pretty="tw::none")
	{
		CopyFieldData(scratch,Element(0),theData,Element(comp));
		scratch *= fluxMask;
		diagnostic.ReportField(filename,scratch,0,units,pretty);
	};

	// Mass, Charge, Energy

	scratch = 0.0;
	for (auto g : group)
	{
		g->LoadMassDensity(scratch2,state1);
		scratch2 *= fluxMask;
		scratch += scratch2;
	}
	diagnostic.VolumeIntegral("Mass",scratch,0);
	diagnostic.ReportField("massdensity",scratch,0,tw::dims::mass_density,"$\\rho$");

	diagnostic.VolumeIntegral("Charge",rho,0);

	scratch = 0.0;
	for (auto g : group)
		for (auto cell : InteriorCellRange(*this))
			scratch(cell) += state1(cell,g->hidx.u);
	diagnostic.VolumeIntegral("Energy",scratch,0);

	// Momentum

	for (tw::Int ax=1;ax<=3;ax++)
	{
		scratch = 0.0;
		for (auto g : group)
			for (auto cell : InteriorCellRange(*this))
				scratch(cell) += state1(cell,g->hidx.npx+ax-1);
		diagnostic.VolumeIntegral("P"+xyz[ax],scratch,0);
	}

	// Dipole moment

	diagnostic.FirstMoment("Dx",rho,0,dipoleCenter,tw::grid::x);
	diagnostic.FirstMoment("Dy",rho,0,dipoleCenter,tw::grid::y);
	diagnostic.FirstMoment("Dz",rho,0,dipoleCenter,tw::grid::z);

	// Impermeable Region

	diagnostic.ReportField("impermeable",fluxMask,0);

	// Collision Diagnostic

	diagnostic.ReportField("collisionFreq",nu_e,0,tw::dims::frequency,"$\\nu_e$");

	// Radiation Diagnostic

	diagnostic.ReportField("rad-intensity",radiationIntensity,0,tw::dims::intensity,"$I$");
	diagnostic.ReportField("rad-ereal",laserAmplitude,0,tw::dims::electric_field,"$\\Re E$");
	diagnostic.ReportField("rad-eimag",laserAmplitude,1,tw::dims::electric_field,"$\\Im E$");
	diagnostic.ReportField("rad-losses",radiativeLosses,0,tw::dims::power_density,"${\\cal P}$");
	diagnostic.ReportField("rad-nreal",refractiveIndex,0);
	diagnostic.ReportField("rad-nimag",refractiveIndex,1);
	diagnostic.ReportField("me_eff",me_eff,0);

	// Electrostatics

	if (electrons)
	{
		diagnostic.ReportField("chem-phi",phi,0,tw::dims::scalar_potential,"$\\phi$");
		diagnostic.ReportField("chem-rho",rho,0,tw::dims::charge_density,"$\\rho$");
	}

	// Constituents and groups

	for (auto g : group)
	{
		for (auto chem : g->chemical)
			WriteSubmoduleData(state1,chem->indexInState,chem->name,tw::dims::density,"$n[$"+chem->name+"$]$");

		WriteSubmoduleData(eos1,g->eidx.T,"T_" + g->name,tw::dims::temperature,"$T$");
		WriteSubmoduleData(eos1,g->eidx.P,"P_" + g->name,tw::dims::pressure,"$P$");
		WriteSubmoduleData(eos1,g->eidx.nmcv,"nmcv_" + g->name,tw::dims::density,"$nmc_v$");
		WriteSubmoduleData(eos1,g->eidx.K,"K_" + g->name,tw::dims::thermal_conductivity,"$K$");
		WriteSubmoduleData(eos1,g->eidx.Tv,"Tv_" + g->name,tw::dims::temperature,"$T_v$");

		for (tw::Int ax=1;ax<=3;ax++)
		{
			g->LoadVelocity(scratch,state1,ax);
			scratch *= fluxMask;
			diagnostic.ReportField("v"+xyz[ax]+"_"+g->name,scratch,0,tw::dims::velocity,"$v_"+xyz[ax]+"$");
		}
	}
}

void sparc::HydroManager::StatusMessage(std::ostream *dest)
{
	std::println(*dest,"{}",statusMessage.str());
}
