module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_test.h"
#include "tw_logger.h"

export module hydro;
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

export import hydro_primitives;
export import :chemical;
export import :group;

using namespace tw::bc;

export namespace sparc
{
struct HydroManager:Driver
{
	std::vector<EquilibriumGroup*> group; // explicitly type submodule list
	std::vector<Reaction*> reaction;
	std::vector<Excitation*> excitation;
	std::vector<Collision*> collision;
	Field state0,state1,creationRate,destructionRate,eos0,eos1;

	tw::Waves waves;
	tw::Conductors conductors;
	std::shared_ptr<ParabolicSolver> parabolicSolver;
	std::shared_ptr<EllipticSolver> ellipticSolver;
	std::shared_ptr<IsotropicPropagator> laserPropagator;
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

	HydroManager(const std::string& name,MetricSpace *ms,Task *tsk);
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

sparc::HydroManager::HydroManager(const std::string& name,MetricSpace *ms, Task *tsk) : Driver(name,ms,tsk)
{
	if (native.unit_system!=tw::units::plasma)
		throw tw::FatalError("HydroManager module requires <native units = plasma>");

	epsilonFactor = 1e-4;
	laserFrequency = 1.0;
	backgroundDensity = backgroundTemperature = 0.0;

	// some arrays must be initialized later since we don't know how many elements they have yet
	scratch.Initialize(*space,task);
	scratch2.Initialize(*space,task);
	fluxMask.Initialize(*space,task);
	rho0.Initialize(*space,task);
	rho.Initialize(*space,task);
	phi.Initialize(*space,task);
	nu_e.Initialize(*space,task);
	me_eff.Initialize(*space,task);
	laserAmplitude.Initialize(*space,task);
	radiativeLosses.Initialize(*space,task);
	radiationIntensity.Initialize(*space,task);
	refractiveIndex.Initialize(*space,task);

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
	r = 0; // running index
	for (auto grp : group)
		r = grp->eidx.Load(r);

	// Hydro indexing
	r = 0; // running index
	for (auto grp : group)
		r = grp->hidx.Load(r,grp->chemical.size());

	// Constituent indexing, and find electrons
	for (auto grp : group) {
		for (auto i=0; i<grp->chemical.size(); i++) {
			grp->chemical[i]->indexInState = grp->hidx.first + i;
			if (grp->chemical[i]->mat.mass == 1.0) {
				ie = grp->hidx.first + i;
				electrons = grp->chemical[i];
				grp->forceFilter = 0.0;
			}
		}
	}

	// Pass indexing down to lower level objects.
	for (auto grp : group)
		grp->SetupIndexing();

	// Setup lambdas to help index microscopics (photoionization is handled by groups)

	auto GetHydroSet = [&] (const std::string& name)
	{
		sparc::hydro_set ans;
		Chemical *chem = (Chemical*)Root()->FindDriver(name,true);
		if (chem==NULL) {
			throw tw::FatalError(std::format("could not find {}",name));
		}
		ans = ((EquilibriumGroup*)chem->super)->hidx;
		ans.ni = chem->indexInState;
		return ans;
	};

	auto GetEOSSet = [&] (const std::string& name)
	{
		Chemical *chem = (Chemical*)Root()->FindDriver(name,true);
		if (chem==NULL) {
			throw tw::FatalError(std::format("could not find {}",name));
		}
		return ((EquilibriumGroup*)chem->super)->eidx;
	};

	auto GetMaterial = [&] (const std::string& name)
	{
		Chemical *chem = (Chemical*)Root()->FindDriver(name,true);
		if (chem==NULL) {
			throw tw::FatalError(std::format("could not find {}",name));
		}
		return chem->mat;
	};

	logger::TRACE("reaction indexing");

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
	logger::TRACE("collision indexing");

	for (auto coll : collision)
	{
		coll->h1 = GetHydroSet(coll->name1);
		coll->e1 = GetEOSSet(coll->name1);
		coll->h2 = GetHydroSet(coll->name2);
		coll->e2 = GetEOSSet(coll->name2);
		coll->m1 = GetMaterial(coll->name1);
		coll->m2 = GetMaterial(coll->name2);
	}
	logger::TRACE("excitation indexing");

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
	Driver::VerifyInput();

	if (backgroundDensity!=0.0)
	{
		if (backgroundDensity<=0.0)
			throw tw::FatalError("Hydro module background density must be positive.");
		if (backgroundTemperature<=0.0)
			throw tw::FatalError("Hydro module background temperature must be positive.");
	}

	// Search submodule list for EquilibriumGroup modules
	for (auto sub : sub_drivers)
	{
		EquilibriumGroup *grp = dynamic_cast<EquilibriumGroup*>(sub);
		if (grp!=NULL) {
			group.push_back(grp);
			grp->VerifyInput();
		}
	}

	// Find existing tools
	for (auto tool : tools) {
		if (std::dynamic_pointer_cast<EllipticSolver>(tool)) {
			ellipticSolver = std::dynamic_pointer_cast<EllipticSolver>(tool);
		} else if (std::dynamic_pointer_cast<ParabolicSolver>(tool)) {
			parabolicSolver = std::dynamic_pointer_cast<ParabolicSolver>(tool);
		} else if (std::dynamic_pointer_cast<IsotropicPropagator>(tool)) {
			laserPropagator = std::dynamic_pointer_cast<IsotropicPropagator>(tool);
		} else if (std::dynamic_pointer_cast<Wave>(tool)) {
			waves.push_back(std::dynamic_pointer_cast<Wave>(tool));
		} else if (std::dynamic_pointer_cast<Conductor>(tool)) {
			conductors.push_back(std::dynamic_pointer_cast<Conductor>(tool));
		}
	}

	// if tools are missing create a default tool
	if (!parabolicSolver) {
		auto new_tool = CreateTool("default_parabolic_solver",tw::tool_type::generalParabolicPropagator);
		AddTool(new_tool);
		parabolicSolver = std::dynamic_pointer_cast<ParabolicSolver>(new_tool);
	}
	if (!ellipticSolver) {
		auto new_tool = space->SpatialDims()==1 ?
			CreateTool("default_elliptic_solver",tw::tool_type::ellipticSolver1D) :
			CreateTool("default_elliptic_solver",tw::tool_type::iterativePoissonSolver);
		AddTool(new_tool);
		ellipticSolver = std::dynamic_pointer_cast<EllipticSolver>(new_tool);		
		// Default tool needs some default boundaries; if a tool was attached the user will have specified them.
		tw::bc::fld x0 = space->bc0[1]==par::axisymmetric ? fld::neumannWall : fld::dirichletCell;
		tw::bc::fld z0 = fld::dirichletCell;
		tw::bc::fld z1 = fld::neumannWall;
		ellipticSolver->SetBoundaryConditions(x0,z0,z0,z0,z0,z1);
	}
	if (!laserPropagator) {
		auto new_tool = CreateTool("default_laser_propagator",tw::tool_type::isotropicPropagator);
		AddTool(new_tool);
		laserPropagator = std::dynamic_pointer_cast<IsotropicPropagator>(new_tool);
	}
}

void sparc::HydroManager::Initialize()
{
	Driver::Initialize();

	logger::TRACE("initialize laser frequency and refractive index");
	if (waves.size())
	{
		laserFrequency = 0.0;
		for (auto pulse : waves)
			laserFrequency += pulse->w;
		laserFrequency /= tw::Float(waves.size());
	}
	refractiveIndex = tw::Complex(1.0,0.0);

	logger::TRACE("initialize electrostatic tools");
	ellipticSolver->SetFieldsBoundaryConditions(phi,Rng(0));
	rho.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	rho.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	rho.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);
	scratch.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	scratch.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	scratch.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);

	logger::TRACE("initialize reactions");
	// Find number of state variables characterizing the whole system
	// DFG - now we have a static variable to help keep count of these more reliably
	// For more see sparc::hydro_set and sparc::eos_set in physics.h
	tw::Int hydroElements = 0;
	for (auto grp : group)
		hydroElements += grp->chemical.size() + sparc::hydro_set::count;
	tw::Int eosElements = group.size()*sparc::eos_set::count;
	state0.Initialize(hydroElements,*space,task);
	state1.Initialize(hydroElements,*space,task);
	creationRate.Initialize(hydroElements,*space,task);
	destructionRate.Initialize(hydroElements,*space,task);
	eos0.Initialize(eosElements,*space,task);
	eos1.Initialize(eosElements,*space,task);

	logger::TRACE("initialize boundary conditions");
	// variables defining normal component boundary conditions
	tw::bc::fld bc0[4],bc1[4];
	for (tw::Int ax=1;ax<=3;ax++)
	{
		bc0[ax] = (space->bc0[ax] == par::reflecting || space->bc0[ax] == par::axisymmetric) ? fld::dirichletWall : fld::neumannWall;
		bc1[ax] = (space->bc1[ax] == par::reflecting || space->bc1[ax] == par::axisymmetric) ? fld::dirichletWall : fld::neumannWall;
	}

	logger::TRACE("setup impermeable regions");
	// setup the flux mask used to make conductors impermeable
	// also used for reflecting boundary conditions at simulation walls
	fluxMask = 1.0;
	for (auto cell : EntireCellRange(*space,1))
		for (auto c : conductors)
			if (c->Inside(space->Pos4(cell)))
				fluxMask(cell) = 0.0;
	logger::TRACE("setup global boundaries");
	for (tw::Int ax=1;ax<=3;ax++)
		for (auto strip : StripRange(*space,ax,0,1,strongbool::yes))
		{
			if (bc0[ax]==fld::dirichletWall && task->n0[ax]==MPI_PROC_NULL && dim[ax]>1)
				fluxMask(strip,lfg[ax]) = fluxMask(strip,lng[ax]) = 0.0;
			if (bc1[ax]==fld::dirichletWall && task->n1[ax]==MPI_PROC_NULL && dim[ax]>1)
				fluxMask(strip,ufg[ax]) = fluxMask(strip,ung[ax]) = 0.0;
		}

	logger::TRACE("setup basic boundary conditions");
	// Default boundary conditions, refine after setting up indexing.
	state0.SetBoundaryConditions(All(state0),tw::grid::x,fld::neumannWall,fld::neumannWall);
	state0.SetBoundaryConditions(All(state0),tw::grid::y,fld::neumannWall,fld::neumannWall);
	state0.SetBoundaryConditions(All(state0),tw::grid::z,fld::neumannWall,fld::neumannWall);
	state1.SetBoundaryConditions(All(state1),tw::grid::x,fld::neumannWall,fld::neumannWall);
	state1.SetBoundaryConditions(All(state1),tw::grid::y,fld::neumannWall,fld::neumannWall);
	state1.SetBoundaryConditions(All(state1),tw::grid::z,fld::neumannWall,fld::neumannWall);
	eos0.SetBoundaryConditions(All(eos0),tw::grid::x,fld::neumannWall,fld::neumannWall);
	eos0.SetBoundaryConditions(All(eos0),tw::grid::y,fld::neumannWall,fld::neumannWall);
	eos0.SetBoundaryConditions(All(eos0),tw::grid::z,fld::neumannWall,fld::neumannWall);
	eos1.SetBoundaryConditions(All(eos1),tw::grid::x,fld::neumannWall,fld::neumannWall);
	eos1.SetBoundaryConditions(All(eos1),tw::grid::y,fld::neumannWall,fld::neumannWall);
	eos1.SetBoundaryConditions(All(eos1),tw::grid::z,fld::neumannWall,fld::neumannWall);
	creationRate.SetBoundaryConditions(All(creationRate),tw::grid::x,fld::neumannWall,fld::neumannWall);
	creationRate.SetBoundaryConditions(All(creationRate),tw::grid::y,fld::neumannWall,fld::neumannWall);
	creationRate.SetBoundaryConditions(All(creationRate),tw::grid::z,fld::neumannWall,fld::neumannWall);
	destructionRate.SetBoundaryConditions(All(destructionRate),tw::grid::x,fld::neumannWall,fld::neumannWall);
	destructionRate.SetBoundaryConditions(All(destructionRate),tw::grid::y,fld::neumannWall,fld::neumannWall);
	destructionRate.SetBoundaryConditions(All(destructionRate),tw::grid::z,fld::neumannWall,fld::neumannWall);
	nu_e.SetBoundaryConditions(tw::grid::x,fld::neumannWall,fld::neumannWall);
	nu_e.SetBoundaryConditions(tw::grid::y,fld::neumannWall,fld::neumannWall);
	nu_e.SetBoundaryConditions(tw::grid::z,fld::neumannWall,fld::neumannWall);

	// Prepare for loading fluid: indexing, background profiles, and refined boundary conditions

	logger::TRACE("setup the hydro tree");
	SetupIndexing();

	logger::TRACE("setup automatic background vapor");
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

	logger::TRACE("refine boundary conditions by group");
	for (auto grp : group)
	{
		Rng e;
		// Temperature
		e = Rng(grp->eidx.T);
		parabolicSolver->SetFieldsBoundaryConditions(eos0,e);
		parabolicSolver->SetFieldsBoundaryConditions(eos1,e);
		// X-Component
		e = Rng(grp->hidx.npx);
		state0.SetBoundaryConditions(e,tw::grid::x,bc0[1],bc1[1]);
		state1.SetBoundaryConditions(e,tw::grid::x,bc0[1],bc1[1]);
		creationRate.SetBoundaryConditions(e,tw::grid::x,bc0[1],bc1[1]);
		destructionRate.SetBoundaryConditions(e,tw::grid::x,bc0[1],bc1[1]);
		// Y-Component
		e = Rng(grp->hidx.npy);
		state0.SetBoundaryConditions(e,tw::grid::y,bc0[2],bc1[2]);
		state1.SetBoundaryConditions(e,tw::grid::y,bc0[2],bc1[2]);
		creationRate.SetBoundaryConditions(e,tw::grid::y,bc0[2],bc1[2]);
		destructionRate.SetBoundaryConditions(e,tw::grid::y,bc0[2],bc1[2]);
		// Z-Component
		e = Rng(grp->hidx.npz);
		state0.SetBoundaryConditions(e,tw::grid::z,bc0[3],bc1[3]);
		state1.SetBoundaryConditions(e,tw::grid::z,bc0[3],bc1[3]);
		creationRate.SetBoundaryConditions(e,tw::grid::z,bc0[3],bc1[3]);
		destructionRate.SetBoundaryConditions(e,tw::grid::z,bc0[3],bc1[3]);
	}

	logger::TRACE("create initial fluid constituents");
	for (auto grp : group)
	{
		grp->forceFilter = grp->mobile ? 1.0 : 0.0;
		grp->GenerateFluid(state1,eos1,scratch,scratch2);
	}

	me_eff = 1.0; // effective mass reserved for future use
	nu_e = 1.0; // put arbitrary value so initial transport coefficients don't blow up

	logger::TRACE("initialize EOS");
	EOSAdvance(0.0); // gets eos0 using state1 only
	eos1 = eos0;
	ComputeElectronCollisionFrequency();
}

void sparc::HydroManager::Reset()
{
	state0 = state1;
	eos0 = eos1;
	rho0 = rho;

	if (!space->IsFirstStep())
	{
		// Generate new fluid due to sources
		// N.b. EquilibriumGroup::GenerateFluid is destructive to the eos field that is passed in.
		bool didGenerate = false;
		for (auto grp : group) {
			didGenerate |= grp->GenerateFluid(state1,eos1,scratch,scratch2);
		}
		if (didGenerate)
		{
			EOSAdvance(space->dX(1,0)); // gets eos0 using state0 and state1
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
		for (auto cell : InteriorCellRange(*this,1))
		{
			T1 = eos1(cell,coll->e1.T) * econv;
			T2 = eos1(cell,coll->e2.T) * econv;
			v12 = std::sqrt(8.0*(T1/m1 + T2/m2)/pi);
			R(cell) = rconv * (4.0/3.0) * v12 * sigma;
		}

	if (coll->type==sparc::coulomb)
		for (auto cell : InteriorCellRange(*this,1))
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
		for (auto cell : InteriorCellRange(*this,1))
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
		for (auto cell : InteriorCellRange(*this,1))
			nu_e(cell) = 0.0;
		for (auto coll : collision)
		{
			LoadCollisionRate(coll,scratch);
			if (coll->h1.ni==ie)
				for (auto cell : InteriorCellRange(*this,1))
					nu_e(cell) += state1(cell,coll->h2.ni) * scratch(cell);
			if (coll->h2.ni==ie)
				for (auto cell : InteriorCellRange(*this,1))
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
			for (auto cell : InteriorCellRange(*this,1))
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
			for (auto cell : InteriorCellRange(*this,1))
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
			for (auto cell : InteriorCellRange(*this,1))
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
		const tw::vec3 R = space->GlobalPhysicalSize().spatial();
		const tw::vec3 Rmks(R.x*tw::dims::length>>native>>mks,R.y*tw::dims::length>>native>>mks,R.z*tw::dims::length>>native>>mks);
		const tw::Float Vmks = space->car*Rmks.x*Rmks.y*Rmks.z + space->cyl*pi*sqr(Rmks.x)*Rmks.z + space->sph*1.33*pi*cub(Rmks.x);
		const tw::Float Smks = space->car*2*(Rmks.x*Rmks.y + Rmks.y*Rmks.z + Rmks.z*Rmks.x) + space->cyl*(2*pi*sqr(Rmks.x)+2*pi*Rmks.x*Rmks.z) + space->sph*4*pi*sqr(Rmks.x);

		// Ohmic heating due to static and laser fields
		if (electrons)
		{
			auto hidx = ((EquilibriumGroup*)electrons->super)->hidx;
			for (auto cell : InteriorCellRange(*this,1))
			{
				const tw::Float sig_AC = nu_e(cell)*state1(cell,ie)/(sqr(laserFrequency) + sqr(nu_e(cell)));
				const tw::Float sig_DC = state1(cell,ie)/nu_e(cell);
				const tw::Float E2_AC = 0.5*norm(laserAmplitude(cell));
				const tw::Float E2_DC = sqr((phi.fwd(cell,0,1) - phi.bak(cell,0,1))/space->dL(cell,1)) +
					sqr((phi.fwd(cell,0,2) - phi.bak(cell,0,2))/space->dL(cell,2)) +
					sqr((phi.fwd(cell,0,3) - phi.bak(cell,0,3))/space->dL(cell,3));
				CreateTotalEnergy(cell,sig_AC*E2_AC,hidx);
				if (electrostaticHeating)
					CreateTotalEnergy(cell,sig_DC*E2_DC,hidx);
			}
		}

		// Photoionization
		for (auto cell : InteriorCellRange(*this,1))
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
			for (auto cell : InteriorCellRange(*this,1))
			{
				tw::Float Ptot=tw::small_pos, ntot=tw::small_pos;
				for (auto grp : group)
				{
					Ptot += eos1(cell,grp->eidx.P);
					ntot += grp->DensitySum(state1,cell);
				}
				const tw::Float TKelvin = (Ptot/ntot)*tw::dims::temperature >> native >> mks;
				const tw::Float Lmks = tw::small_pos + 8.0e-14 * sqr(TKelvin); // mean free path in meters
				const tw::Float Imks = 4.0*stef_boltz*std::pow(TKelvin,4);
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
		const tw::Float Ax0 = space->dS(i,j,k,1);
		const tw::Float Ax1 = space->dS(i+1,j,k,1);
		const tw::Float Ay0 = space->dS(i,j,k,2);
		const tw::Float Ay1 = space->dS(i,j+1,k,2);
		const tw::Float Az0 = space->dS(i,j,k,3);
		const tw::Float Az1 = space->dS(i,j,k+1,3);

		const tw::Float Pc = eos1(1,i,j,k,g->eidx.P);

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
			for (auto chem : g->chemical) {
				chem->PumpFluid(creationRate,destructionRate,g->hidx);
			}
			if (g->chemical[0]==electrons || g->forceFilter==0.0)
			{
				// electrons are not advanced using hydro
				// charge and current are computed implicitly in field advance
				// momentum density is fixed to give the heavy particle velocity (see ApplyEOS)
				// the velocity is then used to convect the energy density
				// mass and momentum density are convected, but then reset to restore quasineutrality and heavy particle velocity

				// take away any collisional sources of momentum
				for (auto cell : InteriorCellRange(*this,1))
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
					for (auto cell : InteriorCellRange(*this,1))
					{
						tw::Float dV,dS0,dS1,dl0,dl1,P0,P1,v0,v1;
						tw::Float forceDensity = 0.0;
						tw::Float powerDensity = 0.0;
						space->GetCellMetrics(cell,ax,&dV,&dS0,&dS1,&dl0,&dl1);

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
				if (space->gridGeometry==tw::grid::cylindrical) {
					for (auto cell : InteriorCellRange(*this,1)) {
						const tw::Float nm = g->DensityWeightedSum(state1,g->matset.mass,cell);
						const tw::Float Pc = eos1(cell,g->eidx.P);
						const tw::vec3 vc = g->Velocity(state1,cell);
						const tw::vec3 pos = space->Pos(cell);
						CreateMomentum(cell,1,fluxMask(cell)*(nm*sqr(vc.y) + Pc)/pos.x,g->hidx);
						DestroyMomentum(cell,2,fluxMask(cell)*nm*vc.x*vc.y/pos.x,g->hidx);
					}
				}
				if (space->gridGeometry==tw::grid::spherical) {
					for (auto cell : InteriorCellRange(*this,1)) {
						const tw::Float nm = g->DensityWeightedSum(state1,g->matset.mass,cell);
						const tw::Float Pc = eos1(cell,g->eidx.P);
						const tw::vec3 vc = g->Velocity(state1,cell);
						const tw::vec3 pos = space->Pos(cell);
						const tw::Float tanz = std::tan(pos.z);
						CreateMomentum(cell,1,fluxMask(cell)*(nm*(sqr(vc.y) + sqr(vc.z)) + 2.0*Pc)/pos.x,g->hidx);
						DestroyMomentum(cell,2,fluxMask(cell)*(nm*(vc.x*vc.z - vc.y*vc.y/tanz) - Pc/tanz)/pos.x,g->hidx);
						DestroyMomentum(cell,3,fluxMask(cell)*nm*(vc.x*vc.y + vc.y*vc.z/tanz)/pos.x,g->hidx);
					}
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

	creationRate.CopyFromNeighbors(All(creationRate));
	creationRate.ApplyBoundaryCondition(All(creationRate));

	destructionRate.CopyFromNeighbors(All(destructionRate));
	destructionRate.ApplyBoundaryCondition(All(destructionRate));
}

void sparc::HydroManager::LaserAdvance(tw::Float dt)
{
	logger::DEBUG("laser advance");
	// Currently SPARC ignores differences in frequency between injected pulses.
	// The average frequency of all the pulses is used throughout (see also HydroManager::Initialize).

	auto vacuumE = [&] (tw::Complex A,tw::Float k,tw::Float z)
	{
		// take enveloped vector potential and get space resolved electric field
		// only works in vacuum (e.g., at boundaries)
		return std::exp(ii*k*z)*ii*std::fabs(k)*A;
	};

	if (waves.size() && lasModel==sparc::vacuum) // use a prescribed field
	{
		#pragma omp parallel
		{
			for (auto cell : InteriorCellRange(*this,1))
			{
				tw::Complex fwd = 0.0;
				tw::Complex bak = 0.0;
				const tw::vec3 pos = space->Pos(cell);
				for (auto pulse : waves)
				{
					if (pulse->direction.z > 0.0)
						fwd += pulse->VectorPotentialEnvelope(space->WindowPos(0),pos,laserFrequency);
					else
						bak += pulse->VectorPotentialEnvelope(space->WindowPos(0),pos,laserFrequency);
				}
				laserAmplitude.Pack(cell, vacuumE(fwd,laserFrequency,pos.z) + vacuumE(bak,-laserFrequency,pos.z));
				radiationIntensity(cell) = 0.5*norm(laserAmplitude(cell));
			}
		}
	}

	if (electrons && waves.size() && lasModel==sparc::isotropic && dim[3]>1)
	{
		#pragma omp parallel
		{
			tw::Complex a0,a1; // amplitudes of incoming waves on left
			tw::Complex aN,aN1; // amplitudes of incoming waves on right
			for (auto strip : StripRange(*this,3,0,1,strongbool::no))
			{
				a0 = a1 = aN = aN1 = 0.0;
				const tw::Float z0 = space->Pos(strip,0).z;
				const tw::Float z1 = space->Pos(strip,1).z;
				const tw::Float zN = space->Pos(strip,dim[3]).z;
				const tw::Float zN1 = space->Pos(strip,dim[3]+1).z;
				for (auto pulse : waves)
				{
					if (pulse->direction.z > 0.0)
					{
						a0 += pulse->VectorPotentialEnvelope(space->WindowPos(0),z0,laserFrequency);
						a1 += pulse->VectorPotentialEnvelope(space->WindowPos(0),z1,laserFrequency);
					}
					else
					{
						aN += pulse->VectorPotentialEnvelope(space->WindowPos(0),zN,laserFrequency);
						aN1 += pulse->VectorPotentialEnvelope(space->WindowPos(0),zN1,laserFrequency);
					}
				}
				a0 = vacuumE(a0,laserFrequency,z0); a1 = vacuumE(a1,laserFrequency,z1);
				aN = vacuumE(aN,-laserFrequency,zN); aN1 = vacuumE(aN1,-laserFrequency,zN1);
				laserPropagator->SetupIncomingWaveLeft(strip,laserAmplitude,a0,a1,laserFrequency);
				laserPropagator->SetupIncomingWaveRight(strip,laserAmplitude,aN,aN1,laserFrequency);

				for (tw::Int k=1;k<=dim[3];k++)
				{
					refractiveIndex.Pack(strip,k, 0.0);
					// Add dispersionless part of susceptibility
					for (auto grp : group)
						for (auto chem : grp->chemical)
						{
							tw::Complex susceptibility(chem->mat.eps[0] - 1.0,chem->mat.eps[1]);
							refractiveIndex(strip,k,0) += state1(strip,k,chem->indexInState) * susceptibility.real();
							refractiveIndex(strip,k,1) += state1(strip,k,chem->indexInState) * susceptibility.imag();
						}
					// Add plasma contribution to susceptibility
					const tw::Float plas_real = -state1(strip,k,ie)/(sqr(laserFrequency) + sqr(nu_e(strip,k)));
					refractiveIndex(strip,k,0) += plas_real;
					refractiveIndex(strip,k,1) -= plas_real*nu_e(strip,k)/laserFrequency;
					// Convert susceptibility to refractive index
					refractiveIndex.Pack(strip,k, std::sqrt(one + refractiveIndex(strip,k)));
				}
			}
		}

		laserPropagator->Advance(laserAmplitude,refractiveIndex,nu_e,laserFrequency,dt);

		#pragma omp parallel
		{
			for (auto cell : InteriorCellRange(*this,1))
				radiationIntensity(cell) = real(refractiveIndex(cell))*0.5*norm(laserAmplitude(cell));
		}
	}
}

tw::Float sparc::HydroManager::EstimateTimeStep()
{
	logger::DEBUG("estimate time step");
	std::vector<tw::Float> dtMax(tw::GetOMPMaxThreads());
	tw::Float dtMaxAllThreads,dtMaxAllNodes;

	#pragma omp parallel
	{
		const tw::Int tid = tw::GetOMPThreadNum();
		const tw::Float sqrt_eps = std::sqrt(epsilonFactor);
		const tw::Float creationDominance = 100.0;
		tw::Float dts;

		if (space->IsFirstStep())
			dtMax[tid] = space->dx(0);
		else
			dtMax[tid] = space->MaxSpacing(0);

		// Courant condition

		for (auto cell : InteriorCellRange(*this,1))
			for (tw::Int s=0;s<group.size();s++)
			{
				const tw::vec3 vel = group[s]->Velocity(state1,cell);
				if (Norm(vel)!=0.0)
				{
					const tw::Float dxi2 = sqr(0.5*space->dL(cell,1));
					const tw::Float deta2 = sqr(0.5*space->dL(cell,2));
					const tw::Float dzeta2 = sqr(0.5*space->dL(cell,3));
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
				for (auto cell : InteriorCellRange(*this,1))
					dtMax[tid] = AsymptoticStep(cell,chem->indexInState);

			for (auto cell : InteriorCellRange(*this,1))
				dtMax[tid] = AsymptoticStep(cell,g->hidx.u);

			for (auto cell : InteriorCellRange(*this,1))
				dtMax[tid] = AsymptoticStep(cell,g->hidx.x);
		}
	} // end parallel region

	// Choose the smallest maximum step from all the threads
	dtMaxAllThreads = *std::min_element(std::begin(dtMax),std::end(dtMax));
	// Choose the smallest maximum step from all the nodes
	dtMaxAllNodes = task->strip[0].GetMin(dtMaxAllThreads);

	// Ensure that step size is sufficently small to resolve each laser pulse.
	tw::Float dtMaxAllPulses = space->MaxSpacing(0);
	for (auto pulse : waves)
	{
		const PulseShape& shape = pulse->pulseShape;
		const tw::Float pulseDuration = shape.t4 - shape.t1;
		// For upcoming pulses...
		if ( space->WindowPos(0) < shape.delay )
		{
			// Find closest pulse
			if ( dtMaxAllPulses > shape.delay - space->WindowPos(0) )
				dtMaxAllPulses = shape.delay - space->WindowPos(0);
		}
		// While inside a pulse
		if (space->WindowPos(0)>=shape.delay && space->WindowPos(0)<shape.delay+pulseDuration)
		{
			// Ensure Step resolves shortest pulse
			if ( dtMaxAllPulses > std::sqrt(epsilonFactor)*pulseDuration )
				dtMaxAllPulses = std::sqrt(epsilonFactor)*pulseDuration;
		}
		if ( dtMaxAllNodes > dtMaxAllPulses )
			dtMaxAllNodes = dtMaxAllPulses;
	}

	// If the critical time step is reached, switch to fixed time step and hope for the best!
	if (dtMaxAllNodes < space->CriticalSpacing(0))
	{
		dtMaxAllNodes = space->dx(0); // use the starting time step
		space->adaptiveTimestep = false; // no more adaptive stepping after this
	}
	// Don't let it fall below the minimum time step
	return dtMaxAllNodes < space->MinSpacing(0) ? space->MinSpacing(0) : dtMaxAllNodes;
}

void sparc::HydroManager::DiffusionAdvance(tw::Float dt)
{
	logger::DEBUG("diffusion advance");
	for (auto g : group)
	{
		// HEAT CONDUCTION

		for (auto c : conductors)
			parabolicSolver->FixTemperature(eos1,Rng(g->eidx.T),c->theRgn,c->Temperature(space->WindowPos(0)));
		g->LoadMassDensity(scratch,state1);
		CopyFieldData(scratch2,Rng(0),eos1,Rng(g->eidx.T));
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
				for (auto cell : InteriorCellRange(*this,1))
					state1(cell,g->hidx.npx+ax-1) = scratch(cell)*scratch2(cell);
			}
		}
	}

	state1.CopyFromNeighbors(All(state1));
	state1.ApplyBoundaryCondition(All(state1));
	eos1.CopyFromNeighbors(All(eos1));
	eos1.ApplyBoundaryCondition(All(eos1));
}

void sparc::HydroManager::FieldAdvance(tw::Float dt)
{
	logger::DEBUG("field advance");
	if (!electrons)
		return;

	const tw::Float q0 = electrons->mat.charge;
	const tw::Float m0 = electrons->mat.mass;
	const tw::Int Pe = ((EquilibriumGroup*)electrons->super)->eidx.P;

	// set up source and coefficient
	if (plasModel==sparc::quasineutral)
	{
		#pragma omp parallel
		{
			tw::Float D1,D2,P0,P1,mu0,mu1,dV,dS0,dS1,dl0,dl1;
			for (auto cell : InteriorCellRange(*this,1))
			{
				rho(cell) = rho0(cell);
				for (tw::Int ax=1;ax<=3;ax++)
				{
					space->GetCellMetrics(cell,ax,&dV,&dS0,&dS1,&dl0,&dl1);
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
	for (auto c : conductors)
		ellipticSolver->FixPotential(phi,c->theRgn,c->Voltage(space->WindowPos(0)));
	ellipticSolver->SetCoefficients(&scratch);
	ellipticSolver->Solve(phi,rho,-1.0);

	// Compute the new charge density
	// TODO: #pragma omp parallel for collapse(3) schedule(static)
	for (tw::Int i=1;i<=dim[1];i++)
		for (tw::Int j=1;j<=dim[2];j++)
			for (tw::Int k=1;k<=dim[3];k++)
			{
				rho(i,j,k) = (phi(i,j,k) - phi(i-1,j,k))*space->dS(i,j,k,1) / space->dl(i,j,k,1);
				rho(i,j,k) += (phi(i,j,k) - phi(i+1,j,k))*space->dS(i+1,j,k,1) / space->dl(i+1,j,k,1);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j-1,k))*space->dS(i,j,k,2) / space->dl(i,j,k,2);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j+1,k))*space->dS(i,j+1,k,2) / space->dl(i,j+1,k,2);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j,k-1))*space->dS(i,j,k,3) / space->dl(i,j,k,3);
				rho(i,j,k) += (phi(i,j,k) - phi(i,j,k+1))*space->dS(i,j,k+1,3) / space->dl(i,j,k+1,3);
				rho(i,j,k) /= space->dS(i,j,k,0);
			}
	rho.CopyFromNeighbors();
	rho.ApplyBoundaryCondition();
}

void sparc::HydroManager::HydroAdvance(const tw::grid::axis& axis,tw::Float dt)
{
	tw::Int ax = tw::grid::naxis(axis);
	logger::DEBUG(std::format("convection advance axis {}",ax));

	if (dim[ax] > 1)
	{
		FCT_Driver convector(&state0,&state1,&scratch,&fluxMask,space);
		tw::bc::fld bc0 = (space->bc0[ax] == par::reflecting || space->bc0[ax] == par::axisymmetric) ? fld::dirichletCell : fld::neumannWall;
		tw::bc::fld bc1 = space->bc1[ax] == par::reflecting ? fld::dirichletCell : fld::neumannWall;

		for (auto g : group)
		{
			logger::TRACE(std::format("group {}",g->name));
			convector.SetDensityElements(Rng(g->hidx.first,g->hidx.last+1));
			logger::TRACE("load velocity");
			convector.SetVelocityElement(0);
			#pragma omp parallel
			{
				g->LoadVelocity(scratch,state1,ax);
				for (auto cell : EntireCellRange(*this,1))
					scratch(cell) *= fluxMask(cell);
			}
			logger::TRACE("convect");
			convector.Convect(axis,bc0,bc1,dt);
		}
		state0.ApplyBoundaryCondition(All(state0));
	}
}

void sparc::HydroManager::ChemAdvance(tw::Float dt)
{
	logger::DEBUG("chemistry advance");
	creationRate *= dt;
	state0 += creationRate;

	destructionRate *= dt;
	state0 -= destructionRate;

	state0.CopyFromNeighbors(All(state0));
	state0.ApplyBoundaryCondition(All(state0));
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
		sparc::hydro_set hidx;
		if (electrons) {
			hidx = ((EquilibriumGroup*)electrons->super)->hidx;
		}

		for (auto cell : EntireCellRange(*this,1))
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
				state1(cell,hidx.npx) = ne*electrons->mat.mass*ionVelocity.x/tw::Float(group.size()-1);
				state1(cell,hidx.npy) = ne*electrons->mat.mass*ionVelocity.y/tw::Float(group.size()-1);
				state1(cell,hidx.npz) = ne*electrons->mat.mass*ionVelocity.z/tw::Float(group.size()-1);
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
	for (auto cell : EntireCellRange(*this,1))
		for (tw::Int c=0;c<state1.Components();c++)
			badCells += std::isnan(state1(cell,c));
	if (badCells)
		throw tw::FatalError("Encountered NaN in hydrodynamic state");
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
	logger::DEBUG("main update");

	if (!space->adaptiveTimestep)
	{
		FirstOrderAdvance(0.5*dx(0),true);
		FirstOrderAdvance(dx(0),true);
	}
	else
	{
		ComputeElectronCollisionFrequency();
		ComputeSources();
		space->ChangeStepSize(EstimateTimeStep());
		FirstOrderAdvance(0.5*space->dX(1,0),false);
		FirstOrderAdvance(space->dX(1,0),true);
	}

	LaserAdvance(space->dX(1,0));
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
	Driver::ReadCheckpoint(inFile);
	eos1.ReadCheckpoint(inFile);
	state1.ReadCheckpoint(inFile);
	rho.ReadCheckpoint(inFile);
	phi.ReadCheckpoint(inFile);
}

void sparc::HydroManager::WriteCheckpoint(std::ofstream& outFile)
{
	Driver::WriteCheckpoint(outFile);
	eos1.WriteCheckpoint(outFile);
	state1.WriteCheckpoint(outFile);
	rho.WriteCheckpoint(outFile);
	phi.WriteCheckpoint(outFile);
}

void sparc::HydroManager::Report(Diagnostic& diagnostic)
{
	Driver::Report(diagnostic);

	std::map<tw::Int,std::string> xyz = {{1,"x"},{2,"y"},{3,"z"}};

	auto WriteSubmoduleData = [&] (Field& theData,tw::Int comp,const std::string& filename,const tw::dims units=tw::dims::none,const std::string& pretty="tw::none")
	{
		CopyFieldData(scratch,Rng(0),theData,Rng(comp));
		scratch *= fluxMask;
		diagnostic.ReportField(filename,scratch,1,0,units,pretty);
	};

	// Mass, Charge, Energy

	scratch = 0.0;
	for (auto g : group)
	{
		g->LoadMassDensity(scratch2,state1);
		scratch2 *= fluxMask;
		scratch += scratch2;
	}
	diagnostic.VolumeIntegral("Mass",scratch,1,0);
	diagnostic.ReportField("massdensity",scratch,1,0,tw::dims::mass_density,"$\\rho$");

	diagnostic.VolumeIntegral("Charge",rho,1,0);

	scratch = 0.0;
	for (auto g : group)
		for (auto cell : InteriorCellRange(*this,1))
			scratch(cell) += state1(cell,g->hidx.u);
	diagnostic.VolumeIntegral("Energy",scratch,1,0);

	// Momentum

	for (tw::Int ax=1;ax<=3;ax++)
	{
		scratch = 0.0;
		for (auto g : group)
			for (auto cell : InteriorCellRange(*this,1))
				scratch(cell) += state1(cell,g->hidx.npx+ax-1);
		diagnostic.VolumeIntegral("P"+xyz[ax],scratch,1,0);
	}

	// Dipole moment

	diagnostic.FirstMoment("Dx",rho,1,0,dipoleCenter,tw::grid::x);
	diagnostic.FirstMoment("Dy",rho,1,0,dipoleCenter,tw::grid::y);
	diagnostic.FirstMoment("Dz",rho,1,0,dipoleCenter,tw::grid::z);

	// Impermeable Region

	diagnostic.ReportField("impermeable",fluxMask,1,0);

	// Collision Diagnostic

	diagnostic.ReportField("collisionFreq",nu_e,1,0,tw::dims::frequency,"$\\nu_e$");

	// Radiation Diagnostic

	diagnostic.ReportField("rad-intensity",radiationIntensity,1,0,tw::dims::intensity,"$I$");
	diagnostic.ReportField("rad-ereal",laserAmplitude,1,0,tw::dims::electric_field,"$\\Re E$");
	diagnostic.ReportField("rad-eimag",laserAmplitude,1,1,tw::dims::electric_field,"$\\Im E$");
	diagnostic.ReportField("rad-losses",radiativeLosses,1,0,tw::dims::power_density,"${\\cal P}$");
	diagnostic.ReportField("rad-nreal",refractiveIndex,1,0);
	diagnostic.ReportField("rad-nimag",refractiveIndex,1,1);
	diagnostic.ReportField("me_eff",me_eff,1,0);

	// Electrostatics

	if (electrons)
	{
		diagnostic.ReportField("chem-phi",phi,1,0,tw::dims::scalar_potential,"$\\phi$");
		diagnostic.ReportField("chem-rho",rho,1,0,tw::dims::charge_density,"$\\rho$");
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
			diagnostic.ReportField("v"+xyz[ax]+"_"+g->name,scratch,1,0,tw::dims::velocity,"$v_"+xyz[ax]+"$");
		}
	}
}

void sparc::HydroManager::StatusMessage(std::ostream *dest)
{
	std::println(*dest,"{}",statusMessage.str());
}
