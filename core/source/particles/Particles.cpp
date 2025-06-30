module;

#include <algorithm>
#include <tree_sitter/api.h>
#include "tw_includes.h"

export module particles;
import input;
import compute_tool;
import twmodule;
import fields;
import parabolic;
import mover;
import physics;
import qed;
import region;
import injection;
import diagnostics;
import laser_solve;

using namespace tw::bc;

// Particle struct is defined in discreteSpace.h
// Heavy lifting is done by classes in particles_*.h

struct LoadingData
{
	tw::Int timeLevel;
	tw::cell cell;
	tw::Float C0,C1,C2,densToAdd,densNow,particleDensity;
	tw::vec3 thermalMomentum,driftMomentum;
	std::valarray<tw::vec3> subGrid;
	tw::Int pointsInSubGrid;
	bool neutralize;

	LoadingData(const Simulation& sim,const tw::vec3& distributionInCell,const tw::cell& c);
	tw::Float GeometryFactor(tw::Float r,tw::Float r0) const
	{
		return fabs(C0 + SafeDiv(C1*r,r0) + sqr(SafeDiv(sqrt(C2)*r,r0)));
	}
};

export struct Species:Module
{
	tw::Float restMass,charge;
	tw::vec3 emissionTemp;
	tw::vec3 distributionInCell;
	tw::Float targetDensity;
	tw::Float minimumDensity;
	tw::Float accelerationTime,accelerationImpulse,accelerationForceNow;
	tw::Float numberCreated;
	tw::bc::par bc0[4],bc1[4];
	bool mobile,radiationDamping;
	tw::Float meanFreePath;
	tw::Int count,sortPeriod;

	std::vector<Particle> particle;
	std::vector<TransferParticle> transfer;

	Mover* mover;
	Ionizer* ionizer;
	QED* qed;

	Field* EM; // Ex,Ey,Ez,Bx,By,Bz
	Field* sources; // rho,Jx,Jy,Jz
	Field* laser; // F0x,F0y,F0z,F1x,F1y,F1z,aa0,aa1
	ComplexField* chi;
	tw::Float* carrierFrequency;
	tw_polarization_type* polarizationType;
	ScalarField* rho00;

	Field *qo_j4; // 4-current from quantum optics modules

	Species(const std::string& name,Simulation *sim);
	virtual ~Species();
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void VerifyInput();
	virtual void Initialize();
	void AddParticle(const float& number,const Primitive& q,const tw::vec4& p,const tw::vec4& s,const tw::Float& Qparam);
	void AddParticle(const TransferParticle& newParticle);
	void CleanParticleList();

	void ApplyGlobalBoundaryConditions();

	void ComputeTransferParticleDestinations();
	void PrepareTransfer(std::vector<TransferParticle>& accumulator,std::vector<tw::Int>& tally,tw::Int axis,tw::Int displ);
	void FinishTransfer(TransferParticle* inBuffer,tw::Int number,tw::Int axis,tw::Int displ);
	void CollectTransfers();

	void BeginMoveWindow();
	void FinishMoveWindow();

	void GenerateParticles(bool init);
	tw::Float AddDensity(const LoadingData& theData);
	tw::Float AddDensityRandom(const LoadingData& theData);
	void DepositInitialCharge(const tw::vec4& pos,tw::Float macroCharge);
	void CalculateDensity(ScalarField& dens);
	void CalculateEnergyDensity(ScalarField& dens);
	void CalculateQuantumParameter(ScalarField& dens);
	tw::Float KineticEnergy(const Region& theRgn);
	tw::Float ParticleNumber(const Region& theRgn);

	virtual bool ReadInputFileDirective(const TSTreeCursor *curs,const std::string& src);
	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
	virtual void Report(Diagnostic&);
	virtual void WarningMessage();

	virtual bool Test(tw::Int& id);
	void EncodingTest();
	void ReflectionTest();
	void MoveWindowTest();
};

export struct Kinetics:Module
{
	std::vector<Species*> species; // explicitly typed submodules

	ScalarField rho00;
	Field* sources;
	ComplexField* chi;
	ScalarField* ESRho;

	Kinetics(const std::string& name,Simulation *sim);
	virtual void Initialize();
	virtual void ExchangeResources();
	virtual bool InspectResource(void* resource,const std::string& description);
	virtual void Update();
	virtual void MoveWindow();
	void Ionize();
	void ProcessQED();

	void TransferParticles();
	tw::Float KineticEnergy(const Region& theRgn);
	virtual void Report(Diagnostic&);

	virtual void ReadCheckpoint(std::ifstream& inFile);
	virtual void WriteCheckpoint(std::ofstream& outFile);
};

////////////////////
// KINETICS CLASS //
////////////////////
// This is the top of the particle containment hierarchy
// Kinetics <- Species <- Particle

Kinetics::Kinetics(const std::string& name,Simulation *sim) : Module(name,sim)
{
	if (native.native!=tw::units::plasma && native.native!=tw::units::atomic)
		throw tw::FatalError("Kinetics module requires <native units = plasma> or <native units = atomic>.");

	rho00.Initialize(*this,owner);
	sources = NULL;
	chi = NULL;
}

void Kinetics::Initialize()
{
	Module::Initialize();
	for (auto sub : submodule)
	{
		Species *sp = dynamic_cast<Species*>(sub);
		if (sp!=NULL)
			species.push_back(sp);
	}
}

void Kinetics::ExchangeResources()
{
	PublishResource(&rho00,"kinetics:rho00");
}

bool Kinetics::InspectResource(void* resource,const std::string& description)
{
	if (description=="electromagnetic:sources")
	{
		sources = (Field*)resource;
		return true;
	}

	if (description=="laser:chi")
	{
		chi = (ComplexField*)resource;
		return true;
	}

	return false;
}

void Kinetics::MoveWindow()
{
	tw::Int i;
	Module::MoveWindow();

	// rho00 must be prepared to receive deposition from incoming particles
	// rho00 array is always left in an unrefined post-deposition state
	// i.e., charge from a given particle should be on only one node
	// therefore charge has to be zeroed on one node after being sent to another
	rho00.DownwardDeposit(tw::grid::z,1);
	for (auto s : StripRange(*this,3,strongbool::yes))
	{
		rho00(s,0) = 0.0;
		rho00(s,1) = 0.0;
		rho00.Shift(s,-1,0.0);
	}

	for (i=0;i<species.size();i++)
		species[i]->BeginMoveWindow();
	TransferParticles();
	for (i=0;i<species.size();i++)
		species[i]->FinishMoveWindow();
}

void Kinetics::Update()
{
	tw::Int i;

	// Assume source arrays are weighted by volume or zero

	Ionize();
	ProcessQED();
	for (i=0;i<species.size();i++)
	{
		if (owner->StepNow() % species[i]->sortPeriod == 0)
			std::sort(species[i]->particle.begin(),species[i]->particle.end());
		species[i]->mover->Advance();
		species[i]->ApplyGlobalBoundaryConditions();
	}
	// This barrier helps keep the stack trace clean for debugging purposes.
	// N.b. barriers are not implemented in TW_MPI.
	owner->strip[0].Barrier();

	TransferParticles();

	for (i=0;i<species.size();i++)
	{
		species[i]->GenerateParticles(false);
		species[i]->CleanParticleList();
	}

	if (sources && owner->neutralize)
		for (auto cell : EntireCellRange(*this))
			(*sources)(cell,0) -= rho00(cell);

	// Take care of boundaries in source arrays within field solvers
	// otherwise every module that deposits something would have to do message passing
}

void Kinetics::TransferParticles()
{
	tw::Int dst,src;
	tw::Int numToSend,sendSize,recvSize;
	tw::Int offset;
	bool odd;

	std::valarray<char> inBuffer,outBuffer;
	TransferParticle *parPtr;

	std::vector<TransferParticle> accumulator;
	std::vector<tw::Int> tally;

	for (auto sp : species)
		sp->ComputeTransferParticleDestinations();

	for (auto a=1;a<=3;a++)
	{
		if (owner->Dim(a)>1)
		{
			for (auto displ=-1;displ<=1;displ+=2)
			{
				accumulator.clear();
				tally.clear();

				owner->strip[a].Shift(1,displ,&src,&dst);
				odd = owner->strip[a].Get_rank() % 2;

				// Load particles and tallies into standard containers for all species
				// These are the transfer particles that need to move along the given axis

				for (auto sp: species)
					sp->PrepareTransfer(accumulator,tally,a,displ);

				// Exchange message sizes

				numToSend = accumulator.size();
				sendSize = sizeof(tw::Int)*species.size() + sizeof(TransferParticle)*numToSend;
				recvSize = sizeof(tw::Int)*species.size();
				if (odd)
				{
					owner->strip[a].Recv(&recvSize,sizeof(tw::Int),src);
					owner->strip[a].Send(&sendSize,sizeof(tw::Int),dst);
				}
				else
				{
					owner->strip[a].Send(&sendSize,sizeof(tw::Int),dst);
					owner->strip[a].Recv(&recvSize,sizeof(tw::Int),src);
				}

				// Pack data into the output buffer
				inBuffer.resize(recvSize);
				outBuffer.resize(sendSize);
				for (auto i=0;i<species.size();i++)
					((tw::Int*)&outBuffer[0])[i] = tally[i];
				if (numToSend)
				{
					parPtr = (TransferParticle*)(&outBuffer[species.size()*sizeof(tw::Int)]);
					for (auto i=0;i<numToSend;i++)
						parPtr[i] = accumulator[i];
				}

				// Exchange particles

				if (odd)
				{
					owner->strip[a].Recv(&inBuffer[0],recvSize,src);
					owner->strip[a].Send(&outBuffer[0],sendSize,dst);
				}
				else
				{
					owner->strip[a].Send(&outBuffer[0],sendSize,dst);
					owner->strip[a].Recv(&inBuffer[0],recvSize,src);
				}

				// Add particles to each species transfer list

				if (recvSize>sizeof(tw::Int)*species.size())
				{
					offset = 0;
					parPtr = (TransferParticle*)(&inBuffer[species.size()*sizeof(tw::Int)]);
					for (auto i=0;i<species.size();i++)
					{
						tally[i] = ((tw::Int*)&inBuffer[0])[i];
						if (tally[i])
							species[i]->FinishTransfer(&parPtr[offset],tally[i],a,displ);
						offset += tally[i];
					}
				}
			}
		}
	}

	for (auto sp: species)
		sp->CollectTransfers();
}

void Kinetics::Ionize()
{
	// There is a diplacement of the bound particle before it turns free, satisfying qE.dr = Uion.
	// This accounts for the work done to remove the particle from the binding potential.
	// Ionization due to wake and laser are simply added (OK if they do not overlap much)
	// Time centering:
	// Fp[6] is a2 at t = 0, carrier resolved field also known at t = 0
	// Ionized particle is being advanced in momentum from t = -1/2 to t = 1/2, and space from t = 0 to t = 1

	weights_3D weights;
	tw::Float w0,a2,probability;
	std::valarray<tw::Float> temp(6),Fp(8);
	tw::vec3 E,vel,momentum;
	tw::Float gamma,m0,q0;
	Species *s1,*s2;

	// Find Laser Solver

	LaserSolver *theLaserSolver = NULL;
	for (tw::Int i=0;i<owner->module.size();i++)
		if (dynamic_cast<LaserSolver*>(owner->module[i]))
			theLaserSolver = (LaserSolver*)owner->module[i];

	// Photoionization

	for (tw::Int s=0;s<species.size();s++)
	{
		Ionizer *ionizer = species[s]->ionizer;

		if (ionizer!=NULL)
		{
			std::vector<Particle>& particle = species[s]->particle;
			s1 = (Species*)owner->GetModule(ionizer->ion_name);
			s2 = (Species*)owner->GetModule(ionizer->electron_name);
			m0 = species[s]->restMass;
			q0 = species[s]->charge;
			for (tw::Int i=0;i<particle.size();i++)
			{
				Particle& curr = particle[i]; // curr = particle being ionized
				owner->StaticSpace::GetWeights(&weights,curr.q);
				probability = a2 = 0.0;
				if (species[s]->laser)
				{
					species[s]->laser->Interpolate(Fp,Element(6),weights);
					a2 = Fp[6];
					w0 = *species[s]->carrierFrequency;
					if (theLaserSolver->polarizationType==circularPolarization)
						probability += this->dx(0)*ionizer->InstantRate(w0,w0*sqrt(0.5*a2));
					else
						probability += this->dx(0)*ionizer->AverageRate(w0,w0*sqrt(a2));
				}
				species[s]->EM->Interpolate(temp,Element(0,2),weights);
				E.x = temp[0]; E.y = temp[1]; E.z = temp[2];
				probability += this->dx(0)*ionizer->InstantRate(1e-6,Magnitude(E));
				gamma = sqrt(sqr(curr.p[0]/m0) + 0.5*sqr(q0/m0)*a2);
				// Starting velocity is that of the neutral
				vel = curr.p.spatial()/(gamma*m0);
				if (probability > owner->uniformDeviate->Next())
				{
					momentum = s1->restMass*gamma*vel;
					if (theLaserSolver)
						momentum += theLaserSolver->GetIonizationKick(a2,s1->charge,s1->restMass);
					s1->AddParticle(curr.number,curr.q,tw::vec4(gamma*s1->restMass,momentum),tw::vec4(0.0),0.0);

					// For the electron, account for depletion of the field due to ionization energy.
					// This comes from a displacement satisfying dr.QE = dU, where dU = ionization energy.
					// const tw::Float dU = ionizer->ionizationPotential;
					// const tw::Float failsafe = sqr(0.01*ionizer->ThresholdEstimate());
					// const tw::vec3 dr = E*dU/(s2->charge*curr.number*(Norm(E)+failsafe));

					momentum = s2->restMass*gamma*vel;
					if (theLaserSolver)
						momentum += theLaserSolver->GetIonizationKick(a2,s2->charge,s2->restMass);
					s2->AddParticle(curr.number,curr.q,tw::vec4(gamma*s2->restMass,momentum),tw::vec4(0.0),0.0);

					// Throw away the neutral
					particle[i] = particle.back();
					particle.pop_back();
					--i;
				}
			}
		}
	}
}

void Kinetics::ProcessQED()
{
	Species *s1, *s2;
	weights_3D weights;

	tw::vec3 E, B, vel;
	std::valarray<tw::Float> temp(6);

	tw::Float eta, chi, Py, Pf;
	tw::Float gamma, omega, energy;
	tw::Float frac, rate, probability;

	tw::vec3 k3u;
	tw::vec4 k4;

	tw::Float prefactor = (cgs::qe*cgs::hbar)/sqr(cgs::me)/sqr(sqr(cgs::c));

	for (tw::Int s=0;s<species.size();s++)
	{
		QED *qed = species[s]->qed;

		// Photon generation

		if (qed!=NULL && !qed->photon_name.empty())
		{
			std::vector<Particle>& fermion = species[s]->particle;
			s1 = (Species*)owner->GetModule(qed->photon_name);

			for (tw::Int i=0;i<fermion.size();i++)
			{
				// Gather momentum & local field components of current fermion
				Particle& part = fermion[i];

				gamma = sqrt(1.0 + Norm(part.p.spatial()));
				vel = part.p.spatial()/gamma/species[s]->restMass;

				owner->StaticSpace::GetWeights(&weights,part.q);
				species[s]->EM->Interpolate(temp,Element(0,5),weights);

				for (tw::Int n=0;n<3;n++)
					vel[n] = vel[n] * tw::dims::velocity >> native >> cgs;
				for (tw::Int n=0;n<3;n++)
					temp[n] = temp[n] * tw::dims::electric_field >> native >> cgs;
				for (tw::Int n=3;n<6;n++)
					temp[n] = temp[n] * tw::dims::magnetic_field >> native >> cgs;

				E.x = temp[0]; E.y = temp[1]; E.z = temp[2];
				B.x = temp[3]; B.y = temp[4]; B.z = temp[5];

				// Calculate fermion quantum parameter
				part.Qparam = prefactor*gamma*sqrt(Norm((cgs::c)*E + (vel|B)) - sqr(vel^E));

				// Calculate instantaneous photo-emission rate & probability
				rate = qed->CalculateRate(part.Qparam,gamma);
				probability = rate*this->dx(0);

				// Generate a photon if the emission criterion is satisfied
				if (probability > owner->uniformDeviate->Next())
				{
					do {
						Py = owner->uniformDeviate->Next();
						chi = qed->NewQParameter(part.Qparam,Py);
					} while (chi < 0.0f);

					k3u = part.p.spatial();
					k3u /= Magnitude(k3u) + tw::tiny;

					energy = 2.0*chi*cub(cgs::me)*std::pow(cgs::c,5)/(cgs::qe*cgs::hbar);
					energy /= sqrt(Norm(E + (k3u|B)) - sqr(k3u^E));
					energy = energy * tw::dims::energy >> cgs >> native;

					if (energy > qed->photon_cutoff_energy)
					{
						k4[0] = energy;
						for (int n=1;n<4;n++)
							k4[n] = energy*k3u[n-1];

						s1->AddParticle(part.number,part.q,k4,tw::vec4(0.0),chi);

						// Update fermion momentum (radiation damping)
						part.p -= k4;
						part.p[0] = sqrt(sqr(s1->restMass) + Norm(part.p.spatial()));

						// Update fermion quantum parameter
						gamma = sqrt(1.0 + Norm(part.p.spatial()));
						vel = part.p.spatial()/gamma/species[s]->restMass;
						for (tw::Int n=0;n<3;n++)
							vel[n] = vel[n] * tw::dims::velocity >> native >> cgs;
						part.Qparam = prefactor*gamma*sqrt(Norm((cgs::c)*E + (vel|B)) - sqr(vel^E));
					}
				}
			}
		}

		// Pair creation
		if (qed!=NULL && !qed->electron_name.empty() && !qed->positron_name.empty())
		{
			std::vector<Particle>& photon = species[s]->particle;
			s1 = (Species*)owner->GetModule(qed->electron_name);
			s2 = (Species*)owner->GetModule(qed->positron_name);

			for (tw::Int i=0;i<photon.size();i++)
			{
				// Gather wavevector and local field components of current photon
				Particle& part = photon[i];

				k3u = part.p.spatial();
				k3u /= Magnitude(k3u) + tw::tiny;

				omega = (part.p[0] * tw::dims::energy >> native >> cgs)/cgs::hbar;

				owner->StaticSpace::GetWeights(&weights,part.q);
				species[s]->EM->Interpolate(temp,Element(0,5),weights);

				for (tw::Int n=0;n<3;n++)
					temp[n] = temp[n] * tw::dims::electric_field >> native >> cgs;
				for (tw::Int n=3;n<6;n++)
					temp[n] = temp[n] * tw::dims::magnetic_field >> native >> cgs;

				E.x = temp[0]; E.y = temp[1]; E.z = temp[2];
				B.x = temp[3]; B.y = temp[4]; B.z = temp[5];

				// Calculate photon quantum parameter
				part.Qparam = prefactor*(cgs::hbar/cgs::me/2.0)*(omega/cgs::c);
				part.Qparam *= sqrt(Norm(E + (k3u|B)) - sqr(k3u^E));

				// Calculate instantaneous pair-creation rate & probability
				rate = qed->CalculateRate(part.Qparam,cgs::hbar*omega);
				probability = rate*this->dx(0);

				// Generate an electron-positron pair if the emission criterion is satisfied
				if (probability > owner->uniformDeviate->Next())
				{
					do {
						Pf = owner->uniformDeviate->Next();
						frac = qed->NewQParameter(part.Qparam,Pf);
					} while (frac < 0.0f);

					tw::vec4 p4a = frac*part.p;
					tw::vec4 p4b = (1.0f-frac)*part.p;

					// calculate fermion 1 quantum parameter
					gamma = sqrt(1.0 + Norm(p4a.spatial()));
					vel = p4a.spatial()/gamma;
					for (tw::Int n=0;n<3;n++)
						vel[n] = vel[n] * tw::dims::velocity >> native >> cgs;
					eta = prefactor*gamma*sqrt(Norm((cgs::c)*E + (vel|B)) - sqr(vel^E));
					// create fermion 1
					s1->AddParticle(part.number,part.q,p4a,tw::vec4(0.0),eta);

					// calculate fermion 2 quantum parameter
					gamma = sqrt(1.0 + Norm(p4b.spatial()));
					vel = p4b.spatial()/gamma;
					for (tw::Int n=0;n<3;n++)
						vel[n] = vel[n] * tw::dims::velocity >> native >> cgs;
					eta = prefactor*gamma*sqrt(Norm((cgs::c)*E + (vel|B)) - sqr(vel^E));
					// create fermion 2
					s2->AddParticle(part.number,part.q,p4b,tw::vec4(0.0),eta);

					// Mark the parent photon for annihilation
					photon[i].number = 0.0;
				}
			}
		}
	}
}

tw::Float Kinetics::KineticEnergy(const Region& theRgn)
{
	tw::Float ans = 0.0;
	for (auto sp : species)
	{
		const tw::Float m0 = sp->restMass;
		for (auto par : sp->particle)
			if (theRgn.Inside(owner->PositionFromPrimitive(par.q).spatial(),*owner))
				ans += par.number * (par.p[0] - m0);
	}
	return ans;
}

void Kinetics::Report(Diagnostic& diagnostic)
{
	Module::Report(diagnostic);
	diagnostic.ReportNumber("Kinetic",KineticEnergy(*diagnostic.theRgn),false);
}

void Kinetics::ReadCheckpoint(std::ifstream& inFile)
{
	Module::ReadCheckpoint(inFile);
	rho00.ReadCheckpoint(inFile);
}

void Kinetics::WriteCheckpoint(std::ofstream& outFile)
{
	Module::WriteCheckpoint(outFile);
	rho00.WriteCheckpoint(outFile);
}


///////////////////
// SPECIES CLASS //
///////////////////


Species::Species(const std::string& name,Simulation* sim) : Module(name,sim)
{
	if (native.native!=tw::units::plasma && native.native!=tw::units::atomic)
		throw tw::FatalError("Species module requires <native units = plasma> or <native units = atomic>.");

	restMass = 1.0;
	charge = -1.0;
	distributionInCell = tw::vec3(2.0,2.0,2.0);
	minimumDensity = 0.0;
	targetDensity = 1.0;
	accelerationTime = accelerationImpulse = accelerationForceNow = 0.0;
	numberCreated = 0.0;
	meanFreePath = 0.0; // means infinity
	count = 1;
	sortPeriod = 10;
	emissionTemp = tw::vec3(0.1,0.1,0.1);

	for (tw::Int i=1;i<=3;i++)
	{
		bc0[i] = owner->bc0[i];
		bc1[i] = owner->bc1[i];
	}

	mobile = true;
	radiationDamping = false;

	mover = NULL;
	ionizer = NULL;
	qed = NULL;

	EM = NULL;
	sources = NULL;
	laser = NULL;
	chi = NULL;
	carrierFrequency = NULL;
	polarizationType = NULL;
	rho00 = NULL;
	qo_j4 = NULL;

	directives.Add("xboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[1],&bc1[1]),false);
	directives.Add("yboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[2],&bc1[2]),false);
	directives.Add("zboundary",new tw::input::Enums<tw::bc::par>(tw::bc::par_map(),&bc0[3],&bc1[3]),false);
	directives.Add("sort period",new tw::input::Int(&sortPeriod),false);
	directives.Add("mobile",new tw::input::Bool(&mobile),false);
	directives.Add("radiation damping",new tw::input::Bool(&radiationDamping),false);
	directives.Add("mean free path",new tw::input::Float(&meanFreePath),false);
	directives.Add("mass",new tw::input::Float(&restMass),false);
	directives.Add("charge",new tw::input::Float(&charge),false);
	directives.Add("minimum density",new tw::input::Float(&minimumDensity),false);
	directives.Add("emission temperature",new tw::input::Vec3(&emissionTemp),false);
	directives.Add("acceleration time",new tw::input::Float(&accelerationTime),false);
	directives.Add("acceleration momentum",new tw::input::Float(&accelerationImpulse),false);
	directives.Add("particles per cell",new tw::input::Vec3(&distributionInCell),true);
	directives.Add("when density",new tw::input::Custom);
}

Species::~Species()
{
	if (mover!=NULL)
		owner->RemoveTool(mover);
	if (ionizer!=NULL)
		owner->RemoveTool(ionizer);
	if (qed!=NULL)
		owner->RemoveTool(qed);
}

void Species::VerifyInput()
{
	Module::VerifyInput();
	for (auto tool : moduleTool)
	{
		if (ionizer==NULL)
			ionizer = dynamic_cast<Ionizer*>(tool);
		if (mover==NULL)
			mover = dynamic_cast<Mover*>(tool);
		if (qed==NULL)
			qed = dynamic_cast<QED*>(tool);
	}
}

void Species::Initialize()
{
	Module::Initialize();
	GenerateParticles(true);
	CleanParticleList();
	if (!owner->neutralize)
		CopyFieldData(*sources,0,*rho00,0); // accepting redundant operations

	// Choose a mover if none was specified
	if (mover==NULL)
	{
		if (EM!=NULL && sources!=NULL && laser==NULL && qo_j4==NULL && restMass!=0.0f && charge!=0.0f)
			mover = (Mover*)owner->CreateTool("Boris-mover",tw::tool_type::borisMover);
		if (EM!=NULL && sources!=NULL && laser!=NULL && qo_j4==NULL)
			mover = (Mover*)owner->CreateTool("PGC-mover",tw::tool_type::pgcMover);
		if (qo_j4!=NULL)
			mover = (Mover*)owner->CreateTool("Bohmian-mover",tw::tool_type::bohmianMover);
		if (restMass==0.0f && charge==0.0f)
			mover = (Mover*)owner->CreateTool("Photon-mover",tw::tool_type::photonMover);
	}

	// Copy pointers to the mover tool
	mover->q0 = charge;
	mover->m0 = restMass;
	mover->particle = &particle;
	mover->transfer = &transfer;
	mover->EM = EM;
	mover->sources = sources;
	mover->laser = laser;
	mover->chi = chi;
	mover->qo_j4 = qo_j4;
}

void Species::WarningMessage()
{
	Module::WarningMessage();
	if (bc0[1]==par::absorbing || bc1[1]==par::absorbing ||
		bc0[2]==par::absorbing || bc1[2]==par::absorbing ||
		bc0[3]==par::absorbing || bc1[3]==par::absorbing)
		std::println(std::cout,"{}: {} will be absorbed",term::warning,name);
}

bool Species::InspectResource(void* resource,const std::string& description)
{
	if (description=="electromagnetic:F")
	{
		EM = (Field*)resource;
		return true;
	}

	if (description=="electromagnetic:sources")
	{
		sources = (Field*)resource;
		return true;
	}

	if (description=="kinetics:rho00")
	{
		rho00 = (ScalarField*)resource;
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

	if (description=="laser:polarizationType")
	{
		polarizationType = (tw_polarization_type*)resource;
		return true;
	}

	if (description=="qo:j4")
	{
		qo_j4 = (Field*)resource;
		return true;
	}

	return false;
}

void Species::AddParticle(const float& number,const Primitive& q,const tw::vec4& p,const tw::vec4& s,const tw::Float& Qparam)
{
	numberCreated += number;
	count++;
	particle.emplace_back(number,q,p,s,(count << 32) + owner->strip[0].Get_rank(),Qparam);
}

/// @brief Add the transfer particle to the final destination's particle list
/// @param xfer the transfer particle data
void Species::AddParticle(const TransferParticle& xfer)
{
	float x[4];
	tw::Int ijk[4];
	x[0] = xfer.x[0];
	ijk[0] = xfer.ijk[0];
	for (tw::Int ax=1;ax<4;ax++)
	{
		x[ax] = xfer.x[ax];
		ijk[ax] = xfer.ijk[ax] - dim[ax] * xfer.dst[ax];
		while (ijk[ax]<1)
		{
			ijk[ax]++;
			x[ax] -= 1.0;
		}
		while (ijk[ax]>dim[ax])
		{
			ijk[ax]--;
			x[ax] += 1.0;
		}
	}
	tw::Int cell = EncodeCell(ijk[0],ijk[1],ijk[2],ijk[3]);
	Primitive q(cell,x[0],x[1],x[2],x[3]);
	particle.emplace_back(xfer.number,q,xfer.p,xfer.s,xfer.tag,xfer.Qparam);
	// don't update count because the transfer particle already has its identifier in aux1 and aux2
}

void Species::CleanParticleList()
{
	for (tw::Int i=0;i<particle.size();i++)
	{
		if (particle[i].number==0.0)
		{
			// replace particle[i] with particle.back(), thus throwing out particle[i]
			particle[i] = particle.back();
			// throw out original particle.back() since it is now represented by particle[i]
			particle.pop_back();
			// decrement index so we can test what used to be particle.back()
			i--;
		}
	}
}



////////////////////////////////////
//   GLOBAL BOUNDARY CONDITIONS   //
////////////////////////////////////



/// Must happen before particle message passing and before particle destination calculation.
/// All particles we want to keep should be transformed so they are in the global domain,
/// except for particles crossing a periodic boundary.
/// Result of the transformation can result in movement to new local domain.
/// If absorption is desired leave particle out of global domain (destination calculation will send to MPI_PROC_NULL).
void Species::ApplyGlobalBoundaryConditions()
{
	tw::Int ax,i,displ,src,dst;
	tw::bc::par boundaryCondition;

	for (displ=-1;displ<=1;displ+=2)
		for (ax=1;ax<=3;ax++)
		{
			boundaryCondition = displ<0 ? bc0[ax] : bc1[ax];
			owner->strip[ax].Shift(1,displ,&src,&dst);
			if (dst==MPI_PROC_NULL && owner->Dim(ax)!=1)
			{
				for (i=0;i<transfer.size();i++)
				{
					if (transfer[i].dst[ax]==displ)
					{
						switch (boundaryCondition)
						{
							case par::periodic:
								// must delay particle wrap-around
								// typically this case will never be encountered since dst!=MPI_PROC_NULL
								break;
							case par::absorbing:
								// do nothing; let it leave
								break;
							case par::emitting:
								transfer[i].dst[ax] = 0;
								transfer[i].x[ax] = 0.0;
								transfer[i].ijk[ax] -= displ;
								transfer[i].p[1] = emissionTemp.x*owner->gaussianDeviate->Next();
								transfer[i].p[2] = emissionTemp.y*owner->gaussianDeviate->Next();
								transfer[i].p[3] = emissionTemp.z*owner->gaussianDeviate->Next();
								transfer[i].p[ax] = -tw::Float(displ)*fabs(transfer[i].p[ax]);
								// DFG: add a line to update the energy
								transfer[i].p[0] = sqrt(sqr(restMass) + Norm(transfer[i].p.spatial()));
								break;
							case par::reflecting:
							case par::axisymmetric:
								transfer[i].dst[ax] = 0;
								transfer[i].x[ax] *= -1.0;
								transfer[i].ijk[ax] -= displ;
								transfer[i].p[ax] *= -1.0;
								break;
							case par::none:
								// do nothing, same effect as absorbing
								break;
						}
					}
				}
			}
		}
}



///////////////////////////////////
//   MESSAGE PASSING PARTICLES   //
///////////////////////////////////



/// Given dst[0] = starting node, and dst[ax] = direction of transfer (-1,0,1), set dst[0] = final node.
/// The transfer can be arrested along any axis by setting dst[ax] = 0.
/// This must be called before message passing begins, and after global boundary conditions are imposed.
/// The destination domain for every particle on the transfer list has to be computed.
/// If the particle leaves the entire system the destination rank is set to MPI_PROC_NULL (only possible for absorbing boundaries).
void Species::ComputeTransferParticleDestinations()
{
	tw::Int i,ax,dest_coords[4];

	for (i=0;i<transfer.size();i++)
	{
		for (ax=1;ax<=3;ax++)
		{
			dest_coords[ax] = owner->domainIndex[ax] + transfer[i].dst[ax];
			// handle absorbed particles
			if (transfer[i].dst[ax]<0 && owner->n0[ax]==MPI_PROC_NULL)
				transfer[i].dst[0] = MPI_PROC_NULL;
			if (transfer[i].dst[ax]>0 && owner->n1[ax]==MPI_PROC_NULL)
				transfer[i].dst[0] = MPI_PROC_NULL;
			// handle periodicity
			if (owner->periodic[ax])
			{
				if (dest_coords[ax]<0)
					dest_coords[ax] = owner->domains[ax]-1;
				if (dest_coords[ax]>=owner->domains[ax])
					dest_coords[ax] = 0;
				if (owner->domains[ax]==1)
				{
					transfer[i].ijk[ax] -= dim[ax] * transfer[i].dst[ax];
					transfer[i].dst[ax] = 0;
				}
			}
		}
		if (transfer[i].dst[0]!=MPI_PROC_NULL)
			transfer[i].dst[0] = owner->strip[0].Cart_rank(dest_coords[1],dest_coords[2],dest_coords[3]);
	}
}

void Species::PrepareTransfer(std::vector<TransferParticle>& accumulator,std::vector<tw::Int>& tally,tw::Int ax,tw::Int displ)
{
	// displ is -1 for downward transfer, +1 for upward transfer
	tw::Int numToSend = 0;

	for (tw::Int i=0;i<transfer.size();i++)
	{
		if (transfer[i].dst[ax]==displ)
		{
			numToSend++;
			accumulator.push_back(transfer[i]);
		}
	}
	tally.push_back(numToSend);
}

void Species::FinishTransfer(TransferParticle* inBuffer,tw::Int numToReceive,tw::Int ax,tw::Int displ)
{
	for (tw::Int i=0;i<numToReceive;i++)
	{
		// put the result on the transfer list
		transfer.push_back(inBuffer[i]);
	}
}

void Species::CollectTransfers()
{
	// No transfer particle should ever be deleted except at the end of this routine.
	// The list of transfer particles includes any transfer particle that touched this node,
	// including particles that originated on this node.
	// This gives all nodes the opportunity to process the particle, if desired.
	for (tw::Int i=0;i<transfer.size();i++)
		if (transfer[i].dst[0]==owner->strip[0].Get_rank())
			AddParticle(transfer[i]);
	transfer.clear();
}


/////////////////////////
//  PARTICLE CREATION  //
/////////////////////////




void Species::DepositInitialCharge(const tw::vec4& pos,tw::Float macroCharge)
{
	if (rho00)
	{
		weights_3D weights;
		owner->GetWeights(&weights,pos);
		(*rho00).InterpolateOnto(macroCharge,weights);
	}
}

void Species::GenerateParticles(bool init)
{
	tw::Float add;
	for (auto prof : profile)
		if ( prof->TimeGate(owner->WindowPos(0),&add) )
			for (auto cell : InteriorCellRange(*this))
			{
				LoadingData loadingData(*owner,distributionInCell,cell);
				loadingData.densToAdd = prof->GetValue(owner->Pos(cell),*owner);
				loadingData.driftMomentum = prof->DriftMomentum(restMass);
				loadingData.neutralize = owner->neutralize;
				loadingData.thermalMomentum = prof->thermalMomentum;

				if (loadingData.densToAdd>0.0)
				{
					if (prof->timingMethod==tw::profile::timing::maintained && !init)
					{
						// Hold a constant charge
						// (comment out for constant current)
						if (owner->neutralize) // don't test loadingData.neutralize
							loadingData.densToAdd = -(*sources)(cell,0)/charge; // hold this region neutral
						else
							loadingData.densToAdd -= (*sources)(cell,0)/charge; // hold at target density (assumes no other species)
						loadingData.neutralize = false;
					}
					if (qo_j4!=NULL && prof->loadingMethod==tw::profile::loading::statistical && !prof->variableCharge)
					{
						loadingData.densToAdd = fabs((*qo_j4)(cell,0));
					}

					if (loadingData.densToAdd < 0.0) // constant charge block may lead to n<0
						loadingData.densToAdd = 0.0;

					// INJECT

					if (prof->variableCharge)
						loadingData.particleDensity = loadingData.densToAdd/(distributionInCell.x*distributionInCell.y*distributionInCell.z);
					else
						loadingData.particleDensity = targetDensity*prof->gammaBoost;
					if (prof->loadingMethod==tw::profile::loading::deterministic)
						AddDensity(loadingData);
					else
						AddDensityRandom(loadingData);
				}
			}
}

LoadingData::LoadingData(const Simulation& sim,const tw::vec3& distribution,const tw::cell& c) : cell(c)
{
	timeLevel = sim.StepNow();
	densToAdd = 0.0;
	densNow = 0.0;
	particleDensity = 0.0;
	neutralize = true;
	C0 = 1.0;
	C1 = 0.0;
	C2 = 0.0;

	tw::Int nx = tw::Int(distribution.x);
	tw::Int ny = tw::Int(distribution.y);
	tw::Int nz = tw::Int(distribution.z);

	pointsInSubGrid = nx*ny*nz;
	subGrid.resize(pointsInSubGrid);
	for (tw::Int k=0;k<nz;k++)
		for (tw::Int j=0;j<ny;j++)
			for (tw::Int i=0;i<nx;i++)
				subGrid[k*nx*ny + j*nx + i] = tw::vec3
					(
						-0.5 + 0.5/tw::Float(nx) + tw::Float(i)/tw::Float(nx),
						-0.5 + 0.5/tw::Float(ny) + tw::Float(j)/tw::Float(ny),
						-0.5 + 0.5/tw::Float(nz) + tw::Float(k)/tw::Float(nz)
					);

	const tw::Float Nr2 = nx*nx;
	if (sim.gridGeometry==tw::grid::cylindrical)
	{
		const tw::Int x = cell.dcd1();
		C0 = 0.0;
		C1 = 1.0;
		C2 = 0.0;
		if (sim.n0[1]==MPI_PROC_NULL)
		{
			if (Nr2==1.0)
			{
				C0 = 1.0;
				C1 = 0.0;
				if (x==1)
					C0 = 4208.0/5951.0;
				if (x==2)
					C0 = 18152.0/17853.0;
				if (x==3)
					C0 = 29704.0/29755.0;
				if (x==4)
					C0 = 5952.0/5951.0;

			}
			else
			{
				if (x==1)
				{
					C0 = 0.0;
					C1 = (7.0 - 22.0*Nr2)/(24.0*Nr2*(Nr2-1.0));
					C2 = (5.0/8.0) * (4.0*Nr2*Nr2-1.0) / ((Nr2-1.0)*(4.0*Nr2-1.0));
				}
			}
		}
	}
}

tw::Float Species::AddDensity(const LoadingData& theData)
{
	const tw::vec3 cellCenter = owner->Pos(theData.cell);
	//const tw::vec3 cellSize = owner->dPos(theData.cell);
	const tw::Float cellVolume = owner->dS(theData.cell,0);
	const tw::Float densToAdd = theData.densToAdd;
	const tw::Float densNow = theData.densNow;
	const tw::vec3 thermalMomentum = theData.thermalMomentum;
	const tw::vec3 driftMomentum = theData.driftMomentum;
	//const bool neutralize = theData.neutralize;
	const tw::Int pointsInSubGrid = theData.pointsInSubGrid;

	Primitive q;
	tw::vec4 p;
	tw::Float particleDensity;
	tw::Int numToAdd,numNow;
	std::valarray<tw::vec4> initialMomenta;

	particleDensity = theData.particleDensity;
	if (particleDensity <= minimumDensity)
		return 0.0;

	numNow = TruncateUnlessClose( densNow / particleDensity);
	numToAdd = TruncateUnlessClose( (densNow + densToAdd) / particleDensity) - numNow;

	if (numToAdd<1) return 0.0;

	p = tw::vec4(0.0);
	initialMomenta.resize(numToAdd);
	for (tw::Int i=0;i<numToAdd;i++)
	{
		initialMomenta[i][1] = thermalMomentum.x*owner->gaussianDeviate->Next();
		initialMomenta[i][2] = thermalMomentum.y*owner->gaussianDeviate->Next();
		initialMomenta[i][3] = thermalMomentum.z*owner->gaussianDeviate->Next();
		p += initialMomenta[i];
	}
	p /= tw::Float(numToAdd);
	initialMomenta -= p - tw::vec4(0,driftMomentum);

	for (tw::Int i=0;i<numToAdd;i++)
	{
		tw::Int j = numNow + i;
		while (j>=pointsInSubGrid) j -= pointsInSubGrid;
		p = initialMomenta[i];
		p[0] = sqrt(sqr(restMass) + Norm(p.spatial()));
		q.x[0] = 0.0; // always center of time cell
		q.x[1] = theData.subGrid[j].x;
		q.x[2] = theData.subGrid[j].y;
		q.x[3] = theData.subGrid[j].z;
		// At the beginning, timeLevel = 0.  So the particle's time will be t=-0.5*dt which can be viewed
		// as in the center of the temporal ghost cell.  The field data is known at t = 0, i.e., on the
		// high wall of the temporal ghost cell.
		q.cell = EncodeCell(theData.timeLevel,theData.cell.dcd1(),theData.cell.dcd2(),theData.cell.dcd3());
		const tw::vec4 x = PositionFromPrimitive(q);
		const tw::Float N = theData.GeometryFactor(x[1],cellCenter.x)*particleDensity*cellVolume;
		AddParticle(N,q,p,tw::vec4(0.0),0.0);
		DepositInitialCharge(x,N*charge);
	}

	return particleDensity*tw::Float(numToAdd);
}

tw::Float Species::AddDensityRandom(const LoadingData& theData)
{
	const tw::vec3 cellCenter = owner->Pos(theData.cell);
	//const tw::vec3 cellSize = owner->dPos(theData.cell);
	const tw::Float cellVolume = owner->dS(theData.cell,0);
	const tw::Float densToAdd = theData.densToAdd;
	const tw::vec3 thermalMomentum = theData.thermalMomentum;
	const tw::vec3 driftMomentum = theData.driftMomentum;
	//const bool neutralize = theData.neutralize;

	tw::vec3 r_primitive;
	Primitive q;
	tw::vec4 p;
	tw::Float particleDensity,fractionalNum;
	tw::Int numToAdd;

	particleDensity = theData.particleDensity;
	if (particleDensity <= minimumDensity)
		return 0.0;

	fractionalNum = densToAdd / particleDensity;
	numToAdd = tw::Int(fractionalNum);
	if (owner->uniformDeviate->Next()<(fractionalNum - tw::Float(numToAdd)))
		numToAdd++;

	if (numToAdd<1) return 0.0;

	for (tw::Int i=0;i<numToAdd;i++)
	{
		r_primitive.x = owner->uniformDeviate->Next() - 0.5;
		r_primitive.y = owner->uniformDeviate->Next() - 0.5;
		r_primitive.z = owner->uniformDeviate->Next() - 0.5;

		p[1] = thermalMomentum.x*owner->gaussianDeviate->Next();
		p[2] = thermalMomentum.y*owner->gaussianDeviate->Next();
		p[3] = thermalMomentum.z*owner->gaussianDeviate->Next();
		p += tw::vec4(0.0,driftMomentum);
		p[0] = sqrt(sqr(restMass) + Norm(p.spatial()));
		q.x[0] = 0.0; // always center of time cell
		q.x[1] = r_primitive.x;
		q.x[2] = r_primitive.y;
		q.x[3] = r_primitive.z;
		// At the beginning, timeLevel = 0.  So the particle's time will be t=-0.5*dt which can be viewed
		// as in the center of the temporal ghost cell.  The field data is known at t = 0, i.e., on the
		// high wall of the temporal ghost cell.
		q.cell = EncodeCell(theData.timeLevel,theData.cell.dcd1(),theData.cell.dcd2(),theData.cell.dcd3());
		const tw::vec4 x = PositionFromPrimitive(q);
		const tw::Float N = theData.GeometryFactor(x[1],cellCenter.x)*particleDensity*cellVolume;
		AddParticle(N,q,p,tw::vec4(0.0),0.0);
		DepositInitialCharge(x,N*charge);
	}

	return particleDensity*tw::Float(numToAdd);
}

void Species::BeginMoveWindow()
{
	// Throw out and transfer particles leaving on the left

	for (tw::Int i=0;i<particle.size();i++)
	{
		tw::Int ijk[4];
		DecodeCell(particle[i].q,ijk);
		ijk[3]--;
		particle[i].q.cell = EncodeCell(ijk[0],ijk[1],ijk[2],ijk[3]);
		if (!RefCellInSpatialDomain(particle[i].q))
		{
			if (owner->n0[3]!=MPI_PROC_NULL) // need to transfer particle
				mover->AddTransferParticle(particle[i]);
			particle[i] = particle.back();
			particle.pop_back();
			i--;
		}
	}
}

/// @brief create new particles on the right
void Species::FinishMoveWindow()
{
	if (owner->n1[3]==MPI_PROC_NULL)
		for (auto prof : profile)
			for (tw::Int j=1;j<=dim[2];j++)
				for (tw::Int i=1;i<=dim[1];i++)
				{
					LoadingData loadingData(*owner,distributionInCell,tw::cell(*this,i,j,dim[3]));
					loadingData.densToAdd = prof->GetValue(owner->Pos(i,j,dim[3]),*owner);
					loadingData.neutralize = owner->neutralize;
					loadingData.thermalMomentum = prof->thermalMomentum;
					loadingData.driftMomentum = prof->DriftMomentum(restMass);
					if (prof->variableCharge)
						loadingData.particleDensity = loadingData.densToAdd/(distributionInCell.x*distributionInCell.y*distributionInCell.z);
					else
						loadingData.particleDensity = targetDensity*prof->gammaBoost;
					if (prof->loadingMethod==tw::profile::loading::deterministic)
						AddDensity(loadingData);
					else
						AddDensityRandom(loadingData);
				}
}

bool Species::ReadInputFileDirective(const TSTreeCursor *curs0,const std::string& src)
{
	if (Module::ReadInputFileDirective(curs0,src))
		return true;

	if (tw::input::node_kind(curs0)=="assignment") {
		auto curs = tw::input::Cursor(curs0);
		ts_tree_cursor_goto_first_child(curs.get());
		if (tw::input::node_text(curs.get(),src) == "when density") {
			// needs to be run together as `particles per cell = 3 3 1 when density = 1.0`
			// so that the distributionInCell is known.
			tw::input::Float directive(&targetDensity);
			directive.Read(curs.get(),src,"when density",native);
			targetDensity /= (distributionInCell.x*distributionInCell.y*distributionInCell.z);
			return true;
		}
	}
	return false;	
}

void Species::ReadCheckpoint(std::ifstream& inFile)
{
	Module::ReadCheckpoint(inFile);
	inFile.read((char *)&numberCreated,sizeof(tw::Float));
	inFile.read((char *)&accelerationForceNow,sizeof(tw::Float));
	inFile.read((char *)&count,sizeof(tw::Int));

	particle.clear();

	tw::Int num;
	inFile.read((char *)&num,sizeof(tw::Int));
	for (tw::Int i=0;i<num;i++)
	{
		particle.emplace_back(0.0,Primitive(0,0.0,0.0,0.0,0.0),tw::vec4(0.0),tw::vec4(0.0),1,0.0);
		particle.back().ReadCheckpoint(inFile);
	}
}

void Species::WriteCheckpoint(std::ofstream& outFile)
{
	Module::WriteCheckpoint(outFile);
	outFile.write((char *)&numberCreated,sizeof(tw::Float));
	outFile.write((char *)&accelerationForceNow,sizeof(tw::Float));
	outFile.write((char *)&count,sizeof(tw::Int));

	tw::Int num = particle.size();
	outFile.write((char *)&num,sizeof(tw::Int));
	for (auto p : particle)
		p.WriteCheckpoint(outFile);
}


///////////////////////////
//  SPECIES DIAGNOSTICS  //
///////////////////////////


void Species::CalculateDensity(ScalarField& dens)
{
	Primitive q;
	weights_3D w;
	dens.Initialize(*this,owner);
	for (auto par : particle)
	{
		if (RefCellInSpatialDomain(par.q))
		{
			dens.StaticSpace::GetWeights(&w,par.q);
			dens.InterpolateOnto(par.number,w);
		}
	}
	transfer.clear();
	dens.SetBoundaryConditions(tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	dens.SetBoundaryConditions(tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	dens.SetBoundaryConditions(tw::grid::z,fld::dirichletCell,fld::dirichletCell);
	dens.DepositFromNeighbors();
	dens.ApplyFoldingCondition();
	dens.DivideCellVolume(*owner);
	dens.ApplyBoundaryCondition();
	dens.Smooth(*owner,smoothing,compensation);
}

void Species::CalculateEnergyDensity(ScalarField& dens)
{
	const tw::Float m0 = restMass;
	Primitive q;
	weights_3D w;
	dens.Initialize(*this,owner);
	for (auto par : particle)
	{
		if (RefCellInSpatialDomain(par.q))
		{
			dens.StaticSpace::GetWeights(&w,par.q);
			dens.InterpolateOnto(par.number * (par.p[0] - m0),w);
		}
	}
	transfer.clear();
	dens.SetBoundaryConditions(tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	dens.SetBoundaryConditions(tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	dens.SetBoundaryConditions(tw::grid::z,fld::dirichletCell,fld::dirichletCell);
	dens.DepositFromNeighbors();
	dens.ApplyFoldingCondition();
	dens.DivideCellVolume(*owner);
	dens.ApplyBoundaryCondition();
	dens.Smooth(*owner,smoothing,compensation);
}

void Species::CalculateQuantumParameter(ScalarField& dens)
{
	Primitive q;
	weights_3D w;
	dens.Initialize(*this,owner);
	for (auto par : particle)
	{
		if (RefCellInSpatialDomain(par.q))
		{
			dens.StaticSpace::GetWeights(&w,par.q);
			dens.InterpolateOnto(par.number*par.Qparam,w);
		}
	}
	transfer.clear();
	dens.SetBoundaryConditions(tw::grid::x,fld::dirichletCell,fld::dirichletCell);
	dens.SetBoundaryConditions(tw::grid::y,fld::dirichletCell,fld::dirichletCell);
	dens.SetBoundaryConditions(tw::grid::z,fld::dirichletCell,fld::dirichletCell);
	dens.DepositFromNeighbors();
	dens.ApplyFoldingCondition();
	dens.DivideCellVolume(*owner);
	dens.ApplyBoundaryCondition();
	dens.Smooth(*owner,smoothing,compensation);
}

tw::Float Species::KineticEnergy(const Region& theRgn)
{
	tw::Float ans = 0.0;
	const tw::Float m0 = restMass;
	for (auto par : particle)
		if (theRgn.Inside(owner->PositionFromPrimitive(par.q).spatial(),*owner))
			ans += par.number * (par.p[0] - m0);
	return ans;
}

tw::Float Species::ParticleNumber(const Region& theRgn)
{
	tw::Float ans = 0.0;
	for (auto par : particle)
		if (theRgn.Inside(owner->PositionFromPrimitive(par.q).spatial(),*owner))
			ans += par.number;
	return ans;
}

void Species::Report(Diagnostic& diagnostic)
{
	Module::Report(diagnostic);

	ScalarField temp;

	diagnostic.ReportNumber("N_"+name,particle.size(),false);
	diagnostic.ReportNumber("Kinetic_"+name,KineticEnergy(*diagnostic.theRgn),false);
	diagnostic.ReportNumber("Number_"+name,ParticleNumber(*diagnostic.theRgn),false);
	diagnostic.ReportNumber("Creation_"+name,numberCreated,false);
	CalculateDensity(temp);
	diagnostic.ReportField(name,temp,0,tw::dims::density,"$n_{\\rm "+name+"}$");
	CalculateEnergyDensity(temp);
	diagnostic.ReportField("u_"+name,temp,0,tw::dims::energy_density,"$u_{\\rm "+name+"}$");
	CalculateQuantumParameter(temp);
	diagnostic.ReportField("Q_"+name,temp,0,tw::dims::none,"$Q_{\\rm "+name+"}$");

	if (qo_j4!=NULL)
	{
		temp *= restMass;
		diagnostic.VolumeIntegral("tot_mass",temp,0);
		temp /= restMass;

		for (auto cell : InteriorCellRange(*this))
			temp(cell) = sqrt(fabs((*qo_j4)(cell,0) * temp(cell)));
		diagnostic.VolumeIntegral("overlap",temp,0); // square in post-processing to get population
	}

	for (auto par : particle)
		diagnostic.ReportParticle(par,restMass);
}
