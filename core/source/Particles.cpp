#include "sim.h"
#include "particles.h"
#include "fieldSolve.h"
#include "laserSolve.h"

////////////////////
// KINETICS CLASS //
////////////////////
// This is the top of the particle containment heierarchy
// Kinetics <- Species <- Particle

Kinetics::Kinetics(const std::string& name,Grid* theGrid) : Module(name,theGrid)
{
	typeCode = tw::module_type::kinetics;
	rho00.Initialize(*this,owner);
	sources = NULL;
	chi = NULL;
}

void Kinetics::Initialize()
{
	Module::Initialize();
	for (tw::Int i=0;i<submodule.size();i++)
		species.push_back((Species*)submodule[i]);
}

bool Kinetics::ValidSubmodule(Module* sub)
{
	return sub->typeCode==tw::module_type::species;
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
	tw::Int i,j;

	// rho00 must be prepared to receive deposition from incoming particles
	// rho00 array is always left in an unrefined post-deposition state
	// i.e., charge from a given particle should be on only one node
	// therefore charge has to be zeroed on one node after being sent to another
	rho00.DownwardDeposit(zAxis,1);
	for (StripIterator s(*this,3,strongbool::yes);s<s.end();++s)
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

	for (i=0;i<species.size();i++)
	{
		if (owner->stepNow%species[i]->sortPeriod==0)
			std::sort(species[i]->particle.begin(),species[i]->particle.end());
		species[i]->DispatchPush();
		species[i]->ApplyGlobalBoundaryConditions();
	}

	TransferParticles();

	Ionize();

	for (i=0;i<species.size();i++)
	{
		species[i]->GenerateParticles(false);
		species[i]->CleanParticleList();
	}

	if (sources && owner->neutralize)
		for (CellIterator cell(*this,true);cell<cell.end();++cell)
			(*sources)(cell,0) -= rho00(cell);

	// Take care of boundaries in source arrays within field solvers
	// otherwise every module that deposits something would have to do message passing
}

void Kinetics::TransferParticles()
{
	tw::Int i,j,displ,a;

	tw::Int dst,src;
	tw::Int numToSend,sendSize,recvSize;
	tw::Int offset;
	bool odd;

	std::valarray<char> inBuffer,outBuffer;
	TransferParticle *parPtr;

	std::vector<TransferParticle> accumulator;
	std::vector<tw::Int> tally;

	for (i=0;i<species.size();i++)
		species[i]->ComputeTransferParticleDestinations();

	for (a=1;a<=3;a++)
	{
		if (owner->localCells[a]>1)
		{
			for (displ=-1;displ<=1;displ+=2)
			{
				accumulator.clear();
				tally.clear();

				owner->strip[a].Shift(1,displ,&src,&dst);
				odd = owner->strip[a].Get_rank() % 2;

				// Load particles and tallies into standard containers for all species
				// These are the transfer particles that need to move along the given axis

				for (i=0;i<species.size();i++)
					species[i]->PrepareTransfer(accumulator,tally,a,displ);

				// Exchange message sizes

				numToSend = accumulator.size();
				sendSize =  sizeof(tw::Int)*species.size() + sizeof(TransferParticle)*numToSend;
				recvSize = 0;
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
				for (i=0;i<species.size();i++)
					((tw::Int*)&outBuffer[0])[i] = tally[i];
				if (numToSend)
				{
					parPtr = (TransferParticle*)(&outBuffer[species.size()*sizeof(tw::Int)]);
					for (i=0;i<numToSend;i++)
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
					for (i=0;i<species.size();i++)
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

	for (i=0;i<species.size();i++)
		species[i]->CollectTransfers();
}

void Kinetics::Ionize()
{
	// Call at start of pusher, because:
	// Consider N = particle number to be known at t = -1/2.  Use E and R at t = 0 to advance N to t = 1/2.
	// Then we push the new and old particles, and deposit J at t = 1/2 using the new particle number.
	// Here, we suppose the rates due to the envelope and wake are additive

	tw::Int i,s;
	weights_3D weights;
	tw::Float instant,peak,a2,probability;
	std::valarray<tw::Float> temp(6),Fp(8);
	std::valarray<tw::Float> jV(4);
	tw::vec3 r,E,vel,momentum;
	tw::Float gamma,m0,q0;
	Species *s1,*s2;

	// Find Laser Solver

	LaserSolver *theLaserSolver = NULL;
	for (i=0;i<owner->module.size();i++)
		if (dynamic_cast<LaserSolver*>(owner->module[i]))
			theLaserSolver = (LaserSolver*)owner->module[i];

	// Photoionization

	for (s=0;s<species.size();s++)
	{
		IonizationData& ionization = species[s]->ionization;
		std::vector<Particle>& particle = species[s]->particle;
		s1 = (Species*)owner->module[ionization.ionSpecies];
		s2 = (Species*)owner->module[ionization.electronSpecies];
		m0 = species[s]->restMass;
		q0 = species[s]->charge;

		if (ionization.ionizationModel!=tw::ionization_model::none)
		{
			for (i=0;i<particle.size();i++)
			{
				Particle& curr = particle[i]; // curr = particle being ionized
				r = owner->PositionFromPrimitive(curr.q);
				owner->GetWeights(&weights,curr.q);
				peak = instant = a2 = 0.0;
				if (species[s]->laser)
				{
					species[s]->laser->Interpolate(Fp,Element(6),weights);
					a2 = Fp[6];
					instant = theLaserSolver->polarizationType==circularPolarization ? (*species[s]->carrierFrequency)*sqrt(0.5*a2) : 0.0;
					peak = theLaserSolver->polarizationType==circularPolarization ? 0.0 : (*species[s]->carrierFrequency)*sqrt(a2);
				}
				species[s]->EM->Interpolate(temp,Element(0,2),weights);
				E.x = temp[0]; E.y = temp[1]; E.z = temp[2];
				instant += Magnitude(E);
				probability = ionization.Rate(instant,peak)*owner->dt;
				gamma = sqrt(1.0 + Norm(curr.p)/(m0*m0) + 0.5*sqr(q0)*a2/(m0*m0));
				vel = curr.p/(gamma*m0);
				if (probability > owner->uniformDeviate->Next())
				{
					throw tw::FatalError("Ionization not supported.");

					momentum = s1->restMass*gamma*vel;
					if (theLaserSolver)
						momentum += theLaserSolver->GetIonizationKick(a2,s1->charge,s1->restMass);
					s1->AddParticle(momentum,curr.q,curr.number);

					momentum = s2->restMass*gamma*vel;
					if (theLaserSolver)
						momentum += theLaserSolver->GetIonizationKick(a2,s2->charge,s2->restMass);
					s2->AddParticle(momentum,curr.q,curr.number);

					// compute ionization current (not centered, but this is usually small)
					if (!theLaserSolver)
					{
						jV[1] = 2.66e-5*ionization.ionizationPotential*curr.number/(Magnitude(E)*owner->dt);
						jV[2] = jV[1];
						jV[3] = jV[1];
						Normalize(E);
						jV[1] *= E.x;
						jV[2] *= E.y;
						jV[3] *= E.z;
						s2->sources->InterpolateOnto(jV,Element(1,3),weights);
					}

					particle[i] = particle.back();
					particle.pop_back();
					--i;
				}
			}
		}
	}
}

tw::vec3 LaserSolver::GetIonizationKick(const tw::Float& a2,const tw::Float& q0,const tw::Float& m0)
{
	tw::Float phase;
	tw::vec3 ans;
	if (polarizationType==circularPolarization)
	{
		phase = owner->uniformDeviate->Next()*2.0*pi;
		// remember "a" has been multiplied by sqrt(2) at the beginning
		ans.x = q0*sqrt(0.5*a2)*cos(phase);
		ans.y = q0*sqrt(0.5*a2)*sin(phase);
		ans.z = 0.25*q0*q0*a2/m0;
	}
	else
	{
		// for linear polarization, assume phase is at zero of vector potential (peak of field)
		ans.x = 0.0;
		ans.y = 0.0;
		ans.z = 0.25*q0*q0*a2/m0;
	}
	return ans;
}

void Kinetics::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "Kinetic ";
}

void Kinetics::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	cols.push_back(KineticEnergy(theRgn)); avg.push_back(false);
}

tw::Float Kinetics::KineticEnergy(const Region& theRgn)
{
	tw::Float ans = 0.0;
	tw::Float p2,m0;
	tw::vec3 pos;
	tw::Int i,j;
	Particle par;

	for (i=0;i<species.size();i++)
	{
		m0 = species[i]->restMass;
		for (j=0;j<species[i]->particle.size();j++)
		{
			par = species[i]->particle[j];
			pos = owner->PositionFromPrimitive(par.q);
			if (theRgn.Inside(pos,*owner))
			{
				p2 = Norm(par.p);
				tw::Float gamma = sqrt(1 + p2/(m0*m0));
				ans += par.number * m0 * (gamma - 1.0);
			}
		}
	}
	return ans;
}

void Kinetics::ReadData(std::ifstream& inFile)
{
	Module::ReadData(inFile);
	rho00.ReadData(inFile);
}

void Kinetics::WriteData(std::ofstream& outFile)
{
	Module::WriteData(outFile);
	rho00.WriteData(outFile);
}


////////////////////
// PARTICLE CLASS //
////////////////////


Particle::Particle()
{
	q.x[0] = 0.0;
	q.x[1] = 0.0;
	q.x[2] = 0.0;
	q.cell = 0;
	number = 1.0;
	aux1 = 0.0;
	aux2 = 0.0;
}

Particle::Particle(const tw::vec3& p,const Primitive& q,const float number,const float aux1,const float aux2)
{
	this->p = p;
	this->q = q;
	this->number = number;
	this->aux1 = aux1;
	this->aux2 = aux2;
}

void Particle::ReadData(std::ifstream& inFile)
{
	inFile.read((char *)this,sizeof(Particle));
}

void Particle::WriteData(std::ofstream& outFile)
{
	outFile.write((char *)this,sizeof(Particle));
}


///////////////////
// SPECIES CLASS //
///////////////////


Species::Species(const std::string& name,Grid* theGrid) : Module(name,theGrid)
{
	typeCode = tw::module_type::species;
	restMass = 1.0;
	charge = -1.0;
	distributionInCell = tw::vec3(2.0,2.0,2.0);
	minimumDensity = 0.0;
	targetDensity = 1.0;
	accelerationTime = accelerationImpulse = accelerationForceNow = 0.0;
	meanFreePath = 0.0; // means infinity
	count = 1;
	sortPeriod = 10;
	emissionTemp = tw::vec3(0.1,0.1,0.1);

	for (tw::Int i=1;i<=3;i++)
	{
		bc0[i] = theGrid->bc0[i];
		bc1[i] = theGrid->bc1[i];
	}

	mobile = true;
	radiationDamping = false;

	EM = NULL;
	sources = NULL;
	laser = NULL;
	chi = NULL;
	carrierFrequency = NULL;
	polarizationType = NULL;
	rho00 = NULL;
	ESField = NULL;
	qo_j4 = NULL;
}

void Species::Initialize()
{
	tw::Int i;
	Module::Initialize();
	for (i=0;i<phaseSpacePlot.size();i++)
		phaseSpacePlot[i]->SetupGeometry(owner->gridGeometry);
	if (ionization.ionizationModel!=tw::ionization_model::none)
	{
		ionization.Initialize(owner->unitDensityCGS,carrierFrequency);
		ionization.electronSpecies = owner->FindModule(ionization.electron_name);
		ionization.ionSpecies = owner->FindModule(ionization.ion_name);
	}
	if (!owner->restarted)
	{
		GenerateParticles(true);
		CleanParticleList();
	}
	if (!owner->neutralize)
		CopyFieldData(*sources,0,*rho00,0); // accepting redundant operations
}

void Species::WarningMessage(std::ostream *theStream)
{
	Module::WarningMessage(theStream);
	if (bc0[1]==absorbing || bc1[1]==absorbing ||
		bc0[2]==absorbing || bc1[2]==absorbing ||
		bc0[3]==absorbing || bc1[3]==absorbing)
		(*theStream) << "WARNING: " << name << " will be absorbed." << std::endl;
}

bool Species::InspectResource(void* resource,const std::string& description)
{
	if (description=="electrostatic:E")
	{
		ESField = (Vec3Field*)resource;
		return true;
	}

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

void Species::AddParticle(const tw::vec3& p,const Primitive& q,const float& number)
{
	particle.push_back(Particle(p,q,number,count,owner->strip[0].Get_rank()));
	count++; // this serves as a unique identifier for particles irrespective of any additions or deletions in the vector
}

void Species::AddParticle(const TransferParticle& xfer)
{
	Primitive q;
	tw::vec3 x = tw::vec3(xfer.x[1],xfer.x[2],xfer.x[3]);
	tw::vec3 p = tw::vec3(xfer.p[1],xfer.p[2],xfer.p[3]);
	SetPrimitiveWithPosition(q,x); // particle must be in extended domain
	particle.push_back(Particle(p,q,xfer.number,xfer.aux1,xfer.aux2));
	// don't update count because the transfer particle already has its identifier in aux1 and aux2
}

void Species::AddTransferParticle(const Particle& src)
{
	TransferParticle dest;
	dest.dst[0] = owner->strip[0].Get_rank();
	for (tw::Int i=1;i<=3;i++)
		dest.dst[i] = 0;
	dest.x = tw::vec4(0.0,PositionFromPrimitive(src.q));
	dest.p = tw::vec4(0.0,src.p);
	dest.number = src.number;
	dest.aux1 = src.aux1;
	dest.aux2 = src.aux2;
	transfer.push_back(dest);
}

void Species::CleanParticleList()
{
	tw::Int i;
	for (i=0;i<particle.size();i++)
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



void Species::ApplyGlobalBoundaryConditions()
{
	// Must happen before particle message passing and before particle destination calculation.
	// All particles we want to keep should be transformed so they are in the global domain,
	// except for particles crossing a periodic boundary.
	// Result of the transformation can result in movement to new local domain.
	// If absorption is desired leave particle out of global domain (destination calculation will send to MPI_PROC_NULL).
	bool didLeave;
	tw::Float extremum;
	tw::Int ax,i,displ,src,dst;
	tw_boundary_spec boundaryCondition;

	for (displ=-1;displ<=1;displ+=2)
		for (ax=1;ax<=3;ax++)
		{
			boundaryCondition = displ<0 ? bc0[ax] : bc1[ax];
			owner->strip[ax].Shift(1,displ,&src,&dst);
			if (dst==MPI_PROC_NULL && owner->localCells[ax]!=1)
			{
				for (i=0;i<transfer.size();i++)
				{
					extremum = displ==1 ? corner[ax-1] + size[ax-1] : corner[ax-1];
					didLeave = displ==1 ? transfer[i].x[ax] >= extremum : transfer[i].x[ax] < extremum;
					if (didLeave)
					{
						switch (boundaryCondition)
						{
							case cyclic:
								// must delay particle wrap-around
								// typically this case will never be encountered since dst!=MPI_PROC_NULL
								break;
							case absorbing:
								// do nothing; let it leave
								break;
							case emitting:
								transfer[i].x[ax] = extremum - tw::Float(displ)*0.5*spacing[ax-1];
								transfer[i].p[1] = emissionTemp.x*owner->gaussianDeviate->Next();
								transfer[i].p[2] = emissionTemp.y*owner->gaussianDeviate->Next();
								transfer[i].p[3] = emissionTemp.z*owner->gaussianDeviate->Next();
								transfer[i].p[ax] = -tw::Float(displ)*fabs(transfer[i].p[ax]);
								break;
							case reflecting:
							case axisymmetric:
								transfer[i].x[ax] += 2.0*(extremum - transfer[i].x[ax]);
								transfer[i].p[ax] *= -1.0;
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



void Species::ComputeTransferParticleDestinations()
{
	// This must be called before message passing begins, and after global boundary conditions are imposed.
	// The destination domain for every particle on the transfer list has to be computed.
	// If the particle leaves the entire system the destination rank is set to MPI_PROC_NULL (only possible for absorbing boundaries).
	// It is assumed that the destination domain is set to the "current" domain upon entry.
	// The destination is encoded as (rank,x-motion,y-motion,z-motion) with i-motion one of -1,0,1
	tw::Int i,ax,dest_coords[4];

	for (i=0;i<transfer.size();i++)
	{
		for (ax=1;ax<=3;ax++)
		{
			if (transfer[i].x[ax] < corner[ax-1])
			{
				if (owner->n0[ax]==MPI_PROC_NULL)
					transfer[i].dst[0] = MPI_PROC_NULL;
				else
					transfer[i].dst[ax] = -1;
			}
			if (transfer[i].x[ax] >= corner[ax-1]+size[ax-1])
			{
				if (owner->n1[ax]==MPI_PROC_NULL)
					transfer[i].dst[0] = MPI_PROC_NULL;
				else
					transfer[i].dst[ax] = 1;
			}
			dest_coords[ax] = owner->domainIndex[ax] + transfer[i].dst[ax];
			// allow for periodicity
			if (dest_coords[ax]<0) dest_coords[ax] = owner->domains[ax]-1;
			if (dest_coords[ax]>=owner->domains[ax]) dest_coords[ax] = 0;
			// if destination is this node stop the movement (e.g., single periodic node)
			if (dest_coords[ax]==owner->domainIndex[ax])
				transfer[i].dst[ax] = 0; // do not send
		}
		if (transfer[i].dst[0]!=MPI_PROC_NULL)
			transfer[i].dst[0] = owner->strip[0].Cart_rank(dest_coords[1],dest_coords[2],dest_coords[3]);
	}
}

void Species::PrepareTransfer(std::vector<TransferParticle>& accumulator,std::vector<tw::Int>& tally,tw::Int ax,tw::Int displ)
{
	// displ is -1 for downward transfer, +1 for upward transfer
	tw::Int i;
	tw::Int numToSend = 0;

	for (i=0;i<transfer.size();i++)
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
	tw::Int i;

	for (i=0;i<numToReceive;i++)
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
	tw::Int i,ax;
	for (i=0;i<transfer.size();i++)
		if (transfer[i].dst[0]==owner->strip[0].Get_rank())
		{
			// If global containment was lost, assume periodicity.
			for (ax=1;ax<=3;ax++)
			{
				if (transfer[i].x[ax]<globalCorner[ax-1])
					transfer[i].x[ax] += globalSize[ax-1];
				if (transfer[i].x[ax]>=globalCorner[ax-1]+globalSize[ax-1])
					transfer[i].x[ax] -= globalSize[ax-1];
			}
			AddParticle(transfer[i]);
		}
	transfer.clear();
}


/////////////////////////
//  PARTICLE CREATION  //
/////////////////////////




void Species::DepositInitialCharge(const tw::vec3& pos,tw::Float macroCharge)
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
	tw::Int i,j,k,s;
	bool timeGate;

	LoadingData loadingData(distributionInCell);

	for (s=0;s<profile.size();s++)
	{
		if (profile[s]->timingMethod==triggeredProfile)
			timeGate = owner->elapsedTime>=profile[s]->t0 && !profile[s]->wasTriggered;
		if (profile[s]->timingMethod==maintainedProfile)
			timeGate = owner->elapsedTime>=profile[s]->t0 && owner->elapsedTime<=profile[s]->t1;
		if ( timeGate )
		{
			profile[s]->wasTriggered = true;
			loadingData.driftMomentum = profile[s]->driftMomentum;
			loadingData.neutralize = profile[s]->neutralize;
			loadingData.thermalMomentum = profile[s]->thermalMomentum;
			for (i=1;i<=owner->Dim(1);i++)
				for (j=1;j<=owner->Dim(2);j++)
					for (k=1;k<=owner->Dim(3);k++)
					{
						loadingData.i = i;
						loadingData.j = j;
						loadingData.k = k;
						loadingData.densToAdd = profile[s]->GetValue(owner->Pos(i,j,k),*owner);

						if (loadingData.densToAdd>0.0)
						{
							if (profile[s]->timingMethod==maintainedProfile && !init)
							{
								// Hold a constant charge
								// (comment out for constant current)
								if (profile[s]->neutralize) // don't test loadingData.neutralize
									loadingData.densToAdd = -(*sources)(i,j,k,0)/charge; // hold this region neutral
								else
									loadingData.densToAdd -= (*sources)(i,j,k,0)/charge; // hold at target density (assumes no other species)
								loadingData.neutralize = false;
							}
							if (qo_j4!=NULL && profile[s]->loadingMethod==statistical && !profile[s]->variableCharge)
							{
								loadingData.densToAdd = fabs((*qo_j4)(i,j,k,0));
							}

							if (loadingData.densToAdd < 0.0) // constant charge block may lead to n<0
								loadingData.densToAdd = 0.0;

							// INJECT

							if (profile[s]->variableCharge)
								loadingData.particleDensity = loadingData.densToAdd/(distributionInCell.x*distributionInCell.y*distributionInCell.z);
							else
								loadingData.particleDensity = targetDensity;
							if (profile[s]->loadingMethod==deterministic)
								AddDensity(loadingData);
							else
								AddDensityRandom(loadingData);
						}
					}
		}
	}
}

LoadingData::LoadingData(const tw::vec3& distribution)
{
	i = j = k = 1;
	densToAdd = 0.0;
	densNow = 0.0;
	particleDensity = 0.0;
	neutralize = true;

	tw::Int nx = tw::Int(distribution.x);
	tw::Int ny = tw::Int(distribution.y);
	tw::Int nz = tw::Int(distribution.z);

	pointsInSubGrid = nx*ny*nz;
	subGrid.resize(pointsInSubGrid);
	tw::Int ii,jj,kk;
	for (kk=0;kk<nz;kk++)
		for (jj=0;jj<ny;jj++)
			for (ii=0;ii<nx;ii++)
				subGrid[kk*nx*ny + jj*nx + ii] = tw::vec3
					(
						-0.5 + 0.5/tw::Float(nx) + tw::Float(ii)/tw::Float(nx),
						-0.5 + 0.5/tw::Float(ny) + tw::Float(jj)/tw::Float(ny),
						-0.5 + 0.5/tw::Float(nz) + tw::Float(kk)/tw::Float(nz)
					);
}

tw::Float Species::AddDensity(const LoadingData& theData)
{
	const tw::vec3 cellCenter = owner->Pos(theData.i,theData.j,theData.k);
	const tw::vec3 cellSize = owner->dPos(theData.i,theData.j,theData.k);
	const tw::Float cellVolume = owner->dS(theData.i,theData.j,theData.k,0);
	const tw::Float densToAdd = theData.densToAdd;
	const tw::Float densNow = theData.densNow;
	const tw::vec3 thermalMomentum = theData.thermalMomentum;
	const tw::vec3 driftMomentum = theData.driftMomentum;
	const bool neutralize = theData.neutralize;
	const tw::Int pointsInSubGrid = theData.pointsInSubGrid;

	tw::Int i,j;
	tw::vec3 r_primitive,r,p;
	Primitive q;
	tw::Float particleDensity,geometryFactor,C0=1.0,C1=0.0,C2=0.0,Nr2;
	tw::Int numToAdd,numNow;
	std::valarray<tw::vec3> initialMomenta;

	particleDensity = theData.particleDensity;
	if (particleDensity <= minimumDensity)
		return 0.0;

	numNow = TruncateUnlessClose( densNow / particleDensity);
	numToAdd = TruncateUnlessClose( (densNow + densToAdd) / particleDensity) - numNow;

	if (numToAdd<1) return 0.0;

	p = tw::vec3(0,0,0);
	initialMomenta.resize(numToAdd);
	for (i=0;i<numToAdd;i++)
	{
		initialMomenta[i].x = thermalMomentum.x*owner->gaussianDeviate->Next();
		initialMomenta[i].y = thermalMomentum.y*owner->gaussianDeviate->Next();
		initialMomenta[i].z = thermalMomentum.z*owner->gaussianDeviate->Next();
		p += initialMomenta[i];
	}
	p /= tw::Float(numToAdd);
	initialMomenta -= p - driftMomentum;

	Nr2 = sqr(distributionInCell.x);
	if (owner->gridGeometry==cylindrical)
	{
		C0 = 0.0;
		C1 = 1.0;
		C2 = 0.0;
		if (owner->n0[1]==MPI_PROC_NULL)
		{
			if (Nr2==1.0)
			{
				C0 = 1.0;
				C1 = 0.0;
				if (theData.i==1)
					C0 = 4208.0/5951.0;
				if (theData.i==2)
					C0 = 18152.0/17853.0;
				if (theData.i==3)
					C0 = 29704.0/29755.0;
				if (theData.i==4)
					C0 = 5952.0/5951.0;

			}
			else
			{
				if (theData.i==1)
				{
					C0 = 0.0;
					C1 = (7.0 - 22.0*Nr2)/(24.0*Nr2*(Nr2-1.0));
					C2 = (5.0/8.0) * (4.0*Nr2*Nr2-1.0) / ((Nr2-1.0)*(4.0*Nr2-1.0));
				}
			}
		}
	}

	for (i=0;i<numToAdd;i++)
	{
		j = numNow + i;
		while (j>=pointsInSubGrid) j -= pointsInSubGrid;
		r_primitive = theData.subGrid[j];
		r = r_primitive*cellSize + cellCenter;
		p = initialMomenta[i];

		geometryFactor = fabs(C0 + SafeDiv(C1*r.x,cellCenter.x) + sqr(SafeDiv(sqrt(C2)*r.x,cellCenter.x)));

		q.x[0] = r_primitive.x;
		q.x[1] = r_primitive.y;
		q.x[2] = r_primitive.z;
		q.cell = EncodeCell(theData.i,theData.j,theData.k);
		AddParticle(p,q,geometryFactor*particleDensity*cellVolume);
		DepositInitialCharge(r,geometryFactor*particleDensity*cellVolume*charge);
	}

	return particleDensity*tw::Float(numToAdd);
}

tw::Float Species::AddDensityRandom(const LoadingData& theData)
{
	const tw::vec3 cellCenter = owner->Pos(theData.i,theData.j,theData.k);
	const tw::vec3 cellSize = owner->dPos(theData.i,theData.j,theData.k);
	const tw::Float cellVolume = owner->dS(theData.i,theData.j,theData.k,0);
	const tw::Float densToAdd = theData.densToAdd;
	const tw::vec3 thermalMomentum = theData.thermalMomentum;
	const tw::vec3 driftMomentum = theData.driftMomentum;
	const bool neutralize = theData.neutralize;

	tw::Int i;
	tw::vec3 r_primitive,r,p;
	Primitive q;
	tw::Float particleDensity,fractionalNum,geometryFactor,C0=1.0,C1=0.0,C2=0.0,Nr2;
	tw::Int numToAdd;

	particleDensity = theData.particleDensity;
	if (particleDensity <= minimumDensity)
		return 0.0;

	fractionalNum = densToAdd / particleDensity;
	numToAdd = tw::Int(fractionalNum);
	if (owner->uniformDeviate->Next()<(fractionalNum - tw::Float(numToAdd)))
		numToAdd++;

	if (numToAdd<1) return 0.0;

	Nr2 = sqr(distributionInCell.x);
	if (owner->gridGeometry==cylindrical)
	{
		C0 = 0.0;
		C1 = 1.0;
		C2 = 0.0;
		if (owner->n0[1]==MPI_PROC_NULL)
		{
			if (Nr2==1.0)
			{
				C0 = 1.0;
				C1 = 0.0;
				if (theData.i==1)
					C0 = 4208.0/5951.0;
				if (theData.i==2)
					C0 = 18152.0/17853.0;
				if (theData.i==3)
					C0 = 29704.0/29755.0;
				if (theData.i==4)
					C0 = 5952.0/5951.0;

			}
			else
			{
				if (theData.i==1)
				{
					C0 = 0.0;
					C1 = (7.0 - 22.0*Nr2)/(24.0*Nr2*(Nr2-1.0));
					C2 = (5.0/8.0) * (4.0*Nr2*Nr2-1.0) / ((Nr2-1.0)*(4.0*Nr2-1.0));
				}
			}
		}
	}

	for (i=0;i<numToAdd;i++)
	{
		r_primitive.x = owner->uniformDeviate->Next() - 0.5;
		r_primitive.y = owner->uniformDeviate->Next() - 0.5;
		r_primitive.z = owner->uniformDeviate->Next() - 0.5;
		r = r_primitive*cellSize + cellCenter;

		p.x = thermalMomentum.x*owner->gaussianDeviate->Next();
		p.y = thermalMomentum.y*owner->gaussianDeviate->Next();
		p.z = thermalMomentum.z*owner->gaussianDeviate->Next();
		p += driftMomentum;

		geometryFactor = fabs(C0 + SafeDiv(C1*r.x,cellCenter.x) + sqr(SafeDiv(sqrt(C2)*r.x,cellCenter.x)));

		q.x[0] = r_primitive.x;
		q.x[1] = r_primitive.y;
		q.x[2] = r_primitive.z;
		q.cell = EncodeCell(theData.i,theData.j,theData.k);
		AddParticle(p,q,geometryFactor*particleDensity*cellVolume);
		DepositInitialCharge(r,geometryFactor*particleDensity*cellVolume*charge);
	}

	return particleDensity*tw::Float(numToAdd);
}

void Species::MoveWindow()
{
}

void Species::BeginMoveWindow()
{
	tw::Int i;

	// Throw out and transfer particles leaving on the left

	for (i=0;i<particle.size();i++)
	{
		particle[i].q.x[2] -= 1.0;
		if (PrimitiveInDomain(particle[i].q))
			MinimizePrimitive(particle[i].q);
		else
		{
			if (owner->n0[3]!=MPI_PROC_NULL) // need to transfer particle
				AddTransferParticle(particle[i]);
			particle[i] = particle.back();
			particle.pop_back();
			i--;
		}
	}
}

void Species::FinishMoveWindow()
{
	tw::Int i,j,s;

	for (i=0;i<transfer.size();i++)
	{
		if (transfer[i].dst[0]==owner->strip[0].Get_rank())
			AddParticle(transfer[i]);
	}
	transfer.clear();

	// Create new particles on the right

	if (owner->n1[3]==MPI_PROC_NULL)
	{
		LoadingData loadingData(distributionInCell);
		for (s=0;s<profile.size();s++)
		{
			loadingData.neutralize = profile[s]->neutralize;
			loadingData.thermalMomentum = profile[s]->thermalMomentum;
			loadingData.driftMomentum = profile[s]->driftMomentum;
			for (j=1;j<=dim[2];j++)
				for (i=1;i<=dim[1];i++)
				{
					loadingData.i = i;
					loadingData.j = j;
					loadingData.k = dim[3];
					loadingData.densToAdd = profile[s]->GetValue(owner->Pos(i,j,dim[3]),*owner);
					if (profile[s]->variableCharge)
						loadingData.particleDensity = loadingData.densToAdd/(distributionInCell.x*distributionInCell.y*distributionInCell.z);
					else
						loadingData.particleDensity = targetDensity;
					if (profile[s]->loadingMethod==deterministic)
						AddDensity(loadingData);
					else
						AddDensityRandom(loadingData);
				}
		}
	}
}

void Species::ReadInputFileDirective(std::stringstream& inputString,const std::string& command)
{
	tw::Int i;
	std::string word;

	ionization.ReadInputFileDirective(inputString,command);

	if (command=="sort") // eg, sort period = 100
	{
		inputString >> word >> word >> sortPeriod;
	}
	if (command=="mobile") // eg, mobile = true
	{
		inputString >> word >> word;
		mobile = (word=="true" || word=="yes" || word=="on");
	}
	if (command=="radiation") // eg, radiation damping = true
	{
		inputString >> word >> word >> word;
		radiationDamping = (word=="true" || word=="yes" || word=="on");
	}
	if (command=="mean") // eg, mean free path = 1.0
		inputString >> word >> word >> word >> meanFreePath;
	if (command=="mass") // eg, mass = 1.0
		inputString >> word >> restMass;
	if (command=="charge") // eg, charge = 1.0
		inputString >> word >> charge;
	if (command=="particles") // eg, particles per cell = 3 3 1 when density = 1.0
	{
		inputString >> word >> word >> word;
		inputString >> distributionInCell.x >> distributionInCell.y >> distributionInCell.z;
		inputString >> word >> word >> word;
		inputString >> targetDensity;
		targetDensity /= distributionInCell.x*distributionInCell.y*distributionInCell.z;
	}
	if (command=="minimum") // eg, minimum density = 1e-4
		inputString >> word >> word >> minimumDensity;
	if (command=="emission") // eg, emission temperature = 0.1 0.1 0.1
		inputString >> word >> word >> emissionTemp.x >> emissionTemp.y >> emissionTemp.z;
	if (command=="accelerate") // eg, accelerate to 100.0 in 10.0
		inputString >> word >> accelerationImpulse >> word >> accelerationTime;
	if (command=="xboundary" || command=="yboundary" || command=="zboundary" ) // eg, xboundary = emitting emitting
		tw::input::ReadBoundaryTerm(bc0,bc1,inputString,command);
}

bool Species::ReadQuasitoolBlock(const std::vector<std::string>& preamble,std::stringstream& inputString)
{
	std::string key(preamble[0]);
	std::string requested_name(preamble.back());
	if (key=="phase" && requested_name==name) // eg, new phase space plot for electrons { ... }
	{
		phaseSpacePlot.push_back(new PhaseSpaceDescriptor(owner->clippingRegion));
		phaseSpacePlot.back()->ReadInputFile(inputString);
		return true;
	}
	if (key=="orbit" && requested_name==name) // eg, new orbit diagnostic for electrons { ... }
	{
		orbitDiagnostic.push_back(new ParticleOrbitDescriptor(owner->clippingRegion));
		orbitDiagnostic.back()->ReadInputFile(inputString);
		return true;
	}
	if (key=="detector" && requested_name==name) // eg, new detector diagnostic for electrons { ... }
	{
		detector.push_back(new ParticleDetectorDescriptor(owner->clippingRegion));
		detector.back()->ReadInputFile(inputString);
		return true;
	}
	return false;
}

void Species::ReadData(std::ifstream& inFile)
{
	tw::Int i,num;

	Module::ReadData(inFile);
	ionization.ReadData(inFile);
	inFile.read((char *)&restMass,sizeof(tw::Float));
	inFile.read((char *)&charge,sizeof(tw::Float));
	inFile.read((char *)&emissionTemp,sizeof(tw::vec3));
	inFile.read((char *)&targetDensity,sizeof(tw::Float));
	inFile.read((char *)&minimumDensity,sizeof(tw::Float));
	inFile.read((char *)&accelerationImpulse,sizeof(tw::Float));
	inFile.read((char *)&accelerationTime,sizeof(tw::Float));
	inFile.read((char *)&accelerationForceNow,sizeof(tw::Float));
	inFile.read((char *)&count,sizeof(tw::Int));
	inFile.read((char *)&sortPeriod,sizeof(tw::Int));
	inFile.read((char *)&distributionInCell,sizeof(tw::vec3));
	inFile.read((char *)bc0,sizeof(tw_boundary_spec)*4);
	inFile.read((char *)bc1,sizeof(tw_boundary_spec)*4);
	inFile.read((char *)&mobile,sizeof(bool));
	inFile.read((char *)&radiationDamping,sizeof(bool));
	inFile.read((char *)&meanFreePath,sizeof(tw::Float));

	inFile.read((char *)&num,sizeof(tw::Int));
	(*owner->tw_out) << "Read " << num << " phase space plots" << std::endl;
	for (i=0;i<num;i++)
	{
		phaseSpacePlot.push_back(new PhaseSpaceDescriptor(owner->clippingRegion));
		phaseSpacePlot.back()->ReadData(inFile);
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	(*owner->tw_out) << "Read " << num << " orbit diagnostics" << std::endl;
	for (i=0;i<num;i++)
	{
		orbitDiagnostic.push_back(new ParticleOrbitDescriptor(owner->clippingRegion));
		orbitDiagnostic.back()->ReadData(inFile);
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	(*owner->tw_out) << "Read " << num << " detector diagnostics" << std::endl;
	for (i=0;i<num;i++)
	{
		detector.push_back(new ParticleDetectorDescriptor(owner->clippingRegion));
		detector.back()->ReadData(inFile);
	}

	inFile.read((char *)&num,sizeof(tw::Int));
	(*owner->tw_out) << "Read " << num << " particles" << std::endl;
	Particle newParticle;
	for (i=0;i<num;i++)
	{
		particle.push_back(newParticle);
		particle[i].ReadData(inFile);
	}
}

void Species::WriteData(std::ofstream& outFile)
{
	tw::Int i;

	Module::WriteData(outFile);
	ionization.WriteData(outFile);
	outFile.write((char *)&restMass,sizeof(tw::Float));
	outFile.write((char *)&charge,sizeof(tw::Float));
	outFile.write((char *)&emissionTemp,sizeof(tw::vec3));
	outFile.write((char *)&targetDensity,sizeof(tw::Float));
	outFile.write((char *)&minimumDensity,sizeof(tw::Float));
	outFile.write((char *)&accelerationImpulse,sizeof(tw::Float));
	outFile.write((char *)&accelerationTime,sizeof(tw::Float));
	outFile.write((char *)&accelerationForceNow,sizeof(tw::Float));
	outFile.write((char *)&count,sizeof(tw::Int));
	outFile.write((char *)&sortPeriod,sizeof(tw::Int));
	outFile.write((char *)&distributionInCell,sizeof(tw::vec3));
	outFile.write((char *)bc0,sizeof(tw_boundary_spec)*4);
	outFile.write((char *)bc1,sizeof(tw_boundary_spec)*4);
	outFile.write((char *)&mobile,sizeof(bool));
	outFile.write((char *)&radiationDamping,sizeof(bool));
	outFile.write((char *)&meanFreePath,sizeof(tw::Float));

	i = phaseSpacePlot.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<phaseSpacePlot.size();i++)
		phaseSpacePlot[i]->WriteData(outFile);

	i = orbitDiagnostic.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<orbitDiagnostic.size();i++)
		orbitDiagnostic[i]->WriteData(outFile);

	i = detector.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<detector.size();i++)
		detector[i]->WriteData(outFile);

	i = particle.size();
	outFile.write((char *)&i,sizeof(tw::Int));
	for (i=0;i<particle.size();i++)
		particle[i].WriteData(outFile);
}


///////////////////////////
//  SPECIES DIAGNOSTICS  //
///////////////////////////


void Species::CalculateDensity(ScalarField& dens)
{
	tw::Int i;
	Primitive q;
	weights_3D w;
	dens.Initialize(*this,owner);
	for (i=0;i<particle.size();i++)
	{
		if (PrimitiveInDomain(particle[i].q))
		{
			dens.GetWeights(&w,particle[i].q);
			dens.InterpolateOnto(particle[i].number,w);
		}
	}
	transfer.clear();
	dens.SetBoundaryConditions(xAxis,dirichletCell,dirichletCell);
	dens.SetBoundaryConditions(yAxis,dirichletCell,dirichletCell);
	dens.SetBoundaryConditions(zAxis,dirichletCell,dirichletCell);
	dens.DepositFromNeighbors();
	dens.ApplyFoldingCondition();
	dens.DivideCellVolume(*owner);
	dens.ApplyBoundaryCondition();
	if (owner->smoothing>0)
		dens.Smooth(*owner,owner->smoothing,owner->compensation);
}

void Species::EnergyHeadings(std::ofstream& outFile)
{
	outFile << "N_" << name << " ";
	if (qo_j4!=NULL)
		outFile << "tot_mass overlap ";
}

void Species::EnergyColumns(std::vector<tw::Float>& cols,std::vector<bool>& avg,const Region& theRgn)
{
	cols.push_back(particle.size()); avg.push_back(false);

	tw::Int i,j,k;
	tw::Int x0,x1,y0,y1,z0,z1;
	ScalarField dens;
	tw::Float overlap,tot_mass;

	if (qo_j4!=NULL)
	{
		CalculateDensity(dens);
		tot_mass = overlap = 0.0;
		theRgn.GetLocalCellBounds(&x0,&x1,&y0,&y1,&z0,&z1);
		for (k=z0;k<=z1;k++)
			for (j=y0;j<=y1;j++)
				for (i=x0;i<=x1;i++)
				{
					tot_mass += restMass * dens(i,j,k) * owner->dS(i,j,k,0);
					overlap += sqrt(fabs((*qo_j4)(i,j,k,0) * dens(i,j,k))) * owner->dS(i,j,k,0);
				}
		cols.push_back(tot_mass); avg.push_back(false);
		cols.push_back(overlap); avg.push_back(false); // square in post-processing to get population
	}
}

void Species::BoxDiagnosticHeader(GridDataDescriptor* box)
{
	owner->WriteBoxDataHeader(name,box);
}

void Species::BoxDiagnose(GridDataDescriptor* box)
{
	// Deposit this species onto a scratch array

	ScalarField dens;
	CalculateDensity(dens);

	// Write out the data

	owner->WriteBoxData(name,box,&dens(0,0,0),dens.Stride());
}

tw::Float Species::GetPhaseSpaceData(Particle* p,std::string& what)
{
	tw::vec3 pos,mom;
	pos = owner->PositionFromPrimitive(p->q);
	owner->CurvilinearToCartesian(&pos);
	mom = p->p;

	if (what=="x")
		return pos.x;
	if (what=="y")
		return pos.y;
	if (what=="z")
		return pos.z - (owner->movingWindow ? owner->elapsedTime : 0.0);
	if (what=="px")
		return mom.x;
	if (what=="py")
		return mom.y;
	if (what=="pz")
		return mom.z;
	if (what=="gbx")
		return mom.x/restMass;
	if (what=="gby")
		return mom.y/restMass;
	if (what=="gbz")
		return mom.z/restMass;
	if (what=="mass")
		return sqrt(sqr(restMass) + Norm(mom));
	if (what=="energy")
		return (sqrt(1 + Norm(mom)/sqr(restMass)) - 1.0)*restMass;
	return 0.0;
}

void Species::CustomDiagnose()
{
	float data,gamma;
	tw::Int i,j,s,pts;
	std::stringstream fileName;
	tw::vec3 position,phaseSpaceSize,phaseSpacePos,lastPos,vel;
	weights_3D weights;
	std::vector<float> parData;
	std::valarray<float> parBuffer;

	tw::Int master = 0;

	ScalarField field;

	if (owner->IsFirstStep() && !(owner->appendMode&&owner->restarted) && owner->strip[0].Get_rank()==master)
	{
		for (i=0;i<phaseSpacePlot.size();i++)
		{
			fileName.str("");
			fileName << phaseSpacePlot[i]->filename << ".dvdat";
			phaseSpacePlot[i]->outFile.open(fileName.str().c_str(),std::ios::binary);
			WriteDVHeader(phaseSpacePlot[i]->outFile,2,phaseSpacePlot[i]->hDim,phaseSpacePlot[i]->vDim,1,
				phaseSpacePlot[i]->min.x,phaseSpacePlot[i]->max.x,
				phaseSpacePlot[i]->min.y,phaseSpacePlot[i]->max.y,0.0,1.0);
			phaseSpacePlot[i]->outFile.close();
		}
		for (i=0;i<orbitDiagnostic.size();i++)
		{
			fileName.str("");
			fileName << orbitDiagnostic[i]->filename << ".dvpar";
			orbitDiagnostic[i]->outFile.open(fileName.str().c_str(),std::ios::binary);
			orbitDiagnostic[i]->outFile.close();
		}
		for (i=0;i<detector.size();i++)
		{
			fileName.str("");
			fileName << detector[i]->filename << ".txt";
			detector[i]->outFile.open(fileName.str().c_str(),std::ios::binary);
			detector[i]->outFile.close();
		}
	}

	for (s=0;s<phaseSpacePlot.size();s++)
	{
		if (phaseSpacePlot[s]->WriteThisStep(owner->elapsedTime,owner->dt,owner->stepNow))
		{
			phaseSpaceSize = phaseSpacePlot[s]->max - phaseSpacePlot[s]->min;
			DiscreteSpace plot_layout(phaseSpacePlot[s]->hDim,phaseSpacePlot[s]->vDim,1,phaseSpacePlot[s]->min,phaseSpaceSize,1);
			field.Initialize(plot_layout,owner);
			field.SetBoundaryConditions(xAxis,dirichletCell,dirichletCell);
			field.SetBoundaryConditions(yAxis,dirichletCell,dirichletCell);
			for (i=0;i<particle.size();i++)
			{
				position = owner->PositionFromPrimitive(particle[i].q);
				if (phaseSpacePlot[s]->theRgn->Inside(position,*owner))
				{
					phaseSpacePos.x = GetPhaseSpaceData(&particle[i],phaseSpacePlot[s]->hAxis);
					phaseSpacePos.y = GetPhaseSpaceData(&particle[i],phaseSpacePlot[s]->vAxis);
					phaseSpacePos.z = 0.5;
					if (phaseSpacePos.x>phaseSpacePlot[s]->min.x && phaseSpacePos.x<phaseSpacePlot[s]->max.x)
						if (phaseSpacePos.y>phaseSpacePlot[s]->min.y && phaseSpacePos.y<phaseSpacePlot[s]->max.y)
						{
							field.GetWeights(&weights,phaseSpacePos);
							field.InterpolateOnto( particle[i].number, weights );
							// phase space volume is divided out below
						}
				}
			}
			owner->strip[0].Sum(&field(0,0,0),&field(0,0,0),sizeof(tw::Float)*(phaseSpacePlot[s]->hDim+2)*(phaseSpacePlot[s]->vDim+2),master);
			if (owner->strip[0].Get_rank()==master)
			{
				field.ApplyBoundaryCondition();
				fileName.str("");
				fileName << phaseSpacePlot[s]->filename << ".dvdat";
				phaseSpacePlot[s]->outFile.open(fileName.str().c_str(),std::ios::binary | std::ios::app);
				for (j=1;j<=field.Dim(2);j++)
					for (i=1;i<=field.Dim(1);i++)
					{
						data = field(i,j,1);
						data /= phaseSpacePlot[s]->horVolFactor[i-1] * phaseSpacePlot[s]->verVolFactor[j-1];
						WriteBigEndian((char *)&data,sizeof(float),0,phaseSpacePlot[s]->outFile);
					}
				phaseSpacePlot[s]->outFile.close();
			}
		}
	}

	for (s=0;s<orbitDiagnostic.size();s++)
	{
		if (orbitDiagnostic[s]->WriteThisStep(owner->elapsedTime,owner->dt,owner->stepNow))
		{
			parData.clear();
			for (i=0;i<particle.size();i++)
			{
				position = owner->PositionFromPrimitive(particle[i].q);
				if (orbitDiagnostic[s]->theRgn->Inside(position,*owner))
					if (sqrt(1.0 + Norm(particle[i].p/restMass)) >= orbitDiagnostic[s]->minGamma)
					{
						owner->CurvilinearToCartesian(&position);
						parData.push_back(position.x);
						parData.push_back(particle[i].p.x);
						parData.push_back(position.y);
						parData.push_back(particle[i].p.y);
						parData.push_back(position.z - (owner->movingWindow ? owner->elapsedTime : 0.0));
						parData.push_back(particle[i].p.z);
						parData.push_back(particle[i].aux1);
						parData.push_back(particle[i].aux2);
					}
			}
			parBuffer.resize(parData.size());
			for (i=0;i<parData.size();i++)
				parBuffer[i] = parData[i];
			pts = parData.size();
			if (owner->strip[0].Get_rank()!=master)
			{
				owner->strip[0].Send(&pts,sizeof(tw::Int),master);
				owner->strip[0].Send(&parBuffer[0],sizeof(float)*pts,master);
			}
			if (owner->strip[0].Get_rank()==master)
			{
				fileName.str("");
				fileName << orbitDiagnostic[s]->filename << ".dvpar";
				orbitDiagnostic[s]->outFile.open(fileName.str().c_str(),std::ios::binary | std::ios::app);
				WriteBigEndian((char*)&parBuffer[0],sizeof(float)*pts,sizeof(float),orbitDiagnostic[s]->outFile);
				for (i=0;i<owner->strip[0].Get_size();i++)
				{
					if (i!=master)
					{
						owner->strip[0].Recv(&pts,sizeof(tw::Int),i);
						parBuffer.resize(pts);
						owner->strip[0].Recv(&parBuffer[0],sizeof(float)*pts,i);
						WriteBigEndian((char*)&parBuffer[0],sizeof(float)*pts,sizeof(float),orbitDiagnostic[s]->outFile);
					}
				}
				float zero = 0.0;
				for (i=0;i<8;i++)
					WriteBigEndian((char*)&zero,sizeof(float),sizeof(float),orbitDiagnostic[s]->outFile);
				orbitDiagnostic[s]->outFile.close();
			}
		}
	}

	for (s=0;s<detector.size();s++)
	{
		parData.clear();
		for (i=0;i<particle.size();i++)
		{
			position = owner->PositionFromPrimitive(particle[i].q);
			if (detector[s]->theRgn->Inside(position,*owner))
			{
				owner->CurvilinearToCartesian(&position);
				vel = particle[i].p/sqrt(1.0 + Norm(particle[i].p/restMass));
				lastPos = position - vel*dt;
				gamma = sqrt(1.0 + Norm(particle[i].p/restMass));
				if (lastPos.z < detector[s]->position.z && position.z > detector[s]->position.z)
					if (gamma >= detector[s]->minGamma)
					{
						parData.push_back(position.x);
						parData.push_back(particle[i].p.x/particle[i].p.z);
						parData.push_back(position.y);
						parData.push_back(particle[i].p.y/particle[i].p.z);
						parData.push_back(owner->elapsedTime); // zero order accurate
						parData.push_back(gamma-1.0);
						parData.push_back(particle[i].number);
					}
			}
		}
		parBuffer.resize(parData.size());
		for (i=0;i<parData.size();i++)
			parBuffer[i] = parData[i];
		pts = parData.size();
		if (owner->strip[0].Get_rank()!=master)
		{
			owner->strip[0].Send(&pts,sizeof(tw::Int),master);
			owner->strip[0].Send(&parBuffer[0],sizeof(float)*pts,master);
		}
		if (owner->strip[0].Get_rank()==master)
		{
			fileName.str("");
			fileName << detector[s]->filename << ".txt";
			detector[s]->outFile.open(fileName.str().c_str(),std::ios::binary | std::ios::app);
			detector[s]->WriteRecords(parBuffer);
			for (i=0;i<owner->strip[0].Get_size();i++)
			{
				if (i!=master)
				{
					owner->strip[0].Recv(&pts,sizeof(tw::Int),i);
					parBuffer.resize(pts);
					owner->strip[0].Recv(&parBuffer[0],sizeof(float)*pts,i);
					detector[s]->WriteRecords(parBuffer);
				}
			}
			detector[s]->outFile.close();
		}
	}
}
