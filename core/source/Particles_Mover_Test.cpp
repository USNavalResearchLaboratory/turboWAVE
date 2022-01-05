#include "meta_base.h"
#include "computeTool.h"
#include "particles_bundle.h"
#include "particles_pusher.h"
#include "particles_slicer.h"
#include "particles_tiler.h"
#include "particles_mover.h"

void Mover::InitTest()
{
	q0 = -1.0;
	m0 = 1.0;
	EM = NULL;
	sources = NULL;
	laser = NULL;
	chi = NULL;
	qo_j4 = NULL;
	ESField = NULL;
	particle = new std::vector<Particle>;
	transfer = new std::vector<TransferParticle>;
}

void Mover::CloseTest()
{
	if (EM!=NULL)
		delete EM;
	if (sources!=NULL)
		delete sources;
	if (laser!=NULL)
		delete laser;
	if (chi!=NULL)
		delete chi;
	if (particle!=NULL)
		delete particle;
	if (transfer!=NULL)
		delete transfer;
}

void Mover::TranslationTest()
{
	// At present any type of mover will need the following resources
	assert(EM!=NULL);
	assert(sources!=NULL);
	assert(particle!=NULL);
	assert(transfer!=NULL);
	// Now check the motion
	const tw::Float tolerance = 1e-5;
	Primitive q0;
	tw::vec3 r0(Corner(*space));// - tw::vec3(2*dx(*space),2*dy(*space),2*dz(*space)));
	tw::vec3 p0(0,0,-1);
	tw::Float numDens=1.0,aux1=1.0,aux2=1.0;
	space->SetPrimitiveWithPosition(q0,r0);
	particle->push_back(Particle(p0,q0,numDens,aux1,aux2));
	Advance();
	tw::vec3 r = space->PositionFromPrimitive((*particle)[0].q);
	assert(fabs(r.z-(r0.z-timestep(*space)/sqrt(1+Norm(p0)))) < tolerance);
	assert((*particle)[0].number == 0.0);
	assert(transfer->size() == 1);
	assert((*transfer)[0].number == numDens);
	assert((*transfer)[0].dst[1] == 0);
	assert((*transfer)[0].dst[2] == 0);
	assert((*transfer)[0].dst[3] == -1);
	assert(fabs((*transfer)[0].x[3]-r.z) < tolerance);
}

bool Mover::Test()
{
	InitTest();
	TranslationTest();
	CloseTest();
	return true;
}

void BorisMover::InitTest()
{
	Mover::InitTest();
	EM = new Field;
	sources = new Field;
	EM->Initialize(6,*space,task);
	sources->Initialize(4,*space,task);
}

void PGCMover::InitTest()
{
	Mover::InitTest();
	EM = new Field;
	sources = new Field;
	laser = new Field;
	chi = new Field;
	EM->Initialize(6,*space,task);
	sources->Initialize(4,*space,task);
	laser->Initialize(8,*space,task);
	chi->Initialize(2,*space,task);
}
