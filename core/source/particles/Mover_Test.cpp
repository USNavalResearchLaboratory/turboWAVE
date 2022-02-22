#include "meta_base.h"
#include "computeTool.h"
#include "bundle.h"
#include "pusher.h"
#include "slicer.h"
#include "tiler.h"
#include "mover.h"

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
	// TODO: provide the data to Simulation which must resize the grid before creating the tool
	task->Initialize(tw::idx4(1,1,2).array,tw::idx4(4,1,4).array,tw::idx4(1,1,0).array);
	space->Resize(*task,tw::vec3(0,0,0),tw::vec3(0.8,0.2,0.8),2);
	space->SetupTimeInfo(0.1);
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

void Mover::MinimizePrimitiveScalarTest()
{
	const tw::Float tolerance = 1e-5;
	Primitive q0;
	q0.cell = 1;
	q0.x[0] = 0.0;
	q0.x[1] = 0.0;
	q0.x[2] = 0.6;
	space->MinimizePrimitive(q0);
	TW_MPI_Lock();
	assert(q0.cell == 2);
	assert(q0.x[0] == 0.0);
	assert(q0.x[1] == 0.0);
	assert(fabs(q0.x[2] + 0.4) < tolerance);
	TW_MPI_Unlock();
	MPI_Barrier(MPI_COMM_WORLD);
}

void Mover::MinimizePrimitiveVectorTest()
{
	const tw::Float tolerance = 1e-5;
	const tw::Int AB = tw::vec_align_bytes;
	const tw::Int N = tw::max_bundle_size;
	alignas(AB) tw::Int cell[N];
	alignas(AB) tw::Int ijk[4][N];
	alignas(AB) float x[4][N];
	alignas(AB) float domainMask[N];
	for (tw::Int i=0;i<N;i++)
	{
		cell[i] = space->EncodeCell(1,1,1);
		for (tw::Int ax=0;ax<4;ax++)
			x[ax][i] = 0.0;
	}
	x[1][N-1] = 0.4;
	x[2][N-1] = 0.6;
	x[3][N-1] = -0.7;
	space->MinimizePrimitive(cell,ijk,x,domainMask);
	TW_MPI_Lock();
	assert(ijk[1][N-1] == 1);
	assert(ijk[2][N-1] == 1);
	assert(ijk[3][N-1] == 0);
	assert(fabs(x[1][N-1] - 0.4) < tolerance);
	assert(fabs(x[2][N-1] + 0.4) < tolerance);
	assert(fabs(x[3][N-1] - 0.3) < tolerance);
	TW_MPI_Unlock();
	MPI_Barrier(MPI_COMM_WORLD);
}

void Mover::TranslationTest()
{
	const tw::Float tolerance = 1e-5;
	// At present any type of mover will need the following resources
	assert(EM!=NULL);
	assert(sources!=NULL);
	assert(particle!=NULL);
	assert(transfer!=NULL);
	// Now check the motion
	Primitive q0;
	tw::vec3 r0(Corner(*space) + tw::vec3(0*dx(*space),.01*dy(*space),0*dz(*space)));
	tw::vec3 p0(0,0,-1);
	tw::Float numDens=1.0,aux1=1.0,aux2=1.0;
	space->SetPrimitiveWithPosition(q0,r0);
	particle->push_back(Particle(p0,q0,numDens,aux1,aux2));
	Advance();
	tw::vec3 r = space->PositionFromPrimitive((*particle)[0].q);
	TW_MPI_Lock();
	// std::cout << (*particle)[0].p << std::endl;
	// std::cout << space->Ignorable(1) << " , " << space->Ignorable(2) << std::endl;
	// std::cout << q0 << " , " << (*particle)[0].q << std::endl;
	assert(fabs(r.z-(r0.z-timestep(*space)/sqrt(1+Norm(p0)))) < tolerance);
	assert((*particle)[0].number == 0.0);
	assert(transfer->size() == 1);
	assert((*transfer)[0].number == numDens);
	assert((*transfer)[0].dst[1] == 0);
	assert((*transfer)[0].dst[2] == 0);
	assert((*transfer)[0].dst[3] == -1);
	assert(fabs((*transfer)[0].x[3]-r.z) < tolerance);
	TW_MPI_Unlock();
	MPI_Barrier(MPI_COMM_WORLD);
}

bool Mover::Test()
{
	InitTest();
	MinimizePrimitiveScalarTest();
	MinimizePrimitiveVectorTest();
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
