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

void Mover::EncodingTest()
{
	tw::Int cell,ijk[4];
	cell = 0;
	space->DecodeCell(cell,ijk);
	cell = space->EncodeCell(ijk[0],ijk[1],ijk[2],ijk[3]);
	ASSERT_EQ(cell,0);
	for (tw::Int i=1;i<4;i++)
	{
		ASSERT_EQ(space->Ignorable(i),space->Dim(i)==1?true:false);
		ASSERT_EQ(ijk[i],space->Ignorable(i)?1:-1);
	}
}

void Mover::MinimizePrimitiveScalarTest()
{
	const float tolerance = 1e-5;
	Primitive q0;
	q0.cell = 1;
	q0.x[1] = 0.4;
	q0.x[2] = 0.6;
	q0.x[3] = 0.6;
	space->MinimizePrimitive(q0);
	ASSERT_EQ(q0.cell,2 + (1-space->Ignorable(2))*(2+space->Dim(2)));
	ASSERT_NEAR(q0.x[1] , 0.4 , tolerance);
	ASSERT_NEAR(q0.x[2] , -0.4 , tolerance);
	ASSERT_NEAR(q0.x[3] , -0.4 , tolerance);
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
		cell[i] = space->EncodeCell(1,1,1,1);
		for (tw::Int ax=0;ax<4;ax++)
			x[ax][i] = 0.0;
	}
	x[1][N-1] = 0.4;
	x[2][N-1] = 0.6;
	x[3][N-1] = -0.7;
	space->MinimizePrimitive(cell,ijk,x,domainMask);
	ASSERT_EQ(ijk[1][N-1] , 1);
	ASSERT_EQ(ijk[2][N-1] , space->Dim(2)==1 ? 1 : 2);
	ASSERT_EQ(ijk[3][N-1] , 0);
	ASSERT_NEAR(x[1][N-1] , 0.4 , tolerance);
	ASSERT_NEAR(x[2][N-1] , -0.4 , tolerance);
	ASSERT_NEAR(x[3][N-1] , 0.3 , tolerance);
}

void Mover::TranslationTest()
{
	InitTest();
	const tw::Float tolerance = 1e-5;
	// Now check the motion
	Primitive q0;
	tw::vec4 r0(0.5*timestep(*space),Corner(*space) + tw::vec3(0*dx(*space),.01*dy(*space),0*dz(*space)));
	tw::vec4 p0(sqrt(2),0,0,-1);
	tw::vec4 s0(1,0,0,0);
	tw::Float numDens=1.0;
	space->SetPrimitiveWithPosition(q0,r0);
	particle->push_back(Particle(numDens,q0,p0,s0));
	Advance();
	tw::vec4 x = space->PositionFromPrimitive((*particle)[0].q);
	ASSERT_NEAR(x[3] , r0[3]-timestep(*space)/p0[0] , tolerance);
	ASSERT_EQ((*particle)[0].number , 0.0);
	ASSERT_EQ(transfer->size() , 1);
	ASSERT_NEAR((*transfer)[0].number , numDens, tolerance);
	ASSERT_EQ((*transfer)[0].dst[1] , 0);
	ASSERT_EQ((*transfer)[0].dst[2] , 0);
	ASSERT_EQ((*transfer)[0].dst[3] , -1);
	ASSERT_EQ((*transfer)[0].ijk[1] , 1);
	ASSERT_EQ((*transfer)[0].ijk[2] , 1);
	ASSERT_EQ((*transfer)[0].ijk[3] , 0);
	ASSERT_NEAR((*transfer)[0].x[1] , (*particle)[0].q.x[1] , tolerance);
	ASSERT_NEAR((*transfer)[0].x[2] , (*particle)[0].q.x[2] , tolerance);
	ASSERT_NEAR((*transfer)[0].x[3] , (*particle)[0].q.x[3] , tolerance);
	CloseTest();
}

bool Mover::Test(tw::Int& id)
{
	if (id==1) {
		id++;
		EncodingTest();
		return true;
	} else if (id==2) {
		id++;
		MinimizePrimitiveScalarTest();
		return true;
	} else if (id==3) {
		id++;
		MinimizePrimitiveVectorTest();
		return true;
	} else if (id==4) {
		id++;
		TranslationTest();
		return true;
	} else if (id==5) {
		id++;
		UniformETest();
		return true;
	} else if (id==6) {
		id++;
		UniformBTest();
		return true;
	} else if (id==7) {
		id=0;
		PlaneWaveTest();
		return true;
	}
	return false;
}

void BorisMover::InitTest()
{
	Mover::InitTest();
	EM = new Field;
	sources = new Field;
	EM->Initialize(6,*space,task);
	sources->Initialize(4,*space,task);
}

void UnitaryMover::InitTest()
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

void PhotonMover::InitTest()
{
	Mover::InitTest();
	q0 = 0.0;
	m0 = 0.0;
	EM = new Field;
	sources = new Field;
	EM->Initialize(6,*space,task);
	sources->Initialize(4,*space,task);
}
