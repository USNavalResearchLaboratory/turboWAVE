#include "simulation.h"
#include "particles.h"
#include "fieldSolve.h"
#include "laserSolve.h"
using namespace tw::bc;

bool Species::Test(tw::Int& id)
{
	if (id==1) {
		id++;
		ReflectionTest();
		return true;
	} else if (id==2) {
		id=0;
		MoveWindowTest();
		return true;
	}
	return false;
}

void Species::ReflectionTest()
{
	laser = NULL;
	chi = NULL;
	qo_j4 = NULL;
	EM = new Field;
	sources = new Field;
	rho00 = new ScalarField;
	EM->Initialize(6,*this,owner);
	sources->Initialize(4,*this,owner);
	rho00->Initialize(*this,owner);
	Initialize();

	// Reflecting boundary condition test
	bc0[3] = par::reflecting;
	bc1[3] = par::reflecting;
	const tw::Float tolerance = 1e-5;
	Primitive q0;
	tw::vec3 dr(0,0,0.35*dz(*this));
	tw::vec3 r0(Corner(*this) - dr);
	tw::vec3 rtest(Corner(*this) + dr);
	tw::vec4 p0(sqrt(2),0,0,-1);
	tw::vec4 s0(1,0,0,0);
	SetPrimitiveWithPosition(q0,r0);
	AddParticle(1.0,q0,p0,s0);
	mover->AddTransferParticle(particle[0]);
	ApplyGlobalBoundaryConditions();
	if (owner->strip[3].Get_rank()==0)
	{
		ASSERT_NEAR(transfer[0].x[3], rtest.z, tolerance);
		ASSERT_NEAR(transfer[0].p[3], -p0[3], tolerance);
	}
	else
	{
		ASSERT_NEAR(transfer[0].x[3], r0.z, tolerance);
		ASSERT_NEAR(transfer[0].p[3], p0[3], tolerance);
	}

	delete EM;
	delete sources;
	delete rho00;
}

void Species::MoveWindowTest()
{
	REGISTER_TEST();
	laser = NULL;
	chi = NULL;
	qo_j4 = NULL;
	EM = new Field;
	sources = new Field;
	rho00 = new ScalarField;
	EM->Initialize(6,*this,owner);
	sources->Initialize(4,*this,owner);
	rho00->Initialize(*this,owner);
	Initialize();

	// Begin move window test
	const tw::Float tolerance = 1e-5;
	particle.clear();
	transfer.clear();
	Primitive q0;
	tw::vec3 r0;
	tw::vec4 p0(1,0,0,0),s0(1,0,0,0);
	if (owner->strip[3].Get_rank()==1)
	{
		r0 = tw::vec3(0.0,0.0,Corner(*this).z);
		SetPrimitiveWithPosition(q0,r0);
		AddParticle(1.0,q0,p0,s0);
		BeginMoveWindow();
		ASSERT_EQ(particle.size(),0);
		ASSERT_EQ(transfer.size(),1);
		ASSERT_NEAR(transfer[0].x[3], r0.z - dz(*this), tolerance);
		ASSERT_EQ(transfer[0].dst[0],1);
		ASSERT_EQ(transfer[0].dst[1],0);
		ASSERT_EQ(transfer[0].dst[2],0);
		ASSERT_EQ(transfer[0].dst[3],-1);
	}

	// Simulate message passing the transfer particle

	// Finish move window test

	delete EM;
	delete sources;
	delete rho00;
}
