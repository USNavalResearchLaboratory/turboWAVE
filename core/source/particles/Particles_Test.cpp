// z-cell arrangement for these tests:
//  || cell 0 | cell 1 || cell 2 | cell 3 ||
//  || node 0 | node 0 || node 1 | node 1 ||
//  0.0       0.2      0.4       0.6      0.8

module;

#include "tw_includes.h"
#include <tw_test.h>

module particles;

using namespace tw::bc;

bool Species::Test(tw::Int& id)
{
	if (id==1) {
		id++;
		EncodingTest();
		return true;
	} else if (id==2) {
		id++;
		ReflectionTest();
		return true;
	} else if (id==3) {
		id=0;
		MoveWindowTest();
		return true;
	}
	return false;
}

void Species::EncodingTest()
{
	tw::Int ijk[4];
	tw::Int cell = EncodeCell(1,2,1,1);
	DecodeCell(cell,ijk);
	ASSERT_EQ(ijk[0],1);
	ASSERT_EQ(ijk[1],2);
	ASSERT_EQ(ijk[2],1);
	ASSERT_EQ(ijk[3],1);
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
	tw::vec4 dr(0.5*this->dx(0),0,0,0.35*this->dx(3));
	tw::vec4 r0 = this->Corner() - dr;
	tw::vec4 rexpect = owner->strip[3].Get_rank() == 0 ? this->Corner() + dr : r0;
	tw::vec4 p0(sqrt(2),0,0,-1);
	tw::vec4 s0(1,0,0,0);
	SetPrimitiveWithPosition(q0,r0);
	ASSERT_NEAR(q0.x[0], 0.0, tolerance);
	ASSERT_NEAR(q0.x[3], 0.15, tolerance);
	AddParticle(1.0,q0,p0,s0,0.0);
	mover->AddTransferParticle(particle[0]);
	ASSERT_NEAR(transfer[0].x[0], q0.x[0], tolerance);
	ASSERT_NEAR(transfer[0].x[3], q0.x[3], tolerance);
	ApplyGlobalBoundaryConditions();
	if (owner->strip[3].Get_rank()==0)
	{
		ASSERT_EQ(transfer[0].ijk[3], 1);
		ASSERT_NEAR(transfer[0].x[3], -q0.x[3], tolerance);
		ASSERT_NEAR(transfer[0].p[3], -p0[3], tolerance);
	}
	else
	{
		ASSERT_EQ(transfer[0].ijk[3], 0);
		ASSERT_NEAR(transfer[0].x[3], q0.x[3], tolerance);
		ASSERT_NEAR(transfer[0].p[3], p0[3], tolerance);
	}
	// as an extra check, make sure we recover expected global coordinates
	tw::Int cell = EncodeCell(transfer[0].ijk[0],transfer[0].ijk[1],transfer[0].ijk[2],transfer[0].ijk[3]);
	Primitive qf(cell,transfer[0].x[0],transfer[0].x[1],transfer[0].x[2],transfer[0].x[3]);
	tw::vec4 ractual = PositionFromPrimitive(qf);
	ASSERT_NEAR(ractual[3], rexpect[3], tolerance);
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
	tw::vec4 r0;
	tw::vec4 p0(1,0,0,0),s0(1,0,0,0);
	if (owner->strip[3].Get_rank()==1)
	{
		r0 = tw::vec4(-0.5*this->dx(0),0,0,this->Corner()[3]);
		SetPrimitiveWithPosition(q0,r0);
		AddParticle(1.0,q0,p0,s0,0.0);
		BeginMoveWindow();
		ASSERT_EQ(particle.size(),0);
		ASSERT_EQ(transfer.size(),1);
		ASSERT_NEAR(transfer[0].x[3], -0.5, tolerance);
		ASSERT_EQ(transfer[0].ijk[3],0);
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
