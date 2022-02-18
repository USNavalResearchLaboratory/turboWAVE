#include "../simulation.h"
#include "particles.h"
#include "../solver/fieldSolve.h"
#include "../solver/laserSolve.h"
using namespace tw::bc;

bool Species::Test()
{
	owner->Initialize(tw::idx4(1,1,2).array,tw::idx4(4,4,4).array,tw::idx4(1,1,0).array);
	owner->Resize(*owner,tw::vec3(0,0,0),tw::vec3(0.8,0.8,0.8),2);
	owner->UpdateTimestep(0.1);
	EM = new Field;
	sources = new Field;
	rho00 = new ScalarField;
	EM->Initialize(6,*this,owner);
	sources->Initialize(4,*this,owner);
	rho00->Initialize(*this,owner);
	mover = (Mover*)owner->CreateTool("test_tool",tw::tool_type::borisMover);
	Initialize();

	// Reflecting boundary condition test
	bc0[3] = par::reflecting;
	bc1[3] = par::reflecting;
	const tw::Float tolerance = 1e-5;
	Primitive q0;
	tw::vec3 dr(0,0,0.35*dz(*this));
	tw::vec3 r0(Corner(*this) - dr);
	tw::vec3 rtest(Corner(*this) + dr);
	tw::vec3 p0(0,0,-1);
	SetPrimitiveWithPosition(q0,r0);
	AddParticle(p0,q0,1.0);
	mover->AddTransferParticle(particle[0]);
	ApplyGlobalBoundaryConditions();
	if (owner->strip[3].Get_rank()==0)
	{
		assert(fabs(rtest.z-transfer[0].x[3])<tolerance);
		assert(transfer[0].p[3]==-p0.z);
	}
	else
	{
		assert(fabs(r0.z-transfer[0].x[3])<tolerance);
		assert(transfer[0].p[3]==p0.z);
	}

	delete EM;
	delete sources;
	delete rho00;
	// mover is deleted by Species

	return true;
}
