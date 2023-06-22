#include "meta_base.h"
#include "computeTool.h"
#include "bundle.h"
#include "pusher.h"
#include "slicer.h"
#include "tiler.h"
#include "mover.h"

void Mover::UniformETest()
{
	InitTest();
	const tw::Float tolerance = 1e-3;
	// Create uniform electric field
	for (auto cell : CellRange(*space,true))
		(*EM)(cell,2) = 1.0;
	Primitive q;
	tw::vec3 r0(Corner(*space) + tw::vec3(dx(*space)/2,dy(*space)/2,dz(*space)/2));
	tw::vec3 p0(0,0,0);
	tw::Float numDens=1.0,aux1=1.0,aux2=1.0;
	space->SetPrimitiveWithPosition(q,r0);
	particle->push_back(Particle(p0,q,numDens,aux1,aux2));
	Advance();
	tw::vec3 p = (*particle)[0].p;
	p0 += q0*timestep(*space)/m0;
	ASSERT_NEAR(p.z , p0.z , tolerance);
	CloseTest();
}

void Mover::UniformBTest()
{
	InitTest();
	const tw::Float tolerance = 1e-3;
	// Create uniform magnetic field
	for (auto cell : CellRange(*space,true))
		(*EM)(cell,4) = 1.0;
	Primitive q;
	tw::vec3 r0(Corner(*space) + tw::vec3(dx(*space)/2,dy(*space)/2,dz(*space)/2));
	tw::vec3 p0(0,0,1);
	tw::Float numDens=1.0,aux1=1.0,aux2=1.0;
	space->SetPrimitiveWithPosition(q,r0);
	particle->push_back(Particle(p0,q,numDens,aux1,aux2));
	Advance();
	tw::vec3 p = (*particle)[0].p;
	ASSERT_NEAR(Norm(p), Norm(p0), tolerance);
	tw::Float theta = Magnitude(p | p0);
	tw::Float theta_expected = sqrt(0.5)*timestep(*space);
	ASSERT_NEAR(theta, theta_expected, tolerance);
	CloseTest();
}
