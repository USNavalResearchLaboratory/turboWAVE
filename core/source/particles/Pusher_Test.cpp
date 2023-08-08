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
    tw::Float qmRatio = m0==0.0 ? 0.0 : q0/m0;
    tw::vec4 r0(-0.5*timestep(*space),Corner(*space) + tw::vec3(dx(*space)/2,dy(*space)/2,dz(*space)/2));
    tw::vec4 p0 = m0==0.0 ? tw::vec4(1,0,0,1) : tw::vec4(1,0,0,0);
    tw::vec4 s0(1,0,0,0);
    tw::Float numDens=1.0;
    space->SetPrimitiveWithPosition(q,r0);
    particle->push_back(Particle(numDens,q,p0,s0));
    Advance();
    tw::vec4 p = (*particle)[0].p;
    p0[3] += qmRatio*timestep(*space);
    ASSERT_NEAR(p[3] , p0[3] , tolerance);
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
    tw::vec4 r0(-0.5*timestep(*space),Corner(*space) + tw::vec3(dx(*space)/2,dy(*space)/2,dz(*space)/2));
    tw::vec4 p0 = m0==0.0 ? tw::vec4(1,0,0,1) : tw::vec4(sqrt(2),0,0,1);
    tw::vec4 s0(1,0,0,0);
    tw::Float numDens=1.0;
    space->SetPrimitiveWithPosition(q,r0);
    particle->push_back(Particle(numDens,q,p0,s0));
    Advance();
    tw::vec4 p = (*particle)[0].p;
    ASSERT_NEAR(Norm(p.spatial()), Norm(p0.spatial()), tolerance);
    tw::Float theta = Magnitude(p.spatial() | p0.spatial());
    tw::Float theta_expected = m0==0.0 ? 0.0 : sqrt(0.5)*timestep(*space);
    ASSERT_NEAR(theta, theta_expected, tolerance);
    CloseTest();
}

void Mover::PlaneWaveTest()
{
    InitTest();
    const tw::Float tolerance = 1e-3;
    const tw::Int steps = 2;
    const tw::Float dr = space->dX(0,3);
    const tw::Float dt = 0.25*dr; // should check courant
    const tw::Float k0 = 0.1/dr;
    const tw::Float w0 = 0.1/dt;
    Primitive q;
    tw::vec4 r0(-0.5*timestep(*space),Corner(*space) + tw::vec3(dx(*space)/2,dy(*space)/2,dz(*space)/2));
    tw::vec4 p0 = m0==0.0 ? tw::vec4(1,0,0,1) : tw::vec4(sqrt(2),0,0,1);
    tw::vec4 s0(1,0,0,0);
    tw::Float numDens=1.0;
    space->SetPrimitiveWithPosition(q,r0);
    particle->push_back(Particle(numDens,q,p0,s0));
    for (tw::Int i=0;i<steps;i++)
    {
        // Update plane wave grid analytically
        for (auto cell : CellRange(*space,true))
        {
            (*EM)(cell,0) = 1*sin(w0*dt*i - k0*space->Pos(cell).z);
            (*EM)(cell,4) = 1*sin(w0*dt*i - k0*(space->Pos(cell).z - dr/2));
        }
        Advance();
    }
    // check plane wave invariant
    tw::vec4 p = (*particle)[0].p;
    ASSERT_NEAR(p[0]-p[3],p0[0]-p0[3],tolerance);
    CloseTest();
}