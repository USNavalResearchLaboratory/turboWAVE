module;

#include "tw_includes.h"
#include <tw_test.h>

module fluid;
import driver;

using namespace tw::bc;

void Fluid::AdvectionTest()
{
    auto bkg = std::make_shared<UniformProfile>("hello",space,task);
    bkg->temperature = .01;
    bkg->density = .01;
    bkg->Initialize();
    profiles.push_back(bkg);
	Initialize();

    // setup a step function and velocity profile
    // important to have flux = 0 at boundaries
    for (auto strip : StripRange(*space,3,0,1,strongbool::no))
    {
        for (tw::Int k=1;k<=Dim(3);k++)
        {
            const tw::Int kg = space->GlobalCellIndex(k,3);
            state1(strip,k,0) *= kg <= space->GlobalDim(3)/2 ? 1.0 : 2.0;
            vel(strip,k,3) = kg > 1 && kg < space->GlobalDim(3) ? 0.9 : 0.0;
        }
    }

    // advection
	state0 = state1; // state1 is the predicted state used to advance state0
	FCT_Driver convector(&state0,&state1,&vel,NULL,space);
	convector.SetDensityElements(Rng(0));
	convector.SetVelocityElement(3);
	convector.Convect(tw::grid::z,fld::dirichletCell,fld::dirichletCell,0.5*dx(0));

    // check positivity
    tw::Float pos = 1.0;
    for (auto cell : InteriorCellRange(*space,1))
        pos *= state0(cell,0)>0.0 ? 1.0 : 0.0;

    ASSERT_NEAR(pos,1.0,1e-6);
}

void Fluid::ConservationTest()
{
    auto bkg = std::make_shared<UniformProfile>("hello",space,task);
    bkg->temperature = .01;
    bkg->density = .01;
    bkg->Initialize();
    profiles.push_back(bkg);
	Initialize();

    // setup a step function and velocity profile
    // important to have flux = 0 at boundaries
    for (auto strip : StripRange(*space,3,0,1,strongbool::no))
    {
        for (tw::Int k=1;k<=Dim(3);k++)
        {
            const tw::Int kg = space->GlobalCellIndex(k,3);
            state1(strip,k,0) *= kg <= space->GlobalDim(3)/2 ? 1.0 : 2.0;
            vel(strip,k,3) = kg > 1 && kg < space->GlobalDim(3) ? 0.9 : 0.0;
        }
    }

    // initial mass for conservation check
    tw::Float initialMass = 0.0;
    for (auto cell : InteriorCellRange(*space,1))
        initialMass += state1(cell,0);
    task->strip[0].AllSum(&initialMass,&initialMass,sizeof(tw::Float),0);

    // advection
	state0 = state1; // state1 is the predicted state used to advance state0
	FCT_Driver convector(&state0,&state1,&vel,NULL,space);
	convector.SetDensityElements(Rng(0));
	convector.SetVelocityElement(3);
	convector.Convect(tw::grid::z,fld::dirichletCell,fld::dirichletCell,0.5*dx(0));

    // check conservation
    tw::Float finalMass = 0.0;
    for (auto cell : InteriorCellRange(*space,1))
        finalMass += state0(cell,0);
     task->strip[0].AllSum(&finalMass,&finalMass,sizeof(tw::Float),0);

    ASSERT_NEAR(finalMass,initialMass,initialMass/1e6);
}
