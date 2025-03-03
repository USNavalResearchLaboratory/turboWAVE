module;

#include "tw_includes.h"
#include "tw_test.h"

module injection;

bool HermiteGauss::Test(tw::Int& id)
{
    // test simple Gaussian beam at a few spacetime points
    tw::vec3 a1;
    REGISTER_TEST();
    id = 0;
    direction = tw::vec3(0.0,0.0,1.0);
    focusPosition = tw::vec3(0.0,0.0,0.0);
    a = tw::vec3(1.0,0.0,0.0);
    w = 1.0;
    pulseShape.risetime = 1.0;
    pulseShape.holdtime = 0.0;
    pulseShape.falltime = 1.0;
    modeData.scale[0] = 1.0;
    modeData.scale[1] = 1.0;
    pulseShape.whichProfile = tw::profile::shape::quintic;
    Initialize();

    if (task->strip[3].Get_rank()==0)
    {
        a1 = VectorPotential(1.0,tw::vec3(0.0,0.0,0.0));
        ASSERT_NEAR(a1.x,1.0,1e-4);
        ASSERT_NEAR(a1.y,0.0,1e-4);
        ASSERT_NEAR(a1.z,0.0,1e-4);

        a1 = VectorPotential(1.0,tw::vec3(0.0,1.0,0.0));
        ASSERT_NEAR(a1.x,1.0/exp(1),1e-4);
        ASSERT_NEAR(a1.y,0.0,1e-4);
        ASSERT_NEAR(a1.z,0.0,1e-4);

        a1 = VectorPotential(0.5,tw::vec3(0.0,0.0,0.0));
        ASSERT_NEAR(a1.x,0.5*cos(w*0.5),1e-4);
        ASSERT_NEAR(a1.y,0.0,1e-4);
        ASSERT_NEAR(a1.z,0.0,1e-4);
    }

    return true;
}
