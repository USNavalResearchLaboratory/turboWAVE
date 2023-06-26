#include "meta_base.h"

bool MetricSpace::Test(tw::Int& id)
{
    id = 0;
    REGISTER_TEST();
    tw::Float pos = X(1,1);
    ASSERT_NEAR(pos, corner.x + spacing.x/2, 1e-4);
    return true;
}
