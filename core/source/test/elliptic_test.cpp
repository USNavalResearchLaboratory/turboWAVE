module;

#include "tw_includes.h"
#include "tw_test.h"

module elliptic;

void PoissonSolver::Test()
{
    // test solution for sheet of charge
    // phi(i) = -rho*dz*|i*dz| / 2
    // where i measures cell displacement from the sheet
    const tw::Float chargeElement = 1.0;
    x0 = x1 = y0 = y1 = tw::bc::fld::periodic;
    z0 = z1 = tw::bc::fld::natural;
    ScalarField phi,rho;
    phi.Initialize(*space,task);
    rho.Initialize(*space,task);
    if (task->strip[3].Get_rank()==1)
        for (tw::Int j=1;j<=space->Dim(2);j++)
            for (tw::Int i=1;i<=space->Dim(1);i++)
                rho(i,j,1) = chargeElement;
    Solve(phi,rho,-1.0);
    // this version of test runner assertions have to come from rank 0
    if (task->strip[3].Get_rank()==0)
    {
        ASSERT_NEAR(phi(1,1,2),-0.5*chargeElement*sqr(space->dx(3)),1e-6);
        ASSERT_NEAR(phi(1,1,1),-chargeElement*sqr(space->dx(3)),1e-6);
    }
}
