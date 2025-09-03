module;

#include "tw_includes.h"
#include "tw_test.h"

module elliptic;

const double tol = 1e-10;

void PoissonSolver::RegisterTests() {
    REGISTER(PoissonSolver,SheetChargeTestDirichlet);
    REGISTER(PoissonSolver,SheetChargeTestOpen);
    REGISTER(PoissonSolver,PointChargeTestDirichlet);
    REGISTER(PoissonSolver,PointChargeTestOpen);
}

void PoissonSolver::SheetChargeTestDirichlet() {
    // test solution for sheet of charge:
    // phi(k) = rho*dz^2*(3/2 - k/10 - |k-3|/2) 
    // holds in 4 cell system with sheet in cell 3, and ghost cells 0 and 5 at zero potential
    const tw::Float chargeElement = 1.0;
    x0 = x1 = y0 = y1 = tw::bc::fld::periodic;
    z0 = z1 = tw::bc::fld::dirichletCell;
    ScalarField phi,rho;
    phi.Initialize(*space,task);
    rho.Initialize(*space,task);
    if (task->strip[3].Get_rank()==1)
        for (tw::Int j=1;j<=space->Dim(2);j++)
            for (tw::Int i=1;i<=space->Dim(1);i++)
                rho(i,j,1) = chargeElement;
    Solve(phi,rho,-1.0);
    tw::Float rd2 = chargeElement*sqr(space->dx(3));
    if (task->strip[3].Get_rank()==0) {
        ASSERT_NEAR(phi(1,1,1),rd2*(1.5-0.1-1),tol);
        ASSERT_NEAR(phi(1,1,2),rd2*(1.5-0.2-0.5),tol);
    } else if (task->strip[3].Get_rank()==1) {
        // as of this writing test runner ignores this
        ASSERT_NEAR(phi(1,1,1),rd2*(1.5-0.3),tol);
        ASSERT_NEAR(phi(1,1,2),rd2*(1.5-0.4-0.5),tol);
    }
}

void PoissonSolver::SheetChargeTestOpen() {
    // test solution for sheet of charge
    // phi(i) = -rho*dz^2*|i| / 2
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
    tw::Float rd2 = chargeElement*sqr(space->dx(3));
    if (task->strip[3].Get_rank()==0) {
        ASSERT_NEAR(phi(1,1,1),-rd2,tol);
        ASSERT_NEAR(phi(1,1,2),-rd2/2,tol);
    } else if (task->strip[3].Get_rank()==1) {
        // as of this writing test runner ignores this
        ASSERT_NEAR(phi(1,1,1),0.0,tol);
        ASSERT_NEAR(phi(1,1,2),-rd2/2,tol);
    }
}

void PoissonSolver::PointChargeTestDirichlet() {
    // Test solution for point charge in transversely periodic system.
    // Rather than test an explicit solution we will do a local differential test.
    // Ghost cells are involved in the test both directly and indirectly.
    const tw::Float chargeElement = 1.0;
    x0 = x1 = y0 = y1 = tw::bc::fld::periodic;
    z0 = z1 = tw::bc::fld::dirichletCell;
    ScalarField phi,rho;
    phi.Initialize(*space,task);
    rho.Initialize(*space,task);
    if (task->strip[3].Get_rank()==1) {
        rho(1,1,1) = chargeElement;
    }
    Solve(phi,rho,-1.0);
    if (task->strip[3].Get_rank()==0) {
        tw::Float d2phi = 0.0;
        d2phi += (phi(0,1,1) - 2*phi(1,1,1) + phi(2,1,1)) / sqr(space->dx(1));
        d2phi += (phi(1,0,1) - 2*phi(1,1,1) + phi(1,2,1)) / sqr(space->dx(2));
        d2phi += (phi(1,1,0) - 2*phi(1,1,1) + phi(1,1,2)) / sqr(space->dx(3));
        ASSERT_NEAR(d2phi,0.0,tol);
        ASSERT_NEAR(phi(1,1,0),0,tol);
        ASSERT_NEAR(phi(2,1,0),0,tol);
        ASSERT_NEAR(phi(1,2,0),0,tol);
        ASSERT_NEAR(phi(2,2,0),0,tol);
    } else if (task->strip[3].Get_rank()==1) {
        // as of this writing test runner ignores this
        tw::Float d2phi = 0.0;
        d2phi += (phi(0,1,1) - 2*phi(1,1,1) + phi(2,1,1)) / sqr(space->dx(1));
        d2phi += (phi(1,0,1) - 2*phi(1,1,1) + phi(1,2,1)) / sqr(space->dx(2));
        d2phi += (phi(1,1,0) - 2*phi(1,1,1) + phi(1,1,2)) / sqr(space->dx(3));
        ASSERT_NEAR(d2phi,-chargeElement,tol);
    }
}

void PoissonSolver::PointChargeTestOpen() {
    // Test solution for point charge in transversely periodic system.
    // Rather than test an explicit solution we will do a local differential test.
    // Ghost cells are involved in the test indirectly.
    const tw::Float chargeElement = 1.0;
    x0 = x1 = y0 = y1 = tw::bc::fld::periodic;
    z0 = z1 = tw::bc::fld::natural;
    ScalarField phi,rho;
    phi.Initialize(*space,task);
    rho.Initialize(*space,task);
    if (task->strip[3].Get_rank()==1) {
        rho(1,1,1) = chargeElement;
    }
    Solve(phi,rho,-1.0);
    if (task->strip[3].Get_rank()==0) {
        tw::Float d2phi = 0.0;
        d2phi += (phi(0,1,1) - 2*phi(1,1,1) + phi(2,1,1)) / sqr(space->dx(1));
        d2phi += (phi(1,0,1) - 2*phi(1,1,1) + phi(1,2,1)) / sqr(space->dx(2));
        d2phi += (phi(1,1,0) - 2*phi(1,1,1) + phi(1,1,2)) / sqr(space->dx(3));
        ASSERT_NEAR(d2phi,0.0,1e-10);
    } else if (task->strip[3].Get_rank()==1) {
        // as of this writing test runner ignores this
        tw::Float d2phi = 0.0;
        d2phi += (phi(0,1,1) - 2*phi(1,1,1) + phi(2,1,1)) / sqr(space->dx(1));
        d2phi += (phi(1,0,1) - 2*phi(1,1,1) + phi(1,2,1)) / sqr(space->dx(2));
        d2phi += (phi(1,1,0) - 2*phi(1,1,1) + phi(1,1,2)) / sqr(space->dx(3));
        ASSERT_NEAR(d2phi,-chargeElement,tol);
    }
}