module;

#include "tw_includes.h"
#include "tw_test.h"

module eos;

std::tuple<Field,Field,ScalarField,ScalarField,ScalarField> EOSIdealGas::InitTest(tw::dnum n0, tw::dnum T0)
{
    // We are intializing a hydro state without any help from higher level SPARC modules,
    // therefore there is a lot of work to do.
    auto N2 = sparc::material();
    auto gamma = 1.4;
    auto NAvogadro = 6.02214e23;
    auto cp = 28.99 / NAvogadro; // J/mol/K/NA = J/particle/K
    N2.mass = 28*1836;
    N2.cvm = cp / gamma / mks::kB;
    auto u0 = (cp/gamma) * (n0 >> mks) * (T0 >> mks) * tw::dims::energy_density >> mks;

    // create fields and setup indexing
    auto hidx = sparc::hydro_set();
    auto eidx = sparc::eos_set();
    auto IE = ScalarField();
    auto nm = ScalarField();
    auto nu_e = ScalarField();
    auto hydro = Field();
    auto eos = Field();
    hidx.Load(0,1);
    eidx.Load(0);
    IE.Initialize(*space,task);
    nm.Initialize(*space,task);
    nu_e.Initialize(*space,task);
    hydro.Initialize(hidx.count+1,*space,task);
    eos.Initialize(eidx.count,*space,task);
    
    Initialize();
    SetupIndexing(0,hidx,eidx,N2);
    nu_e = 1.0;
    for (auto cell : EntireCellRange(*space,1)) {
        hydro(cell,hidx.ni) = n0 >> native;
        hydro(cell,hidx.u) = u0 >> native;
        IE(cell) = hydro(cell,hidx.u);
        nm(cell) = hydro(cell,hidx.ni)*N2.mass;
        eos(cell,eidx.T) = IE(cell) / N2.cvm / hydro(cell,hidx.ni);
    }
    return std::tuple<Field,Field,ScalarField,ScalarField,ScalarField>(hydro,eos,IE,nm,nu_e);
}

void EOSIdealGas::PressureTest()
{
    auto n0 = tw::dnum("2.5e19 [/cm3]");
    auto T0 = tw::dnum("300 [K]");
    auto [hydro,eos,IE,nm,nu_e] = InitTest(n0, T0);
    AddPKV(IE,nm,nu_e,hydro,eos);
    if (task->strip[0].Get_rank()==0)
    {
        const tw::Float pressure_mks = eos(1,1,1,1,eidx.P) * tw::dims::pressure >> native >> mks;
        const tw::Float expected = (n0 >> mks) * mks::kB * (T0 >> mks);
        ASSERT_NEAR(pressure_mks , expected , expected/1e8);
    }
}

std::tuple<Field,Field,ScalarField,ScalarField,ScalarField> EOSTillotson::InitTest()
{
    auto al = sparc::material();
    al.mass = 26.98*1836;

    // Tillotson parameters for Aluminum, from original report
    rho0 = tw::dnum("2.7 [g/cm3]") >> native;
    a = 0.5;   // Tillotson Coefficient
    b = 1.63;   // Tillotson Coefficient
    A = tw::dnum("0.752e6 [bar]") >> native;   // Tillotson Coefficient, pressure
    B = tw::dnum("0.65e6 [bar]") >> native;   // Tillotson Coefficient, pressure
    alpha = 5.0;   // Tillotson Coefficient
    beta = 5.0;   // Tillotson Coefficient
    rhoIV = rho0/1.1;   // Incipient vaporization density
    E0 = tw::dnum("0.05e12 [ergs/g]") >> native;   // Reference energy
    EIV = tw::dnum("0.0276e12 [ergs/g]") >> native;   // Incipient vaporization specific energy
    ECV = tw::dnum("0.141e12 [ergs/g]") >> native;   // Complete vaporization specific energy

    // create fields and setup indexing
    auto hidx = sparc::hydro_set();
    auto eidx = sparc::eos_set();
    auto IE = ScalarField();
    auto nm = ScalarField();
    auto nu_e = ScalarField();
    auto hydro = Field();
    auto eos = Field();
    hidx.Load(0,1);
    eidx.Load(0);
    IE.Initialize(*space,task);
    nm.Initialize(*space,task);
    nu_e.Initialize(*space,task);
    hydro.Initialize(hidx.count+1,*space,task);
    eos.Initialize(eidx.count,*space,task);
    
    Initialize();
    SetupIndexing(0,hidx,eidx,al);
    // This EOS does not worry about partial pressures and therefore does not use IE and nm.
    // The expected data is trivially derived from the Hugoniot table given in Tillotson's report.
    nu_e = 1.0;
    for (auto cell : EntireCellRange(*space,1)) {
        hydro(cell,hidx.ni) = tw::dnum("2.659e23 [/cm3]") >> native;
        nm(cell) = hydro(cell,hidx.ni)*al.mass;
        hydro(cell,hidx.u) = tw::dnum("3.597e7 [J/cm3]") >> native;
        IE(cell) = hydro(cell,hidx.u);
    }
    return std::tuple<Field,Field,ScalarField,ScalarField,ScalarField>(hydro,eos,IE,nm,nu_e);
}

void EOSTillotson::PressureTest()
{
    auto [hydro,eos,IE,nm,nu_e] = InitTest();
    AddPKV(IE,nm,nu_e,hydro,eos);
    if (task->strip[0].Get_rank()==0)
    {
        const tw::Float pressure_cgs = eos(1,1,1,1,eidx.P) * tw::dims::pressure >> native >> cgs;
        const tw::Float expected = 208.83e12; // dynes/cm2
        ASSERT_NEAR(pressure_cgs , expected , expected/100);
    }
}
