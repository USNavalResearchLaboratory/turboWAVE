module;

#include "tw_includes.h"

export module factory;
import base;
import driver;
import fields;
import injection;
import parabolic;
import hyperbolic;
import elliptic;
import physics;
import diagnostics;
import solid_state;
import qstate;
import qed;
import mover;
import electrostatic;
import field_solve;
import laser_solve;
import solid_state;
import fluid;
import particles;
import quantum;
import input;
import metric_space_test;
import iterator_test;
import fft_test;

export namespace factory {

/// @brief verify the object key from the preamble and display help if necessary
/// @param raw_key key directly from preamble, white space will be normalized before comparing
/// @param similar_keys comma delimited list of similar keys for display in case of failure
/// @return true if it is known, false if not
/// @throws error if there are multiple matches
bool VerifyKey(const std::string& raw_key,std::string& similar_keys)
{
	std::map<std::string,tw::tool_type> tool_map = ComputeTool::Map();
	
	std::string normalized,word;
	std::stringstream raw_stream(raw_key);
	do {
		raw_stream >> word;
		if (!normalized.empty()) {
			normalized += " ";
		}
		normalized += word;
	} while (!raw_stream.eof());

	if (tool_map.find(normalized) != tool_map.end()) {
		return true;
	} else {
		for (const auto& pair : tool_map) {
			tw::input::BuildSimilar(similar_keys,normalized,pair.first);
		}
		return false;
	}
}

SharedTool CreateToolFromType(const std::string& name,tw::tool_type theType,MetricSpace *ms,Task *tsk)
{
	SharedTool ans;
	switch (theType)
	{
		case tw::tool_type::warp:
			ans = std::make_shared<Warp>(name,ms,tsk);
			break;

		case tw::tool_type::entireRegion:
			ans = std::make_shared<EntireRegion>(name,ms,tsk);
			break;
		case tw::tool_type::rectRegion:
			ans = std::make_shared<RectRegion>(name,ms,tsk);
			break;
		case tw::tool_type::prismRegion:
			ans = std::make_shared<PrismRegion>(name,ms,tsk);
			break;
		case tw::tool_type::ellipsoidRegion:
			ans = std::make_shared<EllipsoidRegion>(name,ms,tsk);
			break;
		case tw::tool_type::circRegion:
			ans = std::make_shared<CircRegion>(name,ms,tsk);
			break;
		case tw::tool_type::cylinderRegion:
			ans = std::make_shared<CylinderRegion>(name,ms,tsk);
			break;
		case tw::tool_type::roundedCylinderRegion:
			ans = std::make_shared<RoundedCylinderRegion>(name,ms,tsk);
			break;
		case tw::tool_type::cylindricalShellRegion:
			ans = std::make_shared<CylindricalShellRegion>(name,ms,tsk);
			break;
		case tw::tool_type::trueSphereRegion:
			ans = std::make_shared<TrueSphere>(name,ms,tsk);
			break;
		case tw::tool_type::boxArrayRegion:
			ans = std::make_shared<BoxArrayRegion>(name,ms,tsk);
			break;
		case tw::tool_type::torusRegion:
			ans = std::make_shared<TorusRegion>(name,ms,tsk);
			break;
		case tw::tool_type::coneRegion:
			ans = std::make_shared<ConeRegion>(name,ms,tsk);
			break;
		case tw::tool_type::tangentOgiveRegion:
			ans = std::make_shared<TangentOgiveRegion>(name,ms,tsk);
			break;

		case tw::tool_type::conductor:
			ans = std::make_shared<Conductor>(name,ms,tsk);
			break;
		case tw::tool_type::planeWave:
			ans = std::make_shared<PlaneWave>(name,ms,tsk);
			break;
		case tw::tool_type::hermiteGauss:
			ans = std::make_shared<HermiteGauss>(name,ms,tsk);
			break;
		case tw::tool_type::laguerreGauss:
			ans = std::make_shared<LaguerreGauss>(name,ms,tsk);
			break;
		case tw::tool_type::besselBeam:
			ans = std::make_shared<BesselBeam>(name,ms,tsk);
			break;
		case tw::tool_type::airyDisc:
			ans = std::make_shared<AiryDisc>(name,ms,tsk);
			break;
		case tw::tool_type::multipole:
			ans = std::make_shared<Multipole>(name,ms,tsk);
			break;
		case tw::tool_type::uniformProfile:
			ans = std::make_shared<UniformProfile>(name,ms,tsk);
			break;
		case tw::tool_type::piecewiseProfile:
			ans = std::make_shared<PiecewiseProfile>(name,ms,tsk);
			break;
		case tw::tool_type::channelProfile:
			ans = std::make_shared<ChannelProfile>(name,ms,tsk);
			break;
		case tw::tool_type::columnProfile:
			ans = std::make_shared<ColumnProfile>(name,ms,tsk);
			break;
		case tw::tool_type::gaussianProfile:
			ans = std::make_shared<GaussianProfile>(name,ms,tsk);
			break;
		case tw::tool_type::corrugatedProfile:
			ans = std::make_shared<CorrugatedProfile>(name,ms,tsk);
			break;
		case tw::tool_type::eigenmodePropagator:
			ans = std::make_shared<EigenmodePropagator>(name,ms,tsk);
			break;
		case tw::tool_type::adiPropagator:
			ans = std::make_shared<ADIPropagator>(name,ms,tsk);
			break;
		case tw::tool_type::isotropicPropagator:
			ans = std::make_shared<IsotropicPropagator>(name,ms,tsk);
			break;
		case tw::tool_type::generalParabolicPropagator:
			ans = std::make_shared<ParabolicSolver>(name,ms,tsk);
			break;
		case tw::tool_type::schroedingerPropagator:
			ans = std::make_shared<SchroedingerPropagator>(name,ms,tsk);
			break;
		case tw::tool_type::iterativePoissonSolver:
			ans = std::make_shared<IterativePoissonSolver>(name,ms,tsk);
			break;
		case tw::tool_type::ellipticSolver1D:
			ans = std::make_shared<EllipticSolver1D>(name,ms,tsk);
			break;
		case tw::tool_type::facrPoissonSolver:
			ans = std::make_shared<PoissonSolver>(name,ms,tsk);
			break;
		case tw::tool_type::eigenmodePoissonSolver:
			ans = std::make_shared<EigenmodePoissonSolver>(name,ms,tsk);
			break;
		case tw::tool_type::yeePropagatorPML:
			ans = std::make_shared<YeePropagatorPML>(name,ms,tsk);
			break;
		case tw::tool_type::lorentzPropagator:
			ans = std::make_shared<LorentzPropagator>(name,ms,tsk);
			break;
		case tw::tool_type::eosData:
			ans = std::make_shared<EOSComponent>(name,ms,tsk);
			break;
		case tw::tool_type::eosIdealGas:
			ans = std::make_shared<EOSIdealGas>(name,ms,tsk);
			break;
		case tw::tool_type::eosHotElectrons:
			ans = std::make_shared<EOSHotElectrons>(name,ms,tsk);
			break;
		case tw::tool_type::eosMixture:
			ans = std::make_shared<EOSMixture>(name,ms,tsk);
			break;
		case tw::tool_type::eosIdealGasMix:
			ans = std::make_shared<EOSIdealGasMix>(name,ms,tsk);
			break;
		case tw::tool_type::eosSimpleMieGruneisen:
			ans = std::make_shared<EOSSimpleMieGruneisen>(name,ms,tsk);
			break;
		case tw::tool_type::eosLinearMieGruneisen:
			ans = std::make_shared<EOSLinearMieGruneisen>(name,ms,tsk);
			break;
		case tw::tool_type::eosTillotson:
			ans = std::make_shared<EOSTillotson>(name,ms,tsk);
			break;
		case tw::tool_type::mpi:
			ans = std::make_shared<Multiphoton>(name,ms,tsk);
			break;
		case tw::tool_type::adk:
			ans = std::make_shared<ADK>(name,ms,tsk);
			break;
		case tw::tool_type::ppt:
			ans = std::make_shared<PPT>(name,ms,tsk);
			break;
		case tw::tool_type::ppt_tunneling:
			ans = std::make_shared<PPT_Tunneling>(name,ms,tsk);
			break;
		case tw::tool_type::kyh:
			ans = std::make_shared<KYH>(name,ms,tsk);
			break;
		case tw::tool_type::pmpb:
			ans = std::make_shared<PMPB>(name,ms,tsk);
			break;
		case tw::tool_type::boxDiagnostic:
			ans = std::make_shared<BoxDiagnostic>(name,ms,tsk);
			break;
		case tw::tool_type::pointDiagnostic:
			ans = std::make_shared<PointDiagnostic>(name,ms,tsk);
			break;
		case tw::tool_type::volumeDiagnostic:
			ans = std::make_shared<VolumeDiagnostic>(name,ms,tsk);
			break;
		case tw::tool_type::particleOrbits:
			ans = std::make_shared<ParticleOrbits>(name,ms,tsk);
			break;
		case tw::tool_type::phaseSpaceDiagnostic:
			ans = std::make_shared<PhaseSpaceDiagnostic>(name,ms,tsk);
			break;
		case tw::tool_type::boundState:
			ans = std::make_shared<BoundState>(name,ms,tsk);
			break;
		case tw::tool_type::freeState:
			ans = std::make_shared<FreeState>(name,ms,tsk);
			break;
		case tw::tool_type::randomState:
			ans = std::make_shared<RandomState>(name,ms,tsk);
			break;
		case tw::tool_type::tabulatedState:
			ans = std::make_shared<TabulatedState>(name,ms,tsk);
			break;
		case tw::tool_type::photonGenerator:
			ans = std::make_shared<PhotonGenerator>(name,ms,tsk);
			break;
		case tw::tool_type::pairCreator:
			ans = std::make_shared<PairCreator>(name,ms,tsk);
			break;
		case tw::tool_type::borisMover:
			ans = std::make_shared<BorisMover>(name,ms,tsk);
			break;
		case tw::tool_type::hcMover:
			ans = std::make_shared<HCMover>(name,ms,tsk);
			break;
		case tw::tool_type::pgcMover:
			ans = std::make_shared<PGCMover>(name,ms,tsk);
			break;
		case tw::tool_type::unitaryMover:
			ans = std::make_shared<UnitaryMover>(name,ms,tsk);
			break;
		case tw::tool_type::bohmianMover:
			ans = std::make_shared<BohmianMover>(name,ms,tsk);
			break;
		case tw::tool_type::photonMover:
			ans = std::make_shared<PhotonMover>(name,ms,tsk);
			break;
		case tw::tool_type::iteratorTest:
			ans = std::make_shared<IteratorTest>(name,ms,tsk);
			break;
		case tw::tool_type::metricSpaceTest:
			ans = std::make_shared<MetricSpaceTest>(name,ms,tsk);
			break;
		case tw::tool_type::fftTest:
			ans = std::make_shared<FFTTest>(name,ms,tsk);
			break;
		default:
			throw tw::FactoryError("unknown tool type");
	}
	return ans;
}

Module* CreateDriverFromType(const std::string& name,tw::tool_type theType,MetricSpace *ms,Task *tsk)
{
	Module *ans;
	switch (theType)
	{
		case tw::tool_type::curvilinearDirectSolver:
			ans = new CurvilinearDirectSolver(name,ms,tsk);
			break;
		case tw::tool_type::electrostatic:
			ans = new Electrostatic(name,ms,tsk);
			break;
		case tw::tool_type::coulombSolver:
			ans = new CoulombSolver(name,ms,tsk);
			break;
		case tw::tool_type::directSolver:
			ans = new DirectSolver(name,ms,tsk);
			break;
		case tw::tool_type::farFieldDiagnostic:
			ans = new FarFieldDiagnostic(name,ms,tsk);
			break;
		case tw::tool_type::qsLaser:
			ans = new QSSolver(name,ms,tsk);
			break;
		case tw::tool_type::pgcLaser:
			ans = new PGCSolver(name,ms,tsk);
			break;
		case tw::tool_type::boundElectrons:
			ans = new BoundElectrons(name,ms,tsk);
			break;
		case tw::tool_type::fluidFields:
			ans = new Fluid(name,ms,tsk);
			break;
		case tw::tool_type::equilibriumGroup:
			ans = new EquilibriumGroup(name,ms,tsk);
			break;
		case tw::tool_type::chemical:
			ans = new Chemical(name,ms,tsk);
			break;
		case tw::tool_type::sparcHydroManager:
			ans = new sparc::HydroManager(name,ms,tsk);
			break;
		case tw::tool_type::species:
			ans = new Species(name,ms,tsk);
			break;
		case tw::tool_type::kinetics:
			ans = new Kinetics(name,ms,tsk);
			break;
		case tw::tool_type::schroedinger:
			ans = new Schroedinger(name,ms,tsk);
			break;
		case tw::tool_type::pauli:
			ans = new Pauli(name,ms,tsk);
			break;
		case tw::tool_type::kleinGordon:
			ans = new KleinGordon(name,ms,tsk);
			break;
		case tw::tool_type::dirac:
			ans = new Dirac(name,ms,tsk);
			break;
		case tw::tool_type::populationDiagnostic:
			ans = new PopulationDiagnostic(name,ms,tsk);
			break;
		default:
			throw tw::FactoryError("unknown driver type");
	}
	return ans;
}

}