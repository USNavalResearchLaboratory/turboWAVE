module;

#include "tw_includes.h"

export module factory;
import base;
import compute_tool;
import twmodule;
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
	std::map<std::string,tw::module_type> module_map = Module::Map();
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

	bool is_module = module_map.find(normalized) != module_map.end();
	bool is_tool = tool_map.find(normalized) != tool_map.end();
	bool is_region = normalized == "region";
	if (is_module && is_tool || is_module && is_region || is_tool && is_region) {
		throw tw::FatalError(std::format("overlapping key {}, this is a bug",normalized));
	}
	if (is_module || is_tool || is_region) {
		return true;
	} else {
		for (const auto& pair : module_map) {
			tw::input::BuildSimilar(similar_keys,normalized,pair.first);
		}
		for (const auto& pair : tool_map) {
			tw::input::BuildSimilar(similar_keys,normalized,pair.first);
		}
		tw::input::BuildSimilar(similar_keys,normalized,"region");
		return false;
	}
}

ComputeTool* CreateToolFromType(const std::string& name,tw::tool_type theType,MetricSpace *ms,Task *tsk)
{
	ComputeTool *ans;
	switch (theType)
	{
		case tw::tool_type::none:
			ans = NULL;
			break;
		case tw::tool_type::warp:
			ans = new Warp(name,ms,tsk);
			break;
		case tw::tool_type::conductor:
			ans = new Conductor(name,ms,tsk);
			break;
		case tw::tool_type::planeWave:
			ans = new PlaneWave(name,ms,tsk);
			break;
		case tw::tool_type::hermiteGauss:
			ans = new HermiteGauss(name,ms,tsk);
			break;
		case tw::tool_type::laguerreGauss:
			ans = new LaguerreGauss(name,ms,tsk);
			break;
		case tw::tool_type::besselBeam:
			ans = new BesselBeam(name,ms,tsk);
			break;
		case tw::tool_type::airyDisc:
			ans = new AiryDisc(name,ms,tsk);
			break;
		case tw::tool_type::multipole:
			ans = new Multipole(name,ms,tsk);
			break;
		case tw::tool_type::uniformProfile:
			ans = new UniformProfile(name,ms,tsk);
			break;
		case tw::tool_type::piecewiseProfile:
			ans = new PiecewiseProfile(name,ms,tsk);
			break;
		case tw::tool_type::channelProfile:
			ans = new ChannelProfile(name,ms,tsk);
			break;
		case tw::tool_type::columnProfile:
			ans = new ColumnProfile(name,ms,tsk);
			break;
		case tw::tool_type::gaussianProfile:
			ans = new GaussianProfile(name,ms,tsk);
			break;
		case tw::tool_type::corrugatedProfile:
			ans = new CorrugatedProfile(name,ms,tsk);
			break;
		case tw::tool_type::eigenmodePropagator:
			ans = new EigenmodePropagator(name,ms,tsk);
			break;
		case tw::tool_type::adiPropagator:
			ans = new ADIPropagator(name,ms,tsk);
			break;
		case tw::tool_type::isotropicPropagator:
			ans = new IsotropicPropagator(name,ms,tsk);
			break;
		case tw::tool_type::generalParabolicPropagator:
			ans = new ParabolicSolver(name,ms,tsk);
			break;
		case tw::tool_type::schroedingerPropagator:
			ans = new SchroedingerPropagator(name,ms,tsk);
			break;
		case tw::tool_type::iterativePoissonSolver:
			ans = new IterativePoissonSolver(name,ms,tsk);
			break;
		case tw::tool_type::ellipticSolver1D:
			ans = new EllipticSolver1D(name,ms,tsk);
			break;
		case tw::tool_type::facrPoissonSolver:
			ans = new PoissonSolver(name,ms,tsk);
			break;
		case tw::tool_type::eigenmodePoissonSolver:
			ans = new EigenmodePoissonSolver(name,ms,tsk);
			break;
		case tw::tool_type::yeePropagatorPML:
			ans = new YeePropagatorPML(name,ms,tsk);
			break;
		case tw::tool_type::lorentzPropagator:
			ans = new LorentzPropagator(name,ms,tsk);
			break;
		case tw::tool_type::eosData:
			ans = new EOSComponent(name,ms,tsk);
			break;
		case tw::tool_type::eosIdealGas:
			ans = new EOSIdealGas(name,ms,tsk);
			break;
		case tw::tool_type::eosHotElectrons:
			ans = new EOSHotElectrons(name,ms,tsk);
			break;
		case tw::tool_type::eosMixture:
			ans = new EOSMixture(name,ms,tsk);
			break;
		case tw::tool_type::eosIdealGasMix:
			ans = new EOSIdealGasMix(name,ms,tsk);
			break;
		case tw::tool_type::eosSimpleMieGruneisen:
			ans = new EOSSimpleMieGruneisen(name,ms,tsk);
			break;
		case tw::tool_type::eosLinearMieGruneisen:
			ans = new EOSLinearMieGruneisen(name,ms,tsk);
			break;
		case tw::tool_type::eosTillotson:
			ans = new EOSTillotson(name,ms,tsk);
			break;
		case tw::tool_type::mpi:
			ans = new Multiphoton(name,ms,tsk);
			break;
		case tw::tool_type::adk:
			ans = new ADK(name,ms,tsk);
			break;
		case tw::tool_type::ppt:
			ans = new PPT(name,ms,tsk);
			break;
		case tw::tool_type::ppt_tunneling:
			ans = new PPT_Tunneling(name,ms,tsk);
			break;
		case tw::tool_type::kyh:
			ans = new KYH(name,ms,tsk);
			break;
		case tw::tool_type::pmpb:
			ans = new PMPB(name,ms,tsk);
			break;
		case tw::tool_type::boxDiagnostic:
			ans = new BoxDiagnostic(name,ms,tsk);
			break;
		case tw::tool_type::pointDiagnostic:
			ans = new PointDiagnostic(name,ms,tsk);
			break;
		case tw::tool_type::volumeDiagnostic:
			ans = new VolumeDiagnostic(name,ms,tsk);
			break;
		case tw::tool_type::particleOrbits:
			ans = new ParticleOrbits(name,ms,tsk);
			break;
		case tw::tool_type::phaseSpaceDiagnostic:
			ans = new PhaseSpaceDiagnostic(name,ms,tsk);
			break;
		case tw::tool_type::boundState:
			ans = new BoundState(name,ms,tsk);
			break;
		case tw::tool_type::freeState:
			ans = new FreeState(name,ms,tsk);
			break;
		case tw::tool_type::randomState:
			ans = new RandomState(name,ms,tsk);
			break;
		case tw::tool_type::tabulatedState:
			ans = new TabulatedState(name,ms,tsk);
			break;
		case tw::tool_type::photonGenerator:
			ans = new PhotonGenerator(name,ms,tsk);
			break;
		case tw::tool_type::pairCreator:
			ans = new PairCreator(name,ms,tsk);
			break;
		case tw::tool_type::borisMover:
			ans = new BorisMover(name,ms,tsk);
			break;
		case tw::tool_type::hcMover:
			ans = new HCMover(name,ms,tsk);
			break;
		case tw::tool_type::pgcMover:
			ans = new PGCMover(name,ms,tsk);
			break;
		case tw::tool_type::unitaryMover:
			ans = new UnitaryMover(name,ms,tsk);
			break;
		case tw::tool_type::bohmianMover:
			ans = new BohmianMover(name,ms,tsk);
			break;
		case tw::tool_type::photonMover:
			ans = new PhotonMover(name,ms,tsk);
			break;
		case tw::tool_type::iteratorTest:
			ans = new IteratorTest(name,ms,tsk);
			break;
		case tw::tool_type::metricSpaceTest:
			ans = new MetricSpaceTest(name,ms,tsk);
			break;
		case tw::tool_type::fftTest:
			ans = new FFTTest(name,ms,tsk);
			break;
	}
	return ans;
}

Module* CreateModuleFromType(const std::string& name,tw::module_type theType,Simulation *sim)
{
	Module *ans;
	switch (theType)
	{
		case tw::module_type::none:
			ans = NULL;
			break;
		case tw::module_type::curvilinearDirectSolver:
			ans = new CurvilinearDirectSolver(name,sim);
			break;
		case tw::module_type::electrostatic:
			ans = new Electrostatic(name,sim);
			break;
		case tw::module_type::coulombSolver:
			ans = new CoulombSolver(name,sim);
			break;
		case tw::module_type::directSolver:
			ans = new DirectSolver(name,sim);
			break;
		case tw::module_type::farFieldDiagnostic:
			ans = new FarFieldDiagnostic(name,sim);
			break;
		case tw::module_type::qsLaser:
			ans = new QSSolver(name,sim);
			break;
		case tw::module_type::pgcLaser:
			ans = new PGCSolver(name,sim);
			break;
		case tw::module_type::boundElectrons:
			ans = new BoundElectrons(name,sim);
			break;
		case tw::module_type::fluidFields:
			ans = new Fluid(name,sim);
			break;
		case tw::module_type::equilibriumGroup:
			ans = new EquilibriumGroup(name,sim);
			break;
		case tw::module_type::chemical:
			ans = new Chemical(name,sim);
			break;
		case tw::module_type::sparcHydroManager:
			ans = new sparc::HydroManager(name,sim);
			break;
		case tw::module_type::species:
			ans = new Species(name,sim);
			break;
		case tw::module_type::kinetics:
			ans = new Kinetics(name,sim);
			break;
		case tw::module_type::schroedinger:
			ans = new Schroedinger(name,sim);
			break;
		case tw::module_type::pauli:
			ans = new Pauli(name,sim);
			break;
		case tw::module_type::kleinGordon:
			ans = new KleinGordon(name,sim);
			break;
		case tw::module_type::dirac:
			ans = new Dirac(name,sim);
			break;
		case tw::module_type::populationDiagnostic:
			ans = new PopulationDiagnostic(name,sim);
			break;
	}
	return ans;
}

}