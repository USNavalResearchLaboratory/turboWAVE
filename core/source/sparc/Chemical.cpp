module;

#include "tw_includes.h"
#include "tw_logger.h"

export module hydro:chemical;
import input;
import driver;
import fields;
import hydro_primitives;
import eos;
import photoionization;
import injection;
import logger;

export struct Chemical:Driver
{
	tw::Profiles profiles;
	std::shared_ptr<Ionizer> ionizer;
	std::shared_ptr<EOSComponent> eosData;
	std::shared_ptr<UniformProfile> background;

	sparc::material mat;
	tw::Int indexInState;

	Chemical(const std::string& name,MetricSpace *ms, Task *tsk) : Driver(name,ms,tsk)
	{
		if (native.unit_system!=tw::units::plasma)
			throw tw::FatalError("Chemical module requires <native units = plasma>");

		mat.charge = -1.0;
		mat.mass = 1.0;
		mat.cvm = 1.5;
		mat.excitationEnergy = 0.0;
		mat.thermometricConductivity = 0.0;
		mat.kinematicViscosity = 0.0;
		mat.eps[0] = 1.0;
		mat.eps[1] = 0.0;
		mat.AddDirectives(directives);
	}

	void VerifyInput()
	{
		Driver::VerifyInput();
		// search tools for EOS, photoionization, and profiles
		for (auto tool : tools) {
			if (std::dynamic_pointer_cast<EOSComponent>(tool)) {
				eosData = std::dynamic_pointer_cast<EOSComponent>(tool);
			} else if (std::dynamic_pointer_cast<Ionizer>(tool)) {
				ionizer = std::dynamic_pointer_cast<Ionizer>(tool);
			} else if (std::dynamic_pointer_cast<Profile>(tool)) {
				profiles.push_back(std::dynamic_pointer_cast<Profile>(tool));
			}
		}
		// If the EOS tool could not be found, create one automatically.
		if (!eosData) {
			auto new_tool = mat.mass==1.0 ?
				CreateTool("default_hot_electrons",tw::tool_type::eosHotElectrons) :
				CreateTool("default_ideal_gas",tw::tool_type::eosIdealGas);
			AddTool(new_tool);
			eosData = std::dynamic_pointer_cast<EOSComponent>(new_tool);
		}
		// Add a uniform profile for the automatic background fluid.
		auto new_tool = CreateTool("auto_background",tw::tool_type::uniformProfile);
		AddTool(new_tool);
		background = std::dynamic_pointer_cast<UniformProfile>(new_tool);
		profiles.push_back(background);
	}

	/// Set EOS indexing, and return this chemical's photoionization object, which the
	/// caller is expected to index, if it is not empty.
	/// Assumes mat and indexInState are valid.
	std::shared_ptr<Ionizer> SetupIndexing(const sparc::hydro_set& hidx,const sparc::eos_set& eidx)
	{
		logger::TRACE(std::format("indexing {} EOS",name));
		eosData->SetupIndexing(indexInState,hidx,eidx,mat);
		return ionizer;
	}

	/// Handle external sources or sinks that have been prescribed as part of the problem
	void PumpFluid(Field& create,Field& destroy,const sparc::hydro_set& hidx) {
		const tw::Int ns = indexInState;
		const tw::Int npx = hidx.npx;
		const tw::Int npy = hidx.npy;
		const tw::Int npz = hidx.npz;
		const tw::Int U = hidx.u;
		const tw::Int Xi = hidx.x;
		tw::Float add;

		for (auto prof : profiles) {
			if ( prof->whichQuantity==tw::profile::quantity::power && prof->TimeGate(space->WindowPos(0),&add) ) {
				for (auto cell : EntireCellRange(*this,1)) {
					create(cell,U) += prof->GetValue(space->Pos4(cell),*space);
				}
			}
		}
	}

	/// Initial loading of mass, momentum, and energy.
	/// Internal energy associated with a temperature or pressure specification is *not* handled herein.
	bool LoadFluid(Field& hydro,const sparc::hydro_set& hidx)
	{
		bool massLoaded = false;
		tw::Float add = 0.0;

		const tw::Int ns = indexInState;
		const tw::Int npx = hidx.npx;
		const tw::Int npy = hidx.npy;
		const tw::Int npz = hidx.npz;
		const tw::Int U = hidx.u;
		const tw::Int Xi = hidx.x;

		for (auto prof : profiles)
		{
			logger::TRACE(std::format("loading profile <{}>",prof->name));
			if ( prof->TimeGate(space->WindowPos(0),&add) )
			{
				const tw::vec3 p0 = prof->DriftMomentum(mat.mass);
				for (auto cell : EntireCellRange(*this,1))
				{
					const tw::Float dens = prof->GetValue(space->Pos4(cell),*space);
					if (prof->whichQuantity==tw::profile::quantity::density && dens>0.0)
					{
						massLoaded = true;
						const tw::Float kT = prof->Temperature(mat.mass);
						const tw::Float kinetic = 0.5*Norm(dens*p0)/(tw::small_pos + mat.mass*dens);
						const tw::Float vibrational = dens*mat.excitationEnergy/(std::fabs(std::exp(mat.excitationEnergy/kT) - 1.0) + tw::small_pos);
						hydro(cell,ns) = add*hydro(cell,ns) + dens;
						hydro(cell,npx) = add*hydro(cell,npx) + dens*p0.x;
						hydro(cell,npy) = add*hydro(cell,npy) + dens*p0.y;
						hydro(cell,npz) = add*hydro(cell,npz) + dens*p0.z;
						hydro(cell,U) = add*hydro(cell,U) + kinetic + vibrational; // internal energy added in subsequent sweep
						hydro(cell,Xi) = add*hydro(cell,Xi) + vibrational;
					} else if (prof->whichQuantity==tw::profile::quantity::energy) {
						hydro(cell,U) = add*hydro(cell,U) + dens;
					} else if (prof->whichQuantity==tw::profile::quantity::power) {
						// do nothing at init time
					} else if (prof->whichQuantity==tw::profile::quantity::px) {
						hydro(cell,npx) = add*hydro(cell,npx) + dens;
					} else if (prof->whichQuantity==tw::profile::quantity::py) {
						hydro(cell,npy) = add*hydro(cell,npy) + dens;
					} else if (prof->whichQuantity==tw::profile::quantity::pz) {
						hydro(cell,npz) = add*hydro(cell,npz) + dens;
					}
				}
			}
		}
		hydro.ApplyBoundaryCondition(Rng(ns,Xi+1));
		return massLoaded;
	}

	/**
	 * @brief Set the target for temperature or pressure for the owning group based on attached profiles.
	 *        Only one or the other of T or P should be nonzero, and overlapping profiles should be consistent.
	 * 
	 * @param[out] eos 
	 * @param[in] eidx 
	 */
	void LoadTargetParameters(Field& eos,const sparc::eos_set& eidx)
	{
		for (auto prof : profiles) {
			if (prof->whichQuantity==tw::profile::quantity::density) {
				// Get the target parameters for this profile
				const tw::Float kT = prof->temperature; // want direct access so 0 is detected
				const tw::Float P = prof->pressure; // want direct access so 0 is detected
				for (auto cell : EntireCellRange(*this,1)) {
					if (prof->GetValue(space->Pos4(cell),*space) > 0) {
						if (kT > 0) {
							eos(cell,eidx.T) = kT;
						} else if (P > 0) {
							eos(cell,eidx.P) = P;
						}
					}
				}
			}
		}
	}
};