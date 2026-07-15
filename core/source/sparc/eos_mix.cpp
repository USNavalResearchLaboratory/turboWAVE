module;

#include "tw_includes.h"
#include "tw_logger.h"

export module eos:eos_mix;
import :eos_component;
import input;
import driver;
import fields;
import functions;
import hydro_primitives;
import numerics;
import logger;

/**
 * @brief Provide a mixing rule for EOS
 * 
 */
export struct EOSMixture:ComputeTool
{
	sparc::characteristic_values tiny;
	sparc::hydro_set hidx;
	sparc::eos_set eidx;
	sparc::material_set matset;
	EOSMixture(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk) {}
	void Setup(const sparc::hydro_set& h,const sparc::eos_set& e,const sparc::material_set& m,const sparc::characteristic_values& tiny)
	{
		this->tiny = tiny;
		hidx = h;
		eidx = e;
		matset = m;
	}
	tw::Float DensitySum(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += f(cell,s+hidx.first);
		return ans;
	}
	tw::Float MassDensity(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += f(cell,s+hidx.first)*matset.mass[s];
		return ans;
	}
	tw::Float MixVibrationalEnergy(const Field& f,const tw::cell& cell, const tw::Float Tv)
	{
		tw::Float n_tot = tiny.u;
		tw::Float eps_n_tot = 0.0;
		for (tw::Int s=0;s<hidx.num;s++) {
			const tw::Float ns = f(cell,s+hidx.first);
			n_tot += ns * (matset.excitationEnergy[s]>0.0);
			eps_n_tot += ns * matset.excitationEnergy[s];
		}
		return eps_n_tot/std::exp((eps_n_tot + tiny.u)/(Tv*n_tot + tiny.u) - 1.0);
	}
	tw::Float MixVibrationalTemperature(const Field& f,const tw::cell& cell)
	{
		tw::Float xi = tiny.u + f(cell,hidx.x);
		tw::Float n_tot = tiny.u;
		tw::Float eps_n_tot = 0.0;
		for (tw::Int s=0;s<hidx.num;s++) {
			const tw::Float ns = f(cell,s+hidx.first);
			n_tot += ns * (matset.excitationEnergy[s]>0.0);
			eps_n_tot += ns*matset.excitationEnergy[s];
		}
		return eps_n_tot / n_tot / std::log(1 + tw::eps_pos + eps_n_tot/xi);
	}
	tw::Float InternalEnergy(const tw::Float& nm,const Field& f,const tw::cell& cell)
	{
		const tw::vec3 np = tw::vec3(f(cell,hidx.npx),f(cell,hidx.npy),f(cell,hidx.npz));
		const tw::vec3 vel = np/(tiny.n + nm);
		const tw::Float KE = 0.5*nm*Norm(vel);
		const tw::Float primitive = f(cell,hidx.u) - KE;
		const tw::Float failsafe = 1e-6*KE + tiny.u;
		const tw::Float sel = tw::Float(primitive>0.0);
		return sel*primitive + (1.0-sel)*failsafe;
	}
	/**
	 * @brief load the mass density and internal energy density associated with a given component
	 * 
	 * @param[in] i index of component's density in hydro field
	 * @param[out] nm mass density result
	 * @param[out] nmE internal energy density result
	 * @param[in] hydro field to evaluate
	 */
	void LoadPartialExtrinsics(tw::Int i, ScalarField& nm, ScalarField& nmE, Field& hydro) {
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				tw::Float nm_tot = MassDensity(hydro,cell);
				nm(cell) = hydro(cell,i) * matset.mass[i - hidx.first];
				nmE(cell) = InternalEnergy(nm_tot,hydro,cell) * nm(cell) / (tiny.n + nm_tot);
			}
		}
	}
	/**
	 * @brief Update the internal energy after there has been a temperature change due to heat transport.
	 *        This is only first order accurate, but usually small.
	 * 
	 * @param[in] T0 temperature before the heat transport
	 * @param[out] hydro field to update
	 * @param[in] eos field with the new temperature
	 */
	virtual void FinishHeatTransport(ScalarField& T0,Field& hydro,Field& eos)
	{
		// Not centered, because cv is not updated.
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1)) {
				hydro(cell,hidx.u) += eos(cell,eidx.nmcv) * (eos(cell,eidx.T) - T0(cell));
			}
		}
	}
	/**
	 * @brief Add energy needed to achieve the given temperature or pressure.
	 * @details Total energy should be kinetic part upon entry.  The nmcv element of the EOS should not be used as
	 * it will not in general be initialized.
	 * 
	 * @param[in] elements the EOS for each constituent of the mixture
	 * @param[in,out] hydro field providing the density, will also receive the energy
	 * @param[in] eos field providing the target intrinsic, the non-zero one will be used
	 */
	virtual void InitEnergyWithIntrinsics(std::vector<EOSComponent*> elements, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			tw::Float guess1 = 0;
			tw::Float guess2 = 0;
			// unitialized far ghost cells can cause trouble for this loop
			for (auto cell : InteriorCellRange(*space,1)) {
				auto nm = tiny.n + MassDensity(hydro,cell);
				auto target_P = eos(cell,eidx.P);
				auto target_T = eos(cell,eidx.T);
				if (target_T > 0) {
					for (auto c : elements) {
						auto nm = c->mat.mass * hydro(cell,c->hidx.ni);
						hydro(cell,hidx.u) += c->InternalEnergy(nm,target_T);
					}
					hydro(cell,hidx.x) += MixVibrationalEnergy(hydro, cell, target_T);
				} else if (target_P > 0) {
					auto merit = [this,&cell,&elements,hydro,nm,target_P] (tw::Float nmE) {
						tw::Float press = 0;
						for (auto c : elements) {
							auto partial_nm = c->mat.mass * hydro(cell,c->hidx.ni);
							auto partial_nmE = nmE * partial_nm / (tiny.n + nm);
							press += c->Pressure(partial_nm,partial_nmE);
						}
						return press - target_P;
					};
					if (guess1 == guess2) {
						guess1 = target_P*1.4;
						guess2 = target_P*1.5;
					}
					guess1 = SecantMethod(merit,guess1,guess2,target_P*1e-8+tiny.u,guess1*1e-5+tiny.u,100);
					hydro(cell,hidx.u) += guess1;
					guess2 = guess1*1.01;
				}
			}
		}
		hydro.CopyFromNeighbors(Rng(hidx.u));
		hydro.ApplyBoundaryCondition(Rng(hidx.u));
	}
	/**
	 * @brief Initialize temperature without using any reference state.
	 * @details This is not always a valid procedure, i.e., in some cases a reference state
	 * may be required. This may be redundant in principle, e.g. in cases where a temperature
	 * was specified and the energy chosen to be consistent.
	 * 
	 * @param[in] elements constituents of this mixture
	 * @param[in] hydro new density, momentum, energy
	 * @param[in,out] eos T, Tv, and nmcv are initialized
	 */
	virtual void InitTemperature(std::vector<EOSComponent*> elements, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				const tw::Float nm1 = MassDensity(hydro,cell);
				const tw::Float nmE1 = InternalEnergy(nm1,hydro,cell);

				tw::Float nmcv = 0.0;
				for (auto c : elements) {
					tw::Float frac = c->mat.mass * hydro(cell,c->hidx.ni) / (tiny.n + nm1);
					nmcv += c->HeatCapacity(frac*nm1,frac*nmE1);
				}

				eos(cell,eidx.nmcv) = nmcv;
				eos(cell,eidx.T) = nmE1 / (tiny.n + nmcv);
				eos(cell,eidx.Tv) = MixVibrationalTemperature(hydro, cell);
			}
		}
	}
	/**
	 * @brief Update temperature using reference data (e.g. previous time level).
	 * @details The default presumes polytropic ideal gases and ignores reference states.
	 * 
	 * @param[in] elements constituents of this mixture
	 * @param[in] hydroRef reference density, momentum, energy
	 * @param[in] hydro new density, momentum, energy
	 * @param[in] eosRef reference temperature, pressure, heat capacity
	 * @param[in,out] eos T, Tv, and nmcv are updated
	 */
	virtual void UpdateTemperature(std::vector<EOSComponent*> elements, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				const tw::Float nm1 = MassDensity(hydro,cell);
				const tw::Float nmE1 = InternalEnergy(nm1,hydro,cell);

				tw::Float nmcv = 0.0;
				for (auto c : elements) {
					tw::Float frac = c->mat.mass * hydro(cell,c->hidx.ni) / (tiny.n + nm1);
					nmcv += c->HeatCapacity(frac*nm1, frac*nmE1);
				}

				eos(cell,eidx.nmcv) = nmcv;
				eos(cell,eidx.T) = nmE1 / (tiny.n + nmcv);
				eos(cell,eidx.Tv) = MixVibrationalTemperature(hydro, cell);
			}
		}
	}
};

export struct EOSIdealGasMix:EOSMixture
{
	EOSIdealGasMix(const std::string& name,MetricSpace *m,Task *tsk) : EOSMixture(name,m,tsk) {}
	/**
	 * @brief Update temperature assuming polytropic ideal gas mix.
	 * 
	 * @param[in] elements constituents of this mixture
	 * @param[in] hydroRef ignored
	 * @param[in] hydro new density, momentum, energy
	 * @param[in] eosRef ignored
	 * @param[in,out] eos T, Tv, and nmcv are updated
	 */
	virtual void UpdateTemperature(std::vector<EOSComponent*> elements, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos)
	{
		EOSMixture::UpdateTemperature(elements, hydroRef, hydro, eosRef, eos);
	}
};

export struct EOSGenericMix:EOSMixture
{
	EOSGenericMix(const std::string& name,MetricSpace *m,Task *tsk) : EOSMixture(name,m,tsk) {}
	/**
	 * @brief Update temperature accounting for cold curve
	 * @details This version does not use reference data, but instead uses a formula similar
	 * to the polytropic ideal gas, except that the cold energy is subtracted out.
	 * 
	 * @param[in] elements constituents of this mixture
	 * @param[in] hydroRef reference density, momentum, energy
	 * @param[in] hydro new density, momentum, energy
	 * @param[in] eosRef reference temperature, pressure, heat capacity
	 * @param[in,out] eos T, Tv, and nmcv are updated
	 */
	virtual void UpdateTemperature(std::vector<EOSComponent*> elements, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				const tw::Float nm1 = MassDensity(hydro,cell);
				const tw::Float nmE1 = InternalEnergy(nm1,hydro,cell);

				tw::Float nmcv = 0.0;
				tw::Float cold = 0.0;
				for (auto c : elements) {
					tw::Float frac = c->mat.mass * hydro(cell,c->hidx.ni) / (tiny.n + nm1);
					nmcv += c->HeatCapacity(frac*nm1, frac*nmE1);
					cold += c->ColdCurveGet(nm1);
				}
				if (cold > nmE1) {
					cold = nmE1 - tiny.u;
				}

				eos(cell,eidx.nmcv) = nmcv;
				eos(cell,eidx.T) = (nmE1 - cold) / (tiny.n + nmcv);
				eos(cell,eidx.Tv) = MixVibrationalTemperature(hydro, cell);
			}
		}
	}
};
