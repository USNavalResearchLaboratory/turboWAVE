module;

#include "tw_includes.h"

export module eos:eos_mix;
import :eos_component;
import input;
import driver;
import fields;
import functions;
import hydro_primitives;
import numerics;

/// Tool for handling non-additive quantities in a mixture.
/// An important task is calculating temperature given various inputs.
/// Additive quantities like pressure tend to be handled by components.
export struct EOSMixture:ComputeTool
{
	sparc::hydro_set hidx;
	sparc::eos_set eidx;
	sparc::material_set matset;
	EOSMixture(const std::string& name,MetricSpace *m,Task *tsk) : ComputeTool(name,m,tsk) {}
	void SetupIndexing(const sparc::hydro_set& h,const sparc::eos_set& e,const sparc::material_set& m)
	{
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
	tw::Float MixVibrationalEnergy(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += f(cell,s+hidx.first)*matset.excitationEnergy[s];
		return ans;
	}
	tw::Float MixVibrationalStates(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=0;s<hidx.num;s++)
			ans += matset.excitationEnergy[s] > 0.0 ? f(cell,s+hidx.first) : 0.0;
		return ans;
	}
	tw::Float InternalEnergy(const tw::Float& nm,const Field& f,const tw::cell& cell)
	{
		const tw::vec3 np = tw::vec3(f(cell,hidx.npx),f(cell,hidx.npy),f(cell,hidx.npz));
		const tw::vec3 vel = np/(tw::small_pos + nm);
		const tw::Float KE = 0.5*nm*Norm(vel);
		const tw::Float primitive = f(cell,hidx.u) - KE;
		const tw::Float failsafe = 1e-6*KE + tw::small_pos;
		const tw::Float sel = tw::Float(primitive>0.0);
		return sel*primitive + (1.0-sel)*failsafe;
	}
	/**
	 * @brief Add energy needed to achieve the given temperature or pressure, internal energy should be zero on entry.
	 * 
	 * @param[in] elements the EOS for each constituent of the mixture
	 * @param[in,out] hydro field providing the density, will also receive the energy
	 * @param[in] eos field providing the target pressure
	 */
	virtual void InitEnergyWithIntrinsics(std::vector<EOSComponent*> elements, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1)) {
				auto nm = tw::small_pos + MassDensity(hydro,cell);
				auto target_press = eos(cell,eidx.P);
				auto target_T = eos(cell,eidx.T);
				if (target_T > 0) {
					for (auto c : elements) {
						auto n = hydro(cell,c->hidx.ni);
						hydro(cell,hidx.u) += c->InternalEnergy(n,target_T);
					}
				} else if (target_press > 0) {
					auto merit = [&cell,&elements,hydro,nm,target_press] (tw::Float IE) {
						tw::Float press = 0;
						for (auto c : elements) {
							auto n = hydro(cell,c->hidx.ni);
							auto partial_IE = IE * c->mat.mass * n / nm;
							press += c->Pressure(partial_IE,n);
						}
						return press - target_press;
					};
					hydro(cell,hidx.u) += SecantMethod(merit,1.4*target_press,1.5*target_press);
				}
			}
		}
		hydro.ApplyBoundaryCondition(Rng(hidx.u));
	}
	/**
	 * @brief Initialize temperature without using reference data (not always valid)
	 * 
	 * @param[out] IE new internal energy
	 * @param[out] nm new mass density 
	 * @param[in] hydro new density, momentum, energy
	 * @param[in,out] eos new values, temperature is updated
	 */
	virtual void InitTemperature(ScalarField& IE, ScalarField& nm, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				nm(cell) = MassDensity(hydro,cell);
				IE(cell) = InternalEnergy(nm(cell),hydro,cell);
				const tw::Float epsvn = MixVibrationalEnergy(hydro,cell);
				const tw::Float nv = MixVibrationalStates(hydro,cell);

				eos(cell,eidx.T) = IE(cell)/(tw::small_pos + eos(cell,eidx.nmcv));
				eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/std::log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
			}
		}
	}
	/**
	 * @brief Update temperature using reference data (e.g. previous time level)
	 * 
	 * @param[out] IE new internal energy
	 * @param[out] nm new mass density 
	 * @param[in] hydroRef reference density, momentum, energy
	 * @param[in] hydro new density, momentum, energy
	 * @param[in] eosRef reference temperature, pressure, heat capacity
	 * @param[in,out] eos new values, temperature is updated
	 */
	virtual void UpdateTemperature(ScalarField& IE, ScalarField& nm, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				const tw::Float nm1 = tw::small_pos + MassDensity(hydro,cell);
				const tw::Float IE1 = InternalEnergy(nm1,hydro,cell);
				const tw::Float nm0 = tw::small_pos + MassDensity(hydroRef,cell);
				const tw::Float IE0 = InternalEnergy(nm0,hydroRef,cell);
				const tw::Float epsvn = MixVibrationalEnergy(hydro,cell);
				const tw::Float nv = MixVibrationalStates(hydro,cell);
				const tw::Float nmcv_sum = tw::small_pos + eosRef(cell,eidx.nmcv) + eos(cell,eidx.nmcv);

				//eos(cell,eidx.T) = eosRef(cell,eidx.T) + 2.0*(IE1 - IE0)/nmcv_sum; // wrong
				eos(cell,eidx.T) = eosRef(cell,eidx.T) + (IE1*(1+nm0/nm1) - IE0*(1+nm1/nm0))/nmcv_sum;
				eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/std::log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
				nm(cell) = nm1;
				IE(cell) = IE1;
			}
		}
	}
	/// Add energy corresponding to a change in temperature only.
	virtual void UpdateEnergy(ScalarField& nm,ScalarField& T0,Field& hydro,Field& eos)
	{
		// Not centered, because cv is not updated.
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1)) {
				hydro(cell,hidx.u) += eos(cell,eidx.nmcv) * (eos(cell,eidx.T) - T0(cell));
			}
		}
	}
};

/// Mixture that uses a polytropic ideal gas caloric EOS, reference states are ignored
export struct EOSIdealGasMix:EOSMixture
{
	EOSIdealGasMix(const std::string& name,MetricSpace *m,Task *tsk) : EOSMixture(name,m,tsk) {}
	virtual void UpdateTemperature(ScalarField& IE, ScalarField& nm, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos) {
		InitTemperature(IE,nm,hydro,eos);
	}
};
