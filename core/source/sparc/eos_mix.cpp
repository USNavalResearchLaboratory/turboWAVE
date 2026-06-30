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
	 * @brief load the mass density and internal energy density associated with a given component
	 * 
	 * @param[in] i index of component's density in hydro field
	 * @param[out] nm mass density result
	 * @param[out] IE internal energy density result
	 * @param[in] hydro field to evaluate
	 */
	void LoadPartialExtrinsics(tw::Int i, ScalarField& nm, ScalarField& IE, Field& hydro) {
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				tw::Float nm_tot = MassDensity(hydro,cell);
				nm(cell) = hydro(cell,i) * matset.mass[i - hidx.first];
				IE(cell) = InternalEnergy(nm_tot,hydro,cell) * nm(cell) / nm_tot;
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
	 * @brief Add energy needed to achieve the given temperature or pressure, internal energy should be zero on entry.
	 * 
	 * @param[in] elements the EOS for each constituent of the mixture
	 * @param[in,out] hydro field providing the density, will also receive the energy
	 * @param[in] eos field providing the target intrinsic, the non-zero one will be used
	 */
	virtual void InitEnergyWithIntrinsics(std::vector<EOSComponent*> elements, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1)) {
				auto nm = tw::small_pos + MassDensity(hydro,cell);
				auto target_P = eos(cell,eidx.P);
				auto target_T = eos(cell,eidx.T);
				if (target_T > 0) {
					for (auto c : elements) {
						auto nm = c->mat.mass * hydro(cell,c->hidx.ni);
						hydro(cell,hidx.u) += c->InternalEnergy(nm,target_T);
					}
				} else if (target_P > 0) {
					auto merit = [&cell,&elements,hydro,nm,target_P] (tw::Float IE) {
						tw::Float press = 0;
						for (auto c : elements) {
							auto n = hydro(cell,c->hidx.ni);
							auto partial_IE = IE * c->mat.mass * n / nm;
							press += c->Pressure(n*c->mat.mass,partial_IE);
						}
						return press - target_P;
					};
					hydro(cell,hidx.u) += SecantMethod(merit,1.4*target_P,1.5*target_P);
				}
			}
		}
		hydro.ApplyBoundaryCondition(Rng(hidx.u));
	}
	/**
	 * @brief Initialize temperature without using any reference state (not always valid).
	 *        This will also update the aggregated heat capacity.
	 * 
	 * @param[in] elements constituents of this mixture
	 * @param[in] hydro new density, momentum, energy
	 * @param[in,out] eos new values, temperature is updated
	 */
	virtual void InitTemperature(std::vector<EOSComponent*> elements, Field& hydro, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				const tw::Float nm1 = MassDensity(hydro,cell);
				const tw::Float IE1 = InternalEnergy(nm1,hydro,cell);
				const tw::Float epsvn = MixVibrationalEnergy(hydro,cell);
				const tw::Float nv = MixVibrationalStates(hydro,cell);

				tw::Float nmcv = 0.0;
				for (auto c : elements) {
					nmcv += c->HeatCapacity(c->mat.mass*hydro(cell,c->hidx.ni),IE1);
				}

				eos(cell,eidx.nmcv) = nmcv;
				eos(cell,eidx.T) = IE1 / (tw::small_pos + nmcv);
				eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/std::log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
			}
		}
	}
	/**
	 * @brief Update temperature using reference data (e.g. previous time level).
	 *        This will also update the aggregated heat capacity.
	 *        The default presumes polytropic ideal gases and ignores reference states.
	 * 
	 * @param[in] elements constituents of this mixture
	 * @param[in] hydroRef reference density, momentum, energy
	 * @param[in] hydro new density, momentum, energy
	 * @param[in] eosRef reference temperature, pressure, heat capacity
	 * @param[in,out] eos new values, temperature is updated
	 */
	virtual void UpdateTemperature(std::vector<EOSComponent*> elements, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				const tw::Float nm1 = MassDensity(hydro,cell);
				const tw::Float IE1 = InternalEnergy(nm1,hydro,cell);
				const tw::Float epsvn = MixVibrationalEnergy(hydro,cell);
				const tw::Float nv = MixVibrationalStates(hydro,cell);

				tw::Float nmcv = 0.0;
				for (auto c : elements) {
					nmcv += c->HeatCapacity(c->mat.mass*hydro(cell,c->hidx.ni),IE1);
				}

				eos(cell,eidx.nmcv) = nmcv;
				eos(cell,eidx.T) = IE1 / (tw::small_pos + nmcv);
				eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/std::log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
			}
		}
	}
};

export struct EOSIdealGasMix:EOSMixture
{
	EOSIdealGasMix(const std::string& name,MetricSpace *m,Task *tsk) : EOSMixture(name,m,tsk) {}
	/**
	 * @brief Update temperature assuming polytropic ideal gas mix.
	 *        This will also update the aggregated heat capacity.
	 * 
	 * @param[in] elements constituents of this mixture
	 * @param[in] hydroRef ignored
	 * @param[in] hydro new density, momentum, energy
	 * @param[in] eosRef ignored
	 * @param[in,out] eos new values, temperature is updated
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
	 * @brief Update temperature using reference data (e.g. previous time level).
	 *        This will also update the aggregated heat capacity.
	 * 
	 * @param[in] elements constituents of this mixture
	 * @param[in] hydroRef reference density, momentum, energy
	 * @param[in] hydro new density, momentum, energy
	 * @param[in] eosRef reference temperature, pressure, heat capacity
	 * @param[in,out] eos new values, temperature is updated
	 */
	virtual void UpdateTemperature(std::vector<EOSComponent*> elements, Field& hydroRef, Field& hydro, Field& eosRef, Field& eos)
	{
		#pragma omp parallel
		{
			for (auto cell : EntireCellRange(*space,1))
			{
				const tw::Float T0 = eosRef(cell,eidx.T);
				const tw::Float nm0 = MassDensity(hydroRef,cell);
				const tw::Float nm1 = MassDensity(hydro,cell);
				const tw::Float IE0 = InternalEnergy(nm0,hydroRef,cell);
				const tw::Float IE1 = InternalEnergy(nm1,hydro,cell);
				const tw::Float epsvn = MixVibrationalEnergy(hydro,cell);
				const tw::Float nv = MixVibrationalStates(hydro,cell);
				const tw::Float dE = std::copysign(IE0*tw::tiny,IE1 - IE0) + IE1 - IE0;
				const tw::Float small_nmcv = eosRef(cell,eidx.nmcv) * tw::tiny;

				auto nmcv = [dE,IE0,IE1,hydro,hydroRef,&cell,&elements] (tw::Float IE,tw::Float T) {
					tw::Float ans = 0.0;
					tw::Float k = (IE - IE0) / dE;
					for (auto c : elements) {
						tw::Float n0 = hydroRef(cell,c->hidx.ni);
						tw::Float n1 = hydro(cell,c->hidx.ni);
						tw::Float nm = c->mat.mass * (n0 + (n1 - n0)) * (std::isfinite(k) && k > 0.0 && k < 1.0 ? k : 1.0);
						ans += c->HeatCapacity(nm,IE);
					}
					return ans;
				};
				auto dTdE = [dE,nmcv,IE0,IE1,nm0,nm1,hydro,hydroRef,&cell,&elements] (tw::Float IE,tw::Float T) {
					// the following form is designed to approximate the polytropic result when cv is constant, in particular,
					// it comes from supposing dE/dt = d/dt(nmcvT)
					tw::Float cv = 2 * nmcv(0.5*(IE0+IE1),T) / (nm0 + nm1);
					tw::Float dnmcvdE = (nm1 - nm0) * cv / dE;
					return (1 - T*dnmcvdE)/(tw::small_pos + nmcv(IE,T));
				};

				// advance dy/dt = f(t,y) where y = T and t = E
				auto T1 = RK4Step<tw::Float>(T0, IE0, IE1 - IE0, dTdE);
				// if (cell.dcd1()==200) {
				// 	tw::Float poly = IE1 / nmcv(IE1,T1);
				// 	logger::WARN(std::format("T0 = {:.5} T1 = {:.8} poly = {:.8}",T0,T1,poly));
				// }
				eos(cell,eidx.nmcv) = nmcv(IE1,T1);
				eos(cell,eidx.T) = T1;
				eos(cell,eidx.Tv) = (epsvn/(nv+tw::small_pos))/std::log(1.0001 + epsvn/(hydro(cell,hidx.x)+tw::small_pos));
			}
		}
	}
};
