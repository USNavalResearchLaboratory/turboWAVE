module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_test.h"
#include "tw_logger.h"

export module hydro:group;
import :chemical;
import input;
import driver;
import fields;
import diagnostics;
import physics;
import injection;
import chemistry;
import parabolic;
import elliptic;
import logger;

export struct EquilibriumGroup:Driver
{
	std::vector<Chemical*> chemical; // explicitly typed submodule list
	std::shared_ptr<EOSMixture> eosMixData;
	bool mobile;

	// The hydro set contains indices into the state vector for this group.
	// The mass density index corresponds to the first chemical in the group.
	sparc::hydro_set hidx;
	// The eos set contains indices into the eos vector for this group.
	sparc::eos_set eidx;
	// The material data is packed and processed for optimization
	sparc::material_set matset;
	// Following is used to limit motion by zeroing forces on this group
	tw::Float forceFilter;

	tw::Float DensitySum(const Field& f,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=hidx.first;s<hidx.first+hidx.num;s++)
			ans += f(cell,s);
		return ans;
	}
	tw::Float DensityWeightedSum(const Field& f,std::valarray<tw::Float>& qty,const tw::cell& cell)
	{
		tw::Float ans = 0.0;
		for (tw::Int s=hidx.first;s<hidx.first+hidx.num;s++)
			ans += f(cell,s)*qty[s-hidx.first];
		return ans;
	}
	void LoadMassDensity(ScalarField& nm,const Field& f)
	{
		for (auto cell : EntireCellRange(*this,1))
			nm(cell) = DensityWeightedSum(f,matset.mass,cell);
	}
	void LoadMassDensityCv(ScalarField& nmcv,const Field& f)
	{
		for (auto cell : EntireCellRange(*this,1))
			nmcv(cell) = DensityWeightedSum(f,matset.cvm,cell);
	}
	tw::vec3 Velocity(const Field& f,const tw::cell& cell)
	{
		tw::Float nm = DensityWeightedSum(f,matset.mass,cell);
		return tw::vec3(f(cell,hidx.npx),f(cell,hidx.npy),f(cell,hidx.npz))/(tw::small_pos + nm);
	}
	void LoadVelocity(ScalarField& vel,const Field& f,tw::Int ax)
	{
		// assumes velocity components appear in order in state vector
		tw::Float nm;
		for (auto cell : EntireCellRange(*this,1))
		{
			nm = DensityWeightedSum(f,matset.mass,cell);
			vel(cell) = f(cell,hidx.npx+ax-1)/(tw::small_pos + nm);
		}
	}

	EquilibriumGroup(const std::string& name,MetricSpace *ms, Task *tsk) : Driver(name,ms,tsk)
	{
		if (native.unit_system!=tw::units::plasma)
			throw tw::FatalError("EquilibriumGroup module requires <native units = plasma>");

		mobile = true;
		forceFilter = 1.0;

		directives.Add("mobile",new tw::input::Bool(&mobile));
	}

	/// Before calling this on each group, the Hydro object needs to setup
	/// hidx and indexInState for all groups and chemicals.
	void SetupIndexing()
	{
		logger::TRACE(std::format("indexing for {}",name));
		matset.Allocate(chemical.size());
		logger::TRACE(std::format("add materials for {}",name));
		for (tw::Int i=0;i<chemical.size();i++)
			matset.AddMaterial(chemical[i]->mat,i);
		logger::TRACE(std::format("EOS indexing for {}",name));
		eosMixData->SetupIndexing(hidx,eidx,matset);
		// Now we can setup indexing for all Chemical modules
		for (auto chem : chemical) {
			auto ionizer = chem->SetupIndexing(hidx,eidx);
			if (ionizer.use_count() > 0)
			{
				// Setup the indexing for photoionization here (ionization tool cannot do it)
				Chemical *echem = (Chemical*)this->Root()->FindDriver(ionizer->electron_name,true);
				Chemical *ichem = (Chemical*)this->Root()->FindDriver(ionizer->ion_name,true);
				ionizer->hgas = hidx;
				ionizer->hgas.ni = chem->indexInState;
				ionizer->he = ((EquilibriumGroup*)echem->super)->hidx;
				ionizer->he.ni = echem->indexInState;
				ionizer->hi = ((EquilibriumGroup*)ichem->super)->hidx;
				ionizer->hi.ni = ichem->indexInState;
			}
		}
		logger::TRACE(std::format("finish index of chems in {}",name));
	}

	void VerifyInput()
	{
		Driver::VerifyInput();
		// Extract Chemical modules from the submodule list
		for (auto sub : sub_drivers) {
			Chemical *chem = dynamic_cast<Chemical*>(sub);
			if (chem!=NULL)
				chemical.push_back(chem);
			chem->VerifyInput();
		}
		// Find an EOSMixture
		for (auto tool : tools) {
			if (std::dynamic_pointer_cast<EOSMixture>(tool)) {
				eosMixData = std::dynamic_pointer_cast<EOSMixture>(tool);
			}
		}
		// If no EOSMixture create one
		if (!eosMixData) {
			auto new_tool = CreateTool("default_eos_mix",tw::tool_type::eosMixture);
			AddTool(new_tool);
			eosMixData = std::dynamic_pointer_cast<EOSMixture>(new_tool);
		}
	}

	bool GenerateFluid(Field& hydro,Field& eos,ScalarField& scratch,ScalarField& scratch2)
	{
		bool didGenerate = false;
		std::vector<bool> massLoaded(chemical.size());
		logger::TRACE("load totals");
		for (tw::Int i=0;i<chemical.size();i++)
			massLoaded[i] = chemical[i]->LoadFluid(hydro,hidx);
		logger::TRACE("load internal");
		for (tw::Int i=0;i<chemical.size();i++)
			if (massLoaded[i])
				chemical[i]->LoadInternalEnergy(hydro,eos,scratch,scratch2,eosMixData,hidx,eidx);
		for (auto loaded : massLoaded)
			didGenerate |= loaded;
		return didGenerate;
	}
};