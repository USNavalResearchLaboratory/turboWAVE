module;

#include <tree_sitter/api.h>
#include "tw_includes.h"
#include "tw_logger.h"

export module driver:region;
import :primitive_region;
import :engine;

import input;
import metric_space;
import tensor;
import logger;

export enum class bool_op { Union, Intersection, Difference };

export struct Region : Engine
{
	// Transformation parameters:
	// boost,origin,orientation,translation are all used during updates.
	// boost3 and euler are only used to capture user inputs.

	/// displace the body before the boost-rotation
	tw::vec4 origin;
	/// boost the body
	tw::vec4 boost;
	/// rotate the body
	tw::basis orientation;
	/// displace the body after the boost-rotation
	tw::vec4 translation;

	/// in the active view euler components get applied as new_vec = Rz(alpha)*Rx(beta)*Rz(gamma)*old_vec,
	/// i.e. the "last angle" gets applied first.
	tw::vec3 euler,boost3;

	bool complement,moveWithWindow;
	/// The operations go before the corresponding elements.
	/// The operation before the first part follows an implicit entire region (so probably not a union).
	std::vector<bool_op> ops;
	/// If a primitive region exists, it takes precedence
	std::unique_ptr<PrimitiveRegion> primitive;

	Region(const std::string& name,MetricSpace *ms,Task *tsk) : Engine(name,ms,tsk) {
		complement = false;
		moveWithWindow = true;
		orientation.u = tw::vec3(1,0,0);
		orientation.v = tw::vec3(0,1,0);
		orientation.w = tw::vec3(0,0,1);
		directives.Add("origin",new tw::input::Vec4(&origin),false);
		directives.Add("boost",new tw::input::Vec3(&boost3),false);
		directives.Add("euler angles",new tw::input::Vec3(&euler),false);
		directives.Add("translation",new tw::input::Vec4(&translation),false);
		directives.Add("move with window",new tw::input::Bool(&moveWithWindow),false);
		directives.Add("complement",new tw::input::Bool(&complement),false);
	}
	virtual void Initialize() {
		Engine::Initialize();
		orientation.SetWithEulerAngles(euler.x, euler.y, euler.z);
		boost = tw::vec4(std::sqrt(1+Norm(boost3)),boost3);
		logger::DEBUG(std::format("{}:\n{},{},{}\n{},{},{}\n{},{},{}",name,orientation.u.x,orientation.u.y,orientation.u.z,
			orientation.v.x,orientation.v.y,orientation.v.z,
			orientation.w.x,orientation.w.y,orientation.w.z));
	}
	/// Start with conditional Galilean translation to moving window, then
	/// translate-rotate-boost-translate from the simulation frame to the profile's frame.
	/// This is in the reverse order compared to the active view.
	/// If there is a boost it will often be an "unboost."
	void TransformPoint(tw::vec4 *pos, int depth) const {
		if (moveWithWindow && depth==0) {
			space->ToStartingWindow(pos);
		}
		*pos -= translation;
		orientation.ExpressInBasis(pos);
		pos->Boost(boost);
		*pos -= origin;
	}
	bool Inside(const tw::vec4& pos,int depth) const
	{
		auto p = pos;
		TransformPoint(&p,depth);
		if (primitive) {
			return complement ^ primitive->Inside(p);
		}
		bool ans = true;
		for (auto i=0; i<elements.size(); i++) {
			auto rgn = std::dynamic_pointer_cast<Region>(elements[i]);
			if (rgn) {
				if (ops[i] == bool_op::Union) {
					ans = ans || rgn->Inside(p,depth+1);
				} else if (ops[i] == bool_op::Intersection) {
					ans = ans && rgn->Inside(p,depth+1);
				} else if (ops[i] == bool_op::Difference) {
					ans = ans && !rgn->Inside(p,depth+1);
				}
			}
		}
		return complement ^ ans;
	}
	/// take boundaries of an aligned hull and form its 8 vertices for further transformation
	std::vector<tw::vec4> BoundingVertices(std::array<tw::Float,6> hull) const {
		std::vector<tw::vec4> ans;
		ans.push_back(tw::vec4(0,hull[0],hull[2],hull[4]));
		ans.push_back(tw::vec4(0,hull[0],hull[2],hull[5]));
		ans.push_back(tw::vec4(0,hull[0],hull[3],hull[4]));
		ans.push_back(tw::vec4(0,hull[0],hull[3],hull[5]));
		ans.push_back(tw::vec4(0,hull[1],hull[2],hull[4]));
		ans.push_back(tw::vec4(0,hull[1],hull[2],hull[5]));
		ans.push_back(tw::vec4(0,hull[1],hull[3],hull[4]));
		ans.push_back(tw::vec4(0,hull[1],hull[3],hull[5]));
		return ans;
	}
	/// take any 8 vertices and return an aligned hull, i.e., take a box that may have been rotated
	/// and form a box around it that is aligned to the coordinate system
	std::array<tw::Float,6> AlignedHull(const std::vector<tw::vec4>& vertices) const {
		std::array<tw::Float,6> ans {tw::big_pos,tw::big_neg,tw::big_pos,tw::big_neg,tw::big_pos,tw::big_neg};
		for (auto i=0;i<8;i++) {
			ans[0] = std::min(ans[0],vertices[i][1]);
			ans[1] = std::max(ans[1],vertices[i][1]);
			ans[2] = std::min(ans[2],vertices[i][2]);
			ans[3] = std::max(ans[3],vertices[i][2]);
			ans[4] = std::min(ans[4],vertices[i][3]);
			ans[5] = std::max(ans[5],vertices[i][3]);
		}
		return ans;
	}
    std::array<tw::Float,6> Bounds(int depth) const {
		std::array<tw::Float,6> entire {tw::big_neg,tw::big_pos,tw::big_neg,tw::big_pos,tw::big_neg,tw::big_pos};
		if (complement) {
			// complement of something bounded is something unbounded
			// (but what if this thing is unbounded?)
			return entire;
		}
		if (primitive) {
			auto v = BoundingVertices(primitive->Bounds());
			for (auto i=0;i<8;i++) {
				TransformPoint(&v[i],depth);
			}
			return AlignedHull(v);
		}
		auto ans = entire;
		for (auto i=0; i<elements.size(); i++) {
			auto rgn = std::dynamic_pointer_cast<Region>(elements[i]);
			if (rgn) {
				if (ops[i] == bool_op::Union) {
					auto curr = rgn->Bounds(depth);
					for (auto j=0;j<3;j++) {
						ans[2*j] = std::min(ans[2*j],curr[2*j]);
						ans[2*j+1] = std::max(ans[2*j+1],curr[2*j+1]);
					}
				} else if (ops[i] == bool_op::Intersection) {
					auto curr = rgn->Bounds(depth);
					for (auto j=0;j<3;j++) {
						ans[2*j] = std::max(ans[2*j],curr[2*j]);
						ans[2*j+1] = std::min(ans[2*j+1],curr[2*j+1]);
						// TODO: can we do something other than error out?
						if (ans[2*j] > ans[2*j+1]) {
							throw tw::FatalError("error intersecting hulls");
						}
					}
				} else if (ops[i] == bool_op::Difference) {
					auto curr = rgn->Bounds(depth);
					for (auto j=0;j<3;j++) {
						ans[2*j] = std::min(ans[2*j],curr[2*j]);
						ans[2*j+1] = std::max(ans[2*j+1],curr[2*j+1]);
					}
				}
			}
		}
		return ans;
	}
	/// Compute [xl,xh,yl,yh,zl,zh] such that this region just fills the box defined by [xl,xh] * [yl,yh] * [zl,zh] (inclusive).
	/// This returns dirty results if the region is partly or fully out of the domain, clean with `GetLocalCellBounds`.
	void GetRawCellBounds(tw::Int rawBounds[6]) const {
		tw::Int dims[3];
		std::vector<tw::Float> temp;

		auto lims = Bounds(0);
		for (auto d=0;d<3;d++) {
			dims[d] = space->Dim(d+1);
		}

		// std::lower_bound returns first element >= to the test data
		// if r0.x < xpos[0], xl = 0   ,	if r1.x < xpos[0], xh = -1
		// if r0.x > xpos[xN1], xl = xN1+1   ,   if r1.x > xpos[xN1], xh = xN1
		for (auto d=0;d<3;d++)
		{
			temp.resize(dims[d]+2);
			for (auto i=0;i<=dims[d]+1;i++)
				temp[i] = space->X(i,d+1);
			rawBounds[2*d] = tw::Int(std::lower_bound(temp.begin(),temp.end(),lims[2*d]) - temp.begin());
			rawBounds[2*d+1] = tw::Int(std::lower_bound(temp.begin(),temp.end(),lims[2*d+1]) - temp.begin())-1;
		}
	}
	/// This takes `GetRawCellBounds` and constrains it such that the low>=1 and the high<=dim.
	/// Then `GetGlobalCellBounds` has enough information to set up a min/max procedure.
	void GetLocalCellBounds(tw::Int loc[6]) const {
		tw::Int raw[6];
		GetRawCellBounds(raw);
		for (auto d = 0; d<3; d++) {
			loc[2*d] = raw[2*d] < 1 ? 1 : raw[2*d];
			loc[2*d+1] = raw[2*d+1] > space->Dim(d+1) ? space->Dim(d+1) : raw[2*d+1];
		}
	}
	/// Some diagnostics use this to step through an index space just within the region boundaries.
	/// This is an expensive call that has to search mesh points and do a remote reduction.
	void GetGlobalCellBounds(tw::Int glob[6]) const {
		tw::Int low,high;
		tw::Int loc[6];
		GetLocalCellBounds(loc);
		for (auto d=0;d<3;d++)
		{
			if (loc[2*d]>=1 && loc[2*d]<=space->Dim(d+1))
				low = space->GlobalCellIndex(loc[2*d],d+1);
			else
				low = space->GlobalDim(d+1);
			glob[2*d] = task->strip[d+1].GetMin(low);
		}

		for (auto d=0;d<3;d++)
		{
			if (loc[2*d+1]>=1 && loc[2*d+1]<=space->Dim(d+1))
				high = space->GlobalCellIndex(loc[2*d+1],d+1);
			else
				high = 1;
			glob[2*d+1] = task->strip[d+1].GetMax(high);
		}
	}
};

export struct SimpleRegion : Region {
	SimpleRegion(const std::string& name,MetricSpace *ms,Task *tsk,std::unique_ptr<PrimitiveRegion> primitive) : Region(name,ms,tsk) {
		this->primitive = std::move(primitive);
	}
	/// we need an override since this object works a little differently, it has to check for assignments
	/// to the primtive's variables, and disallows nesting, custom directives, etc..
	virtual void ReadInputFileBlock(TSTreeCursor *curs,const std::string& src)
	{
		if (tw::input::node_kind(curs)=="block") {
			if (ts_tree_cursor_goto_first_child(curs)) {
				if (tw::input::next_named_node(curs,true)) {
					do
					{
						if (tw::input::node_kind(curs)=="comment") {
							continue;
						}
						try {
							if (!directives.ReadNext(curs,src)) {
								tw::input::ThrowParsingError(curs,src,"unknown directive");
							}
						} catch (tw::FatalError) {
							if (!primitive->directives.ReadNext(curs,src)) {
								tw::input::ThrowParsingError(curs,src,"unknown directive");
							}
						}
					} while (tw::input::next_named_node(curs,false));
				}
			}
		}
		directives.ThrowErrorIfMissingKeys(name);
	}
};

export struct UnionRegion : Region {
	UnionRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {}
	virtual void Initialize() {
		Region::Initialize();
		ops.clear();
		ops.push_back(bool_op::Intersection);
		for (auto i=1;i<elements.size();i++) {
			ops.push_back(bool_op::Union);
		}
	}
};

export struct IntersectionRegion : Region {
	IntersectionRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {}
	virtual void Initialize() {
		Region::Initialize();
		ops.clear();
		ops.push_back(bool_op::Intersection);
		for (auto i=1;i<elements.size();i++) {
			ops.push_back(bool_op::Intersection);
		}
	}
};

export struct DifferenceRegion : Region {
	DifferenceRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {}
	virtual void Initialize() {
		Region::Initialize();
		ops.clear();
		ops.push_back(bool_op::Intersection);
		for (auto i=1;i<elements.size();i++) {
			ops.push_back(bool_op::Difference);
		}
	}
};
