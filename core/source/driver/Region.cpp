module;

#include <tree_sitter/api.h>
#include "tw_includes.h"

export module driver:region;
import :tool;

import input;
import metric_space;
import tensor;

enum class bool_op { Union, Intersection, Difference };

export struct Region : ComputeTool
{
	tw::vec3 rbox;
	tw::basis orientation;
	tw::vec3 translation;
	tw::vec3 euler;
	bool complement,moveWithWindow;
	/// The operations go before the corresponding parts.
	/// The operation before the first part follows an implicit entire region (so probably not a union).
	std::vector<bool_op> ops;
	std::vector<std::shared_ptr<Region>> composite;

	Region(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk) {
		complement = false;
		moveWithWindow = true;
		orientation.u = tw::vec3(1,0,0);
		orientation.v = tw::vec3(0,1,0);
		orientation.w = tw::vec3(0,0,1);
		directives.Add("translation",new tw::input::Vec3(&translation),false);
		directives.Add("euler angles",new tw::input::Vec3(&euler),false);
		directives.Add("move with window",new tw::input::Bool(&moveWithWindow),false);
		directives.Add("complement",new tw::input::Bool(&complement),false);
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		bool ans = true;
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);
		for (auto i=0; i<composite.size(); i++) {
			if (ops[i] == bool_op::Union) {
				ans = ans || composite[i]->Inside(p,ds);
			} else if (ops[i] == bool_op::Intersection) {
				ans = ans && composite[i]->Inside(p,ds);
			} else if (ops[i] == bool_op::Difference) {
				ans = ans && !composite[i]->Inside(p,ds);
			}
		}
		return complement ^ ans;
	}
	virtual void Initialize(const MetricSpace& ds,Task *tsk) {
		orientation.SetWithEulerAngles(euler.x, euler.y, euler.z);
	}
	virtual void Translate(const tw::vec3& dr)
	{
		translation += dr;
	}
	void GetBoxLim(tw::Float *low,tw::Float *high,tw::Int ax)
	{
		*low = translation[ax-1] - rbox[ax-1];
		*high = translation[ax-1] + rbox[ax-1];
	}
	/// Set up the index limits such that the region is in cells [xl,xh] * [yl,yh] * [zl,zh] (inclusive)
	void GetRawCellBounds(tw::Int rawBounds[6],const MetricSpace& ms) const {
		tw::Int dims[3];
		tw::Float lims[6];
		std::vector<tw::Float> temp;

		for (auto d=0;d<3;d++)
		{
			dims[d] = ms.Dim(d+1);
			lims[2*d] = translation[d] - rbox[d];
			lims[2*d+1] = translation[d] + rbox[d];
		}

		// std::lower_bound returns first element >= to the test data
		// if r0.x < xpos[0], xl = 0   ,	if r1.x < xpos[0], xh = -1
		// if r0.x > xpos[xN1], xl = xN1+1   ,   if r1.x > xpos[xN1], xh = xN1
		for (auto d=0;d<3;d++)
		{
			temp.resize(dims[d]+2);
			for (auto i=0;i<=dims[d]+1;i++)
				temp[i] = ms.X(i,d+1);
			rawBounds[2*d] = tw::Int(std::lower_bound(temp.begin(),temp.end(),lims[2*d]) - temp.begin());
			rawBounds[2*d+1] = tw::Int(std::lower_bound(temp.begin(),temp.end(),lims[2*d+1]) - temp.begin())-1;
		}
	}
	void GetLocalCellBounds(tw::Int loc[6],const MetricSpace& ms) const {
		tw::Int raw[6];
		GetRawCellBounds(raw,ms);
		for (auto d = 0; d<3; d++) {
			loc[2*d] = raw[2*d] < 1 ? 1 : raw[2*d];
			loc[2*d+1] = raw[2*d+1] > ms.Dim(d+1) ? ms.Dim(d+1) : raw[2*d+1];
		}
	}
	void GetGlobalCellBounds(tw::Int glob[6],const MetricSpace& ms,Task *tsk) const {
		tw::Int low,high;
		tw::Int loc[6];
		GetLocalCellBounds(loc,ms);
		for (auto d=0;d<3;d++)
		{
			if (loc[2*d]>=1 && loc[2*d]<=ms.Dim(d+1))
				low = ms.GlobalCellIndex(loc[2*d],d+1);
			else
				low = ms.GlobalDim(d+1);
			glob[2*d] = tsk->strip[d+1].GetMin(low);
		}

		for (auto d=0;d<3;d++)
		{
			if (loc[2*d+1]>=1 && loc[2*d+1]<=ms.Dim(d+1))
				high = ms.GlobalCellIndex(loc[2*d+1],d+1);
			else
				high = 1;
			glob[2*d+1] = tsk->strip[d+1].GetMax(high);
		}
	}

	/// @brief read all directives in the block
	/// @param curs can be on block or on first child of block
	/// @param src source document
	virtual void ReadInputFileBlock(TSTreeCursor *curs,const std::string& src)
	{
		if (tw::input::node_kind(curs)=="block") {
			ts_tree_cursor_goto_first_child(curs);
		}
		do
		{
			directives.ReadNext(curs,src);
			// following is unconditional, see comments up top
			auto curs1 = tw::input::Cursor(curs);
			ReadInputFileDirective(curs1.get(),src);
		} while (ts_tree_cursor_goto_next_sibling(curs));
		directives.ThrowErrorIfMissingKeys(name);
	}

	virtual void ReadCheckpoint(std::ifstream& inFile)
	{
		inFile.read((char *)&rbox,sizeof(rbox));
		inFile.read((char *)&translation,sizeof(translation));
		inFile.read((char *)&orientation,sizeof(orientation));
		inFile.read((char *)&moveWithWindow,sizeof(bool));
	}

	virtual void WriteCheckpoint(std::ofstream& outFile)
	{
		outFile << name << " ";
		outFile.write((char *)&rbox,sizeof(rbox));
		outFile.write((char *)&translation,sizeof(translation));
		outFile.write((char *)&orientation,sizeof(orientation));
		outFile.write((char *)&moveWithWindow,sizeof(bool));
	}
};

export typedef std::shared_ptr<Region> SharedRegion;

export struct EntireRegion:Region
{
	EntireRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {}
	virtual void Initialize(const MetricSpace& ds,Task *tsk) {
		Region::Initialize(ds,tsk);
		auto globalSize = ds.GlobalPhysicalSize();
		auto globalCorner = ds.GlobalCorner();
		for (auto i=0; i<4; i++) {
			rbox[i] = globalSize[i+1]/2;
			translation[i] = globalCorner[i+1] + rbox[i];
		}
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		return complement ^ true;
	}
};

export struct RectRegion:Region
{
	tw::Float bounds[6];
	RectRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {
		directives.Add("bounds",new tw::input::Numbers<tw::Float>(&bounds[0],6),true);
	}
	virtual void Initialize(const MetricSpace& ds,Task *tsk) {
		Region::Initialize(ds,tsk);
		rbox.x = std::fabs(bounds[1]-bounds[0])/2;
		rbox.y = std::fabs(bounds[3]-bounds[2])/2;
		rbox.z = std::fabs(bounds[5]-bounds[4])/2;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);
		return complement ^ (std::fabs(p.x)<rbox.x && std::fabs(p.y)<rbox.y && std::fabs(p.z)<rbox.z);
	}
};

export struct PrismRegion:RectRegion
{
	PrismRegion(const std::string& name,MetricSpace *ms,Task *tsk) : RectRegion(name,ms,tsk) {}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);
		return complement ^ (std::fabs(p.x)<rbox.x && std::fabs(p.y)<rbox.y && std::fabs(p.z)<rbox.z*0.5*(rbox.x-p.x)/rbox.x);
	}
};

export struct CircRegion:Region
{
	tw::Float radius;
	CircRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {
		directives.Add("radius",new tw::input::Float(&radius),true);
	}
	virtual void Initialize(const MetricSpace& ds,Task *tsk) {
		Region::Initialize(ds,tsk);
		rbox = radius;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		return complement ^ (Norm(pos-translation)<sqr(rbox.x));
	}
};

export struct CylinderRegion:Region
{
	tw::Float radius,length;
	CylinderRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {
		directives.Add("radius",new tw::input::Float(&radius),true);
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual void Initialize(const MetricSpace& ds,Task *tsk) {
		Region::Initialize(ds,tsk);
		rbox.x = rbox.y = radius;
		rbox.z = length/2;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);
		return complement ^ (sqr(p.x) + sqr(p.y) < sqr(rbox.x) && std::fabs(p.z) < rbox.z);
	}
};

export struct CylindricalShellRegion:Region
{
	tw::Float innerRadius,outerRadius,length;
	CylindricalShellRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {
		directives.Add("inner radius",new tw::input::Float(&innerRadius),true);
		directives.Add("outer radius",new tw::input::Float(&outerRadius),true);
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual void Initialize(const MetricSpace& ds,Task *tsk) {
		Region::Initialize(ds,tsk);
		rbox.x = rbox.y = outerRadius;
		rbox.z = length/2;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho;
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);
		rho = std::sqrt(p.x*p.x + p.y*p.y);
		return complement ^ ( rho < outerRadius && rho > innerRadius && sqr(p.z) < sqr(rbox.z));
	}
};

export struct RoundedCylinderRegion:CylinderRegion
{
	RoundedCylinderRegion(const std::string& name,MetricSpace *ms,Task *tsk) : CylinderRegion(name,ms,tsk) {}
	// TODO: the rbox in this case does not enclose everything
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		bool ans;
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);
		ans = sqr(p.x) + sqr(p.y) < sqr(rbox.x) && std::fabs(p.z) < rbox.z;
		ans = ans || Norm(p - tw::vec3(0,0,rbox.z)) < sqr(rbox.x);
		ans = ans || Norm(p + tw::vec3(0,0,rbox.z)) < sqr(rbox.x);
		return complement ^ ans;
	}
};

export struct EllipsoidRegion:RectRegion
{
	EllipsoidRegion(const std::string& name,MetricSpace *ms,Task *tsk) : RectRegion(name,ms,tsk) {}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);
		return complement ^ (sqr(p.x/rbox.x) + sqr(p.y/rbox.y) + sqr(p.z/rbox.z) < 1.0);
	}
};

export struct TrueSphere:CircRegion
{
	TrueSphere(const std::string& name,MetricSpace *ms,Task *tsk) : CircRegion(name,ms,tsk) {}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 c_cart = translation;
		tw::vec3 p_cart = pos;
		ds.CurvilinearToCartesian(&c_cart);
		ds.CurvilinearToCartesian(&p_cart);
		return complement ^ (Norm(p_cart-c_cart)<sqr(rbox.x));
	}
};

export struct BoxArrayRegion:Region
{
	tw::vec3 size,spacing;

	BoxArrayRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk)	{
		size = tw::vec3(1,1,1);
		spacing = tw::vec3(2,2,2);
		directives.Add("size",new tw::input::Vec3(&size));
		directives.Add("spacing",new tw::input::Vec3(&spacing));
	}
	virtual void Initialize(const MetricSpace& ds,Task *tsk) {
		Region::Initialize(ds,tsk);
		rbox = 0.5*size;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);
		p += 0.5*size;
		tw::Int i = MyFloor(p.x/spacing.x);
		tw::Int j = MyFloor(p.y/spacing.y);
		tw::Int k = MyFloor(p.z/spacing.z);
		bool ans = p.x - tw::Float(i)*spacing.x < size.x && p.y - tw::Float(j)*spacing.y < size.y && p.z - tw::Float(k)*spacing.z < size.z;
		return complement ^ ans;
	}
};

export struct TorusRegion:Region
{
	tw::Float majorRadius,minorRadius,length;
	TorusRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {
		majorRadius = 1.0;
		minorRadius = 0.1;
		directives.Add("minor radius",new tw::input::Float(&minorRadius));
		directives.Add("major radius",new tw::input::Float(&majorRadius));
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual void Initialize(const MetricSpace& ds,Task *tsk) {
		Region::Initialize(ds,tsk);
		rbox = tw::vec3(majorRadius+minorRadius,majorRadius+minorRadius,minorRadius);
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho;
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);
		rho = std::sqrt(p.x*p.x + p.y*p.y);
		return complement ^ ( rho > majorRadius-minorRadius && rho < majorRadius+minorRadius &&
			sqr(p.z) < sqr(minorRadius)-sqr(rho-majorRadius));
	}
};

export struct ConeRegion:Region
{
	tw::Float tipRadius,baseRadius,length;
	ConeRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {
		tipRadius = 0.1;
		baseRadius = 1.0;
		length = 1.0;
		directives.Add("tip radius",new tw::input::Float(&tipRadius));
		directives.Add("base radius",new tw::input::Float(&baseRadius));
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual void Initialize(const MetricSpace& ds,Task *tsk) {
		Region::Initialize(ds,tsk);
		rbox.x = rbox.y = baseRadius;
		rbox.z = length/2;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho,rOfz;
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);
		rho = std::sqrt(p.x*p.x + p.y*p.y);
		rOfz = tipRadius - (p.z - rbox.z)*(baseRadius-tipRadius)/(2.0*rbox.z);
		return complement ^ ( rho < rOfz && sqr(p.z) < sqr(rbox.z));
	}
};

export struct TangentOgiveRegion:Region
{
	tw::Float tipRadius,bodyRadius,length;
	TangentOgiveRegion(const std::string& name,MetricSpace *ms,Task *tsk) : Region(name,ms,tsk) {
		tipRadius = 1.0;
		bodyRadius = 5.0;
		length = 5.0;
		rbox = tw::vec3(bodyRadius,bodyRadius,tw::big_pos);
		directives.Add("tip radius",new tw::input::Float(&tipRadius));
		directives.Add("body radius",new tw::input::Float(&bodyRadius));
		directives.Add("length",new tw::input::Float(&length),true);
	}

	virtual void Initialize(const MetricSpace& ds,Task *tsk) {
		Region::Initialize(ds,tsk);
		rbox.x = rbox.y = bodyRadius;
		rbox.z = length/2;
	}
	virtual bool Inside(const tw::vec3& pos,const MetricSpace& ds) const
	{
		tw::Float rho,rOfz,ogiveRadius,x0,xt,yt;
		tw::vec3 p = pos - translation;
		orientation.ExpressInBasis(&p);

		x0 = 2.0*rbox.z - tipRadius;
		ogiveRadius = (x0*x0+sqr(bodyRadius)-sqr(tipRadius))/(2.0*(bodyRadius-tipRadius));
		yt = tipRadius*(ogiveRadius-bodyRadius)/(ogiveRadius-tipRadius);
		xt = x0 + std::sqrt(sqr(tipRadius)-yt*yt);

		rho = std::sqrt(p.x*p.x + p.y*p.y);
		p.z += rbox.z;
		rOfz = -1.0;
		if (p.z>0.0 && p.z<xt)
			rOfz = std::sqrt(sqr(ogiveRadius)-p.z*p.z) + bodyRadius - ogiveRadius;
		if (p.z>=xt && p.z<2.0*rbox.z)
			rOfz = std::sqrt(sqr(tipRadius)-sqr(p.z-x0));
		return complement ^ (rho < rOfz);
	}
};
