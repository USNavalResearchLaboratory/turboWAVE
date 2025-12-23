module;

#include "tw_includes.h"

export module driver:primitive_region;
import :tool;

import metric_space;
import tensor;

export struct PrimitiveRegion : ComputeTool
{
	PrimitiveRegion(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk) {}
	virtual bool Inside(const tw::vec3& pos) const = 0;
    virtual std::array<tw::Float,6> Bounds() const = 0;
};

export struct EntireRegion : PrimitiveRegion
{
	EntireRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {}
	virtual bool Inside(const tw::vec3& pos) const {
		return true;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-tw::big_pos,tw::big_pos,-tw::big_pos,tw::big_pos,-tw::big_pos,tw::big_pos};
	}
};

export struct RectRegion : PrimitiveRegion
{
	tw::vec3 size;
	RectRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		directives.Add("size",new tw::input::Vec3(&size),true);
	}
	virtual bool Inside(const tw::vec3& p) const {
		return std::fabs(2*p.x)<size.x && std::fabs(2*p.y)<size.y && std::fabs(2*p.z)<size.z;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-size.x/2,size.x/2,-size.y/2,size.y/2,-size.z/2,size.z/2};
	}
};

export struct PrismRegion : RectRegion
{
	PrismRegion(const std::string& name,MetricSpace *ms,Task *tsk) : RectRegion(name,ms,tsk) {}
	virtual bool Inside(const tw::vec3& p) const {
		return std::fabs(2*p.x)<size.x && std::fabs(2*p.y)<size.y && std::fabs(2*p.z)<size.z*0.5*(size.x-2*p.x)/size.x;
	}
};

export struct CircRegion : PrimitiveRegion
{
	tw::Float radius;
	CircRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		radius = 1.0;
		directives.Add("radius",new tw::input::Float(&radius),true);
	}
	virtual bool Inside(const tw::vec3& p) const {
		return Norm(p)<sqr(radius);
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-radius,radius,-radius,radius,-radius,radius};
	}
};

export struct CylinderRegion : PrimitiveRegion
{
	tw::Float radius,length;
	CylinderRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		directives.Add("radius",new tw::input::Float(&radius),true);
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual bool Inside(const tw::vec3& p) const {
		return sqr(p.x) + sqr(p.y) < sqr(radius) && std::fabs(2*p.z) < length;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-radius,radius,-radius,radius,-length/2,length/2};
	}
};

export struct CylindricalShellRegion : PrimitiveRegion
{
	tw::Float innerRadius,outerRadius,length;
	CylindricalShellRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		directives.Add("inner radius",new tw::input::Float(&innerRadius),true);
		directives.Add("outer radius",new tw::input::Float(&outerRadius),true);
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual bool Inside(const tw::vec3& p) const {
		tw::Float rho;
		rho = std::sqrt(p.x*p.x + p.y*p.y);
		return rho < outerRadius && rho > innerRadius && std::fabs(2*p.z) < length;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-outerRadius,outerRadius,-outerRadius,outerRadius,-length/2,length/2};
	}
};

export struct RoundedCylinderRegion : CylinderRegion
{
	RoundedCylinderRegion(const std::string& name,MetricSpace *ms,Task *tsk) : CylinderRegion(name,ms,tsk) {}
	// TODO: the inherited Bounds function does not enclose everything
	virtual bool Inside(const tw::vec3& p) const {
		bool ans;
		ans = sqr(p.x) + sqr(p.y) < sqr(radius) && std::fabs(2*p.z) < length;
		ans = ans || Norm(p - tw::vec3(0,0,length/2)) < sqr(radius);
		ans = ans || Norm(p + tw::vec3(0,0,length/2)) < sqr(radius);
		return ans;
	}
};

export struct EllipsoidRegion : RectRegion
{
	EllipsoidRegion(const std::string& name,MetricSpace *ms,Task *tsk) : RectRegion(name,ms,tsk) {}
	virtual bool Inside(const tw::vec3& p) const {
		return sqr(2*p.x/size.x) + sqr(2*p.y/size.y) + sqr(2*p.z/size.z) < 1.0;
	}
};

// export struct TrueSphere : CircRegion
// {
// 	TrueSphere(const std::string& name,MetricSpace *ms,Task *tsk) : CircRegion(name,ms,tsk) {}
// 	virtual bool Inside(const tw::vec3& pos,int depth) const {
// 		tw::vec3 c_cart = translation;
// 		tw::vec3 p_cart = pos;
// 		space->CurvilinearToCartesian(&c_cart);
// 		space->CurvilinearToCartesian(&p_cart);
// 		if (moveWithWindow && depth==0) {
// 			space->ToStartingWindow(&p_cart);
// 		}
// 		return complement ^ (Norm(p_cart-c_cart)<sqr(rbox.x));
// 	}
// };

export struct BoxArrayRegion : PrimitiveRegion
{
	tw::vec3 size,spacing;

	BoxArrayRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk)	{
		directives.Add("size",new tw::input::Vec3(&size),true);
		directives.Add("spacing",new tw::input::Vec3(&spacing),true);
	}
	virtual bool Inside(const tw::vec3& pos) const {
		// this is an infinite array, but of course can be intersected with something else to produce a finite one
		auto p = pos + 0.5*size;
		tw::Int i = MyFloor(p.x/spacing.x);
		tw::Int j = MyFloor(p.y/spacing.y);
		tw::Int k = MyFloor(p.z/spacing.z);
		return p.x - tw::Float(i)*spacing.x < size.x && p.y - tw::Float(j)*spacing.y < size.y && p.z - tw::Float(k)*spacing.z < size.z;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-size.x/2,size.x/2,-size.y/2,size.y/2,-size.z/2,size.z/2};
	}
};

export struct TorusRegion : PrimitiveRegion
{
	tw::Float majorRadius,minorRadius,length;
	TorusRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		majorRadius = 1.0;
		minorRadius = 0.1;
		directives.Add("minor radius",new tw::input::Float(&minorRadius));
		directives.Add("major radius",new tw::input::Float(&majorRadius));
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual bool Inside(const tw::vec3& p) const {
		tw::Float rho;
		rho = std::sqrt(p.x*p.x + p.y*p.y);
		return rho > majorRadius-minorRadius && rho < majorRadius+minorRadius &&
			sqr(p.z) < sqr(minorRadius)-sqr(rho-majorRadius);
	}
    virtual std::array<tw::Float,6> Bounds() const {
		auto r = tw::vec3(majorRadius+minorRadius,majorRadius+minorRadius,minorRadius);
		return std::array<tw::Float,6> {-r.x,r.x,-r.y,r.y,-r.z,r.z};
	}
};

export struct ConeRegion : PrimitiveRegion
{
	tw::Float tipRadius,baseRadius,length;
	ConeRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		tipRadius = 0.1;
		baseRadius = 1.0;
		length = 1.0;
		directives.Add("tip radius",new tw::input::Float(&tipRadius));
		directives.Add("base radius",new tw::input::Float(&baseRadius));
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual bool Inside(const tw::vec3& p) const {
		tw::Float rho,rOfz;
		rho = std::sqrt(p.x*p.x + p.y*p.y);
		rOfz = tipRadius - (p.z - length/2)*(baseRadius-tipRadius)/length;
		return rho < rOfz && sqr(p.z) < sqr(length/2);
	}
    virtual std::array<tw::Float,6> Bounds() const {
		auto r = tw::vec3(baseRadius,baseRadius,length/2);
		return std::array<tw::Float,6> {-r.x,r.x,-r.y,r.y,-r.z,r.z};
	}
};

export struct TangentOgiveRegion : PrimitiveRegion
{
	tw::Float tipRadius,bodyRadius,length;
	TangentOgiveRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		directives.Add("tip radius",new tw::input::Float(&tipRadius),true);
		directives.Add("body radius",new tw::input::Float(&bodyRadius),true);
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual bool Inside(const tw::vec3& p) const {
		tw::Float rho,rOfz,ogiveRadius,x0,xt,yt;

		x0 = length - tipRadius;
		ogiveRadius = (x0*x0+sqr(bodyRadius)-sqr(tipRadius))/(2.0*(bodyRadius-tipRadius));
		yt = tipRadius*(ogiveRadius-bodyRadius)/(ogiveRadius-tipRadius);
		xt = x0 + std::sqrt(sqr(tipRadius)-yt*yt);

		rho = std::sqrt(p.x*p.x + p.y*p.y);
		auto z = p.z + length/2;
		rOfz = -1.0;
		if (z>0.0 && z<xt)
			rOfz = std::sqrt(sqr(ogiveRadius)-z*z) + bodyRadius - ogiveRadius;
		if (z>=xt && z<2.0*length/2)
			rOfz = std::sqrt(sqr(tipRadius)-sqr(z-x0));
		return rho < rOfz;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		auto r = tw::vec3(bodyRadius,bodyRadius,length/2);
		return std::array<tw::Float,6> {-r.x,r.x,-r.y,r.y,-r.z,r.z};
	}
};
