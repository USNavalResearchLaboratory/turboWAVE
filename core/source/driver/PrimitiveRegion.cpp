module;

#include "tw_includes.h"
#include "tw_test.h"

export module driver:primitive_region;
import :tool;

import metric_space;
import tensor;

export struct PrimitiveRegion : ComputeTool
{
	PrimitiveRegion(const std::string& name,MetricSpace *ms,Task *tsk) : ComputeTool(name,ms,tsk) {}
	virtual bool Inside(const tw::vec4& pos) const = 0;
    virtual std::array<tw::Float,6> Bounds() const = 0;
	// This gets invoked by the owning SimpleRegion's test, because the factory is unaware of primitive regions.
	// We have to *require* these tests, or else the runner will report an empty test as success.
	virtual void SpotCheckInsideTest() = 0;
};

export struct EntireRegion : PrimitiveRegion
{
	EntireRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {}
	virtual bool Inside(const tw::vec4& pos) const {
		return true;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-tw::big_pos,tw::big_pos,-tw::big_pos,tw::big_pos,-tw::big_pos,tw::big_pos};
	}
	virtual void SpotCheckInsideTest();
};

export struct RectRegion : PrimitiveRegion
{
	tw::vec3 size;
	RectRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		directives.Add("size",new tw::input::Vec3(&size),true);
	}
	virtual bool Inside(const tw::vec4& p) const {
		return std::fabs(2*p[1])<size.x && std::fabs(2*p[2])<size.y && std::fabs(2*p[3])<size.z;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-size.x/2,size.x/2,-size.y/2,size.y/2,-size.z/2,size.z/2};
	}
	virtual void SpotCheckInsideTest();
};

export struct PrismRegion : RectRegion
{
	PrismRegion(const std::string& name,MetricSpace *ms,Task *tsk) : RectRegion(name,ms,tsk) {}
	virtual bool Inside(const tw::vec4& p) const {
		return std::fabs(2*p[1])<size.x && std::fabs(2*p[2])<size.y && std::fabs(2*p[3])<size.z*0.5*(size.x-2*p[1])/size.x;
	}
	virtual void SpotCheckInsideTest();
};

export struct CircRegion : PrimitiveRegion
{
	tw::Float radius;
	CircRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		radius = 1.0;
		directives.Add("radius",new tw::input::Float(&radius),true);
	}
	virtual bool Inside(const tw::vec4& p) const {
		return Norm(p)<sqr(radius);
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-radius,radius,-radius,radius,-radius,radius};
	}
	virtual void SpotCheckInsideTest();
};

export struct CylinderRegion : PrimitiveRegion
{
	tw::Float radius,length;
	CylinderRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		directives.Add("radius",new tw::input::Float(&radius),true);
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual bool Inside(const tw::vec4& p) const {
		return sqr(p[1]) + sqr(p[2]) < sqr(radius) && std::fabs(2*p[3]) < length;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-radius,radius,-radius,radius,-length/2,length/2};
	}
	virtual void SpotCheckInsideTest();
};

export struct CylindricalShellRegion : PrimitiveRegion
{
	tw::Float innerRadius,outerRadius,length;
	CylindricalShellRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		directives.Add("inner radius",new tw::input::Float(&innerRadius),true);
		directives.Add("outer radius",new tw::input::Float(&outerRadius),true);
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual bool Inside(const tw::vec4& p) const {
		tw::Float rho;
		rho = std::sqrt(p[1]*p[1] + p[2]*p[2]);
		return rho < outerRadius && rho > innerRadius && std::fabs(2*p[3]) < length;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-outerRadius,outerRadius,-outerRadius,outerRadius,-length/2,length/2};
	}
	virtual void SpotCheckInsideTest();
};

export struct RoundedCylinderRegion : CylinderRegion
{
	RoundedCylinderRegion(const std::string& name,MetricSpace *ms,Task *tsk) : CylinderRegion(name,ms,tsk) {}
	// TODO: the inherited Bounds function does not enclose everything
	virtual bool Inside(const tw::vec4& p) const {
		bool ans;
		ans = sqr(p[1]) + sqr(p[2]) < sqr(radius) && std::fabs(2*p[3]) < length;
		ans = ans || Norm(p - tw::vec4(0,0,0,length/2)) < sqr(radius);
		ans = ans || Norm(p + tw::vec4(0,0,0,length/2)) < sqr(radius);
		return ans;
	}
	virtual void SpotCheckInsideTest();
};

export struct EllipsoidRegion : RectRegion
{
	EllipsoidRegion(const std::string& name,MetricSpace *ms,Task *tsk) : RectRegion(name,ms,tsk) {}
	virtual bool Inside(const tw::vec4& p) const {
		return sqr(2*p[1]/size.x) + sqr(2*p[2]/size.y) + sqr(2*p[3]/size.z) < 1.0;
	}
	virtual void SpotCheckInsideTest();
};

// export struct TrueSphere : CircRegion
// {
// 	TrueSphere(const std::string& name,MetricSpace *ms,Task *tsk) : CircRegion(name,ms,tsk) {}
// 	virtual bool Inside(const tw::vec4& pos,int depth) const {
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
	virtual bool Inside(const tw::vec4& pos) const {
		// this is an infinite array, but of course can be intersected with something else to produce a finite one
		auto p = pos;
		p.Add3(0.5*size);
		tw::Int i = MyFloor(p[1]/spacing.x);
		tw::Int j = MyFloor(p[2]/spacing.y);
		tw::Int k = MyFloor(p[3]/spacing.z);
		return p[1] - tw::Float(i)*spacing.x < size.x && p[2] - tw::Float(j)*spacing.y < size.y && p[3] - tw::Float(k)*spacing.z < size.z;
	}
    virtual std::array<tw::Float,6> Bounds() const {
		return std::array<tw::Float,6> {-size.x/2,size.x/2,-size.y/2,size.y/2,-size.z/2,size.z/2};
	}
	virtual void SpotCheckInsideTest();
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
	virtual bool Inside(const tw::vec4& p) const {
		tw::Float rho;
		rho = std::sqrt(p[1]*p[1] + p[2]*p[2]);
		return rho > majorRadius-minorRadius && rho < majorRadius+minorRadius &&
			sqr(p[3]) < sqr(minorRadius)-sqr(rho-majorRadius);
	}
    virtual std::array<tw::Float,6> Bounds() const {
		auto r = tw::vec3(majorRadius+minorRadius,majorRadius+minorRadius,minorRadius);
		return std::array<tw::Float,6> {-r.x,r.x,-r.y,r.y,-r.z,r.z};
	}
	virtual void SpotCheckInsideTest();
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
	virtual bool Inside(const tw::vec4& p) const {
		tw::Float rho,rOfz;
		rho = std::sqrt(p[1]*p[1] + p[2]*p[2]);
		rOfz = tipRadius - (p[3] - length/2)*(baseRadius-tipRadius)/length;
		return rho < rOfz && sqr(p[3]) < sqr(length/2);
	}
    virtual std::array<tw::Float,6> Bounds() const {
		auto r = tw::vec3(baseRadius,baseRadius,length/2);
		return std::array<tw::Float,6> {-r.x,r.x,-r.y,r.y,-r.z,r.z};
	}
	virtual void SpotCheckInsideTest();
};

export struct TangentOgiveRegion : PrimitiveRegion
{
	tw::Float tipRadius,bodyRadius,length;
	TangentOgiveRegion(const std::string& name,MetricSpace *ms,Task *tsk) : PrimitiveRegion(name,ms,tsk) {
		directives.Add("tip radius",new tw::input::Float(&tipRadius),true);
		directives.Add("body radius",new tw::input::Float(&bodyRadius),true);
		directives.Add("length",new tw::input::Float(&length),true);
	}
	virtual bool Inside(const tw::vec4& p) const {
		tw::Float rho,rOfz,ogiveRadius,x0,xt,yt;

		x0 = length - tipRadius;
		ogiveRadius = (x0*x0+sqr(bodyRadius)-sqr(tipRadius))/(2.0*(bodyRadius-tipRadius));
		yt = tipRadius*(ogiveRadius-bodyRadius)/(ogiveRadius-tipRadius);
		xt = x0 + std::sqrt(sqr(tipRadius)-yt*yt);

		rho = std::sqrt(p[1]*p[1] + p[2]*p[2]);
		auto z = p[3] + length/2;
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
	virtual void SpotCheckInsideTest();
};
