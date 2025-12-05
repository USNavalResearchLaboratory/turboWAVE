module;

#include "tw_includes.h"
#include "tw_test.h"

/// Test of turboWAVE regions.
export module regions_test;
import base;
import metric_space;
import driver;

export struct RegionsTest: ComputeTool {
    RegionsTest(const std::string& name,MetricSpace *ms,Task *tsk): ComputeTool(name,ms,tsk) {}
    virtual void RegisterTests() {
        REGISTER(RegionsTest,SingleContainerTest);
        REGISTER(RegionsTest,UnionTest);
        REGISTER(RegionsTest,IntersectionTest);
        REGISTER(RegionsTest,DifferenceTest);
    }
	void SingleContainerTest() {
        auto rgn = new RectRegion("r1",space,task);
        rgn->bounds = {1,2,-1,1,-1,1};
        rgn->Initialize();
        ASSERT_FALSE(rgn->Inside(tw::vec3(0,0,0),0));
        ASSERT_TRUE(rgn->Inside(tw::vec3(1.5,0,0),0));
        delete rgn;
    }
    void UnionTest() {
        auto r1 = std::make_shared<CircRegion>("c1",space,task);
        auto r2 = std::make_shared<CircRegion>("c2",space,task);
        auto u = new Region("union",space,task);
        tw::Float dz = 0.9;
        r1->radius = 1.0;
        r2->radius = 1.0;
        r1->translation = tw::vec3(0,0,-dz);
        r2->translation = tw::vec3(0,0,dz);
        u->ops.push_back(bool_op::Intersection); // first is with the implicit EntireRegion
        u->composite.push_back(r1);
        u->ops.push_back(bool_op::Union);
        u->composite.push_back(r2);
        r1->Initialize();
        r2->Initialize();
        u->Initialize();
        tw::Float crossingPt = std::sqrt(1 - dz*dz);
        ASSERT_TRUE(u->Inside(tw::vec3(0,0,-dz),0));
        ASSERT_TRUE(u->Inside(tw::vec3(0,0,dz),0));
        ASSERT_TRUE(u->Inside(tw::vec3(crossingPt - .001,0,0),0));
        ASSERT_FALSE(u->Inside(tw::vec3(crossingPt + .001,0,0),0));
    }
	void IntersectionTest() {
        auto r1 = std::make_shared<CircRegion>("c1",space,task);
        auto r2 = std::make_shared<CircRegion>("c2",space,task);
        auto u = new Region("intersection",space,task);
        tw::Float dz = 0.9;
        r1->radius = 1.0;
        r2->radius = 1.0;
        r1->translation = tw::vec3(0,0,-dz);
        r2->translation = tw::vec3(0,0,dz);
        u->ops.push_back(bool_op::Intersection); // first is with the implicit EntireRegion
        u->composite.push_back(r1);
        u->ops.push_back(bool_op::Intersection);
        u->composite.push_back(r2);
        r1->Initialize();
        r2->Initialize();
        u->Initialize();
        tw::Float crossingPt = std::sqrt(1 - dz*dz);
        ASSERT_FALSE(u->Inside(tw::vec3(0,0,-dz),0));
        ASSERT_FALSE(u->Inside(tw::vec3(0,0,dz),0));
        ASSERT_TRUE(u->Inside(tw::vec3(crossingPt - .001,0,0),0));
        ASSERT_FALSE(u->Inside(tw::vec3(crossingPt + .001,0,0),0));
    }
	void DifferenceTest() {
        auto r1 = std::make_shared<CircRegion>("c1",space,task);
        auto r2 = std::make_shared<CircRegion>("c2",space,task);
        auto u = new Region("diff",space,task);
        tw::Float dz = 0.9;
        r1->radius = 1.0;
        r2->radius = 1.0;
        r1->translation = tw::vec3(0,0,-dz);
        r2->translation = tw::vec3(0,0,dz);
        u->ops.push_back(bool_op::Intersection); // first is with the implicit EntireRegion
        u->composite.push_back(r1);
        u->ops.push_back(bool_op::Difference);
        u->composite.push_back(r2);
        r1->Initialize();
        r2->Initialize();
        u->Initialize();
        tw::Float crossingPt = std::sqrt(1 - dz*dz);
        ASSERT_TRUE(u->Inside(tw::vec3(0,0,-dz),0));
        ASSERT_FALSE(u->Inside(tw::vec3(0,0,dz),0));
        ASSERT_FALSE(u->Inside(tw::vec3(crossingPt - .001,0,0),0));
        ASSERT_FALSE(u->Inside(tw::vec3(crossingPt + .001,0,0),0));
    }
};
