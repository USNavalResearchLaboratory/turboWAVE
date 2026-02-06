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
        auto prim = std::make_unique<RectRegion>("r1",space,task);
        prim->size = tw::vec3(1,2,2);
        auto rgn = std::make_unique<SimpleRegion>("r1",space,task,std::move(prim));
        rgn->translation = tw::vec3(1.5,0,0);
        rgn->Initialize();
        ASSERT_FALSE(rgn->Inside(tw::vec4(0,0,0,0),0));
        ASSERT_TRUE(rgn->Inside(tw::vec4(0,1.5,0,0),0));
    }
    void UnionTest() {
        auto r1 = std::make_shared<SimpleRegion>("c1",space,task,
            std::make_unique<CircRegion>("c1",space,task));
        auto r2 = std::make_shared<SimpleRegion>("c2",space,task,
            std::make_unique<CircRegion>("c2",space,task));
        auto u = std::make_unique<UnionRegion>("union",space,task);
        u->elements.push_back(r1);
        u->elements.push_back(r2);
        tw::Float dz = 0.9;
        r1->translation = tw::vec3(0,0,-dz);
        r2->translation = tw::vec3(0,0,dz);
        r1->Initialize();
        r2->Initialize();
        u->Initialize();
        tw::Float crossingPt = std::sqrt(1 - dz*dz);
        ASSERT_TRUE(u->Inside(tw::vec4(0,0,0,-dz),0));
        ASSERT_TRUE(u->Inside(tw::vec4(0,0,0,dz),0));
        ASSERT_TRUE(u->Inside(tw::vec4(0,crossingPt - .001,0,0),0));
        ASSERT_FALSE(u->Inside(tw::vec4(0,crossingPt + .001,0,0),0));
    }
	void IntersectionTest() {
        auto r1 = std::make_shared<SimpleRegion>("c1",space,task,
            std::make_unique<CircRegion>("c1",space,task));
        auto r2 = std::make_shared<SimpleRegion>("c2",space,task,
            std::make_unique<CircRegion>("c2",space,task));
        auto u = std::make_unique<IntersectionRegion>("inter",space,task);
        u->elements.push_back(r1);
        u->elements.push_back(r2);
        tw::Float dz = 0.9;
        r1->translation = tw::vec3(0,0,-dz);
        r2->translation = tw::vec3(0,0,dz);
        r1->Initialize();
        r2->Initialize();
        u->Initialize();
        tw::Float crossingPt = std::sqrt(1 - dz*dz);
        ASSERT_FALSE(u->Inside(tw::vec4(0,0,0,-dz),0));
        ASSERT_FALSE(u->Inside(tw::vec4(0,0,0,dz),0));
        ASSERT_TRUE(u->Inside(tw::vec4(0,crossingPt - .001,0,0),0));
        ASSERT_FALSE(u->Inside(tw::vec4(0,crossingPt + .001,0,0),0));
    }
	void DifferenceTest() {
        auto r1 = std::make_shared<SimpleRegion>("c1",space,task,
            std::make_unique<CircRegion>("c1",space,task));
        auto r2 = std::make_shared<SimpleRegion>("c2",space,task,
            std::make_unique<CircRegion>("c2",space,task));
        auto u = std::make_unique<DifferenceRegion>("diff",space,task);
        u->elements.push_back(r1);
        u->elements.push_back(r2);
        tw::Float dz = 0.9;
        r1->translation = tw::vec3(0,0,-dz);
        r2->translation = tw::vec3(0,0,dz);
        r1->Initialize();
        r2->Initialize();
        u->Initialize();
        tw::Float crossingPt = std::sqrt(1 - dz*dz);
        ASSERT_TRUE(u->Inside(tw::vec4(0,0,0,-dz),0));
        ASSERT_FALSE(u->Inside(tw::vec4(0,0,0,dz),0));
        ASSERT_FALSE(u->Inside(tw::vec4(0,crossingPt - .001,0,0),0));
        ASSERT_FALSE(u->Inside(tw::vec4(0,crossingPt + .001,0,0),0));
    }
};
