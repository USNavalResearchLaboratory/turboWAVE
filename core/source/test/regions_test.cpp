module;

#include "tw_includes.h"
#include "tw_test.h"

module driver;
import base;
import metric_space;

void EntireRegion::SpotCheckInsideTest() {
    // should be true for any input
    ASSERT_TRUE(Inside(tw::vec4(1000,-5e5,4e-5,9e32)));
}

void RectRegion::SpotCheckInsideTest() {
    size = tw::vec3(1,20,2);
    ASSERT_FALSE(Inside(tw::vec4(0,0,-11,0.9)));
    ASSERT_TRUE(Inside(tw::vec4(0,0,-9,0.9)));
}

void PrismRegion::SpotCheckInsideTest() {
    size = tw::vec3(1,20,2);
    ASSERT_FALSE(Inside(tw::vec4(0,0.49,0,.011)));
    ASSERT_TRUE(Inside(tw::vec4(0,0.49,0,0.009)));
}

void CircRegion::SpotCheckInsideTest() {
    radius = 5.0;
    ASSERT_FALSE(Inside(tw::vec4(0,5,5,0)));
    ASSERT_TRUE(Inside(tw::vec4(0,3,3,0)));
}

void CylinderRegion::SpotCheckInsideTest() {
    radius = 1.0;
    length = 2.0;
    ASSERT_FALSE(Inside(tw::vec4(0,0.5,0.5,1.01)));
    ASSERT_TRUE(Inside(tw::vec4(0,0.5,0.5,0.99)));
}

void CylindricalShellRegion::SpotCheckInsideTest() {
    innerRadius = 1.0;
    outerRadius = 1.1;
    length = 2.0;
    ASSERT_FALSE(Inside(tw::vec4(0,0.707,0.706,0)));
    ASSERT_TRUE(Inside(tw::vec4(0,0.707,0.708,0)));
}

void RoundedCylinderRegion::SpotCheckInsideTest() {
    radius = 1.0;
    length = 2.0;
    ASSERT_FALSE(Inside(tw::vec4(0,0,0,2.01)));
    ASSERT_TRUE(Inside(tw::vec4(0,0,0,1.99)));
}

void EllipsoidRegion::SpotCheckInsideTest() {
    size = tw::vec3(1,2,3);
    ASSERT_FALSE(Inside(tw::vec4(0,0.6,0,0)));
    ASSERT_TRUE(Inside(tw::vec4(0,0,0.9,0)));
}

void BoxArrayRegion::SpotCheckInsideTest() {
    spacing = tw::vec3(0,0,10);
    size = tw::vec3(1,1,1);
    ASSERT_FALSE(Inside(tw::vec4(0,0,0,0.6)));
    ASSERT_TRUE(Inside(tw::vec4(0,0,0,30)));
}

void TorusRegion::SpotCheckInsideTest() {
    minorRadius = 1.0;
    majorRadius = 3.0;
    ASSERT_FALSE(Inside(tw::vec4(0,0,0,0)));
    ASSERT_TRUE(Inside(tw::vec4(0,3,0,0)));
}

void ConeRegion::SpotCheckInsideTest() {
    baseRadius = 2.0;
    tipRadius = 1.0;
    length = 2.0;
    ASSERT_FALSE(Inside(tw::vec4(0,0,0,1.01)));
    ASSERT_TRUE(Inside(tw::vec4(0,0.49,0,0.99)));
}

void TangentOgiveRegion::SpotCheckInsideTest() {
    tipRadius = 0.1;
    bodyRadius = 1.0;
    length = 2.0;
    ASSERT_FALSE(Inside(tw::vec4(0,0.0,0.0,-1.01)));
    ASSERT_TRUE(Inside(tw::vec4(0,0,0,0)));
}

void SimpleRegion::SpotCheckInsideTest() {
    primitive->SpotCheckInsideTest();
}

void UnionRegion::SpotCheckInsideTest() {
    auto r1 = std::make_shared<SimpleRegion>("c1",space,task,
        std::make_unique<CircRegion>("c1",space,task));
    auto r2 = std::make_shared<SimpleRegion>("c2",space,task,
        std::make_unique<CircRegion>("c2",space,task));
    auto u = std::make_unique<UnionRegion>("union",space,task);
    u->elements.push_back(r1);
    u->elements.push_back(r2);
    tw::Float dz = 0.9;
    r1->translation = tw::vec4(0,0,0,-dz);
    r2->translation = tw::vec4(0,0,0,dz);
    r1->Initialize();
    r2->Initialize();
    u->Initialize();
    tw::Float crossingPt = std::sqrt(1 - dz*dz);
    ASSERT_TRUE(u->Inside(tw::vec4(0,0,0,-dz),0));
    ASSERT_TRUE(u->Inside(tw::vec4(0,0,0,dz),0));
    ASSERT_TRUE(u->Inside(tw::vec4(0,crossingPt - .001,0,0),0));
    ASSERT_FALSE(u->Inside(tw::vec4(0,crossingPt + .001,0,0),0));
}

void IntersectionRegion::SpotCheckInsideTest() {
    auto r1 = std::make_shared<SimpleRegion>("c1",space,task,
        std::make_unique<CircRegion>("c1",space,task));
    auto r2 = std::make_shared<SimpleRegion>("c2",space,task,
        std::make_unique<CircRegion>("c2",space,task));
    auto u = std::make_unique<IntersectionRegion>("inter",space,task);
    u->elements.push_back(r1);
    u->elements.push_back(r2);
    tw::Float dz = 0.9;
    r1->translation = tw::vec4(0,0,0,-dz);
    r2->translation = tw::vec4(0,0,0,dz);
    r1->Initialize();
    r2->Initialize();
    u->Initialize();
    tw::Float crossingPt = std::sqrt(1 - dz*dz);
    ASSERT_FALSE(u->Inside(tw::vec4(0,0,0,-dz),0));
    ASSERT_FALSE(u->Inside(tw::vec4(0,0,0,dz),0));
    ASSERT_TRUE(u->Inside(tw::vec4(0,crossingPt - .001,0,0),0));
    ASSERT_FALSE(u->Inside(tw::vec4(0,crossingPt + .001,0,0),0));
}

void DifferenceRegion::SpotCheckInsideTest() {
    auto r1 = std::make_shared<SimpleRegion>("c1",space,task,
        std::make_unique<CircRegion>("c1",space,task));
    auto r2 = std::make_shared<SimpleRegion>("c2",space,task,
        std::make_unique<CircRegion>("c2",space,task));
    auto u = std::make_unique<DifferenceRegion>("diff",space,task);
    u->elements.push_back(r1);
    u->elements.push_back(r2);
    tw::Float dz = 0.9;
    r1->translation = tw::vec4(0,0,0,-dz);
    r2->translation = tw::vec4(0,0,0,dz);
    r1->Initialize();
    r2->Initialize();
    u->Initialize();
    tw::Float crossingPt = std::sqrt(1 - dz*dz);
    ASSERT_TRUE(u->Inside(tw::vec4(0,0,0,-dz),0));
    ASSERT_FALSE(u->Inside(tw::vec4(0,0,0,dz),0));
    ASSERT_FALSE(u->Inside(tw::vec4(0,crossingPt - .001,0,0),0));
    ASSERT_FALSE(u->Inside(tw::vec4(0,crossingPt + .001,0,0),0));
}
