module;

#include "tw_includes.h"
#include "tw_test.h"

export module root_finder_test;
import base;
import tw_iterator;
import driver;
import numerics;
import metric_space;

export struct RootFinderTest: ComputeTool {
    RootFinderTest(const std::string& name,MetricSpace *ms,Task *tsk): ComputeTool(name,ms,tsk) {}
    virtual void RegisterTests() {
        REGISTER(RootFinderTest,QuadraticRootTest);
    }
	void QuadraticRootTest() {
        auto f = [] (tw::Float x) {
            return x*x + 4*x + 4;
        };
        ASSERT_NEAR(SecantMethod(f,1,2),-2,1e-6);
    }
};
