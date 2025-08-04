module;

#include "tw_includes.h"
#include "tw_test.h"

export module metric_space_test;
import base;
import tw_iterator;
import compute_tool;
import static_space;
import metric_space;

export struct MetricSpaceTest: ComputeTool {
    MetricSpaceTest(const std::string& name,MetricSpace *ms,Task *tsk): ComputeTool(name,ms,tsk) {}
    virtual void RegisterTests() {
        REGISTER(MetricSpaceTest,Test);
    }
	void Test() {
        tw::Float pos = space->X(1,1);
        ASSERT_NEAR(pos, space->Corner()[1] + space->dx(1)/2, 1e-4);
    }
};
