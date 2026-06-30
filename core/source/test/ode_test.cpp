module;

#include "tw_includes.h"
#include "tw_test.h"

export module ode_test;
import base;
import tw_iterator;
import driver;
import numerics;
import metric_space;

export struct ODETest: ComputeTool {
    ODETest(const std::string& name,MetricSpace *ms,Task *tsk): ComputeTool(name,ms,tsk) {}
    virtual void RegisterTests() {
        REGISTER(ODETest,AdaptiveRK4Test);
        REGISTER(ODETest,BackwardTest);
    }
	void AdaptiveRK4Test() {
        // integrate dy/dt = sin(t)
        auto f = [] (tw::Float t,tw::Float y) {
            return std::sin(t);
        };
        ASSERT_NEAR(RK4Integrate<tw::Float>(0,0,100,1,f,1e-4),1-std::cos(100),1e-8);
    }
    void BackwardTest() {
        // integrate dy/dt = sin(t) backwards
        auto f = [] (tw::Float t,tw::Float y) {
            return std::sin(t);
        };
        ASSERT_NEAR(RK4Integrate<tw::Float>(1-std::cos(100),100,0,1,f,1e-4),0,1e-8);
    }
};
