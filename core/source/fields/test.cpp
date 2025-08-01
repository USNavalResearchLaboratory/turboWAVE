module;

#include "tw_includes.h"
#include "tw_test.h"

export module fields:test;
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

export struct IteratorTest: ComputeTool {
    IteratorTest(const std::string& name,MetricSpace *ms,Task *tsk): ComputeTool(name,ms,tsk) {}
	virtual void RegisterTests() {
        REGISTER(IteratorTest,TestCellRef1Layer);
        REGISTER(IteratorTest,TestStripRef1Layer);
    }
    void TestCellRef1Layer() {
        auto ss = StaticSpace(
            tw::node5{1,8,4,4,1},
            tw::vec4(1,1,1,1),
            std_packing,
            tw::node4{0,1,1,1}
        );

        auto cell = tw::cell(ss,1,1,1,1);
        auto idx = cell.Index(0);
        ASSERT_EQ(idx,1+1*(4+2)+1*(4+2)*(4+2));
        ASSERT_EQ(cell.dcd1(),1);
        ASSERT_EQ(cell.dcd2(),1);
        ASSERT_EQ(cell.dcd3(),1);

        // test out of range but ignorable coord, and ghost cell on axis 3
        cell = tw::cell(ss,5,3,3,0);
        idx = cell.Index(0);
        ASSERT_EQ(idx,0+3*(4+2)+3*(4+2)*(4+2)+5*0);
    }
    void TestStripRef1Layer() {
        auto ss = StaticSpace(
            tw::node5{1,8,4,4,1},
            tw::vec4(1,1,1,1),
            std_packing,
            tw::node4{0,1,1,1}
        );

        auto strip = tw::strip(ss,3,0,tw::node4{1,2,3,1});
        auto idx = strip.Index(4,0);
        ASSERT_EQ(idx,4 + 3*(4+2) + 2*(4+2)*(4+2));
        ASSERT_EQ(strip.dcd1(4),2);
        ASSERT_EQ(strip.dcd2(4),3);
        ASSERT_EQ(strip.dcd3(4),4);
    }
};