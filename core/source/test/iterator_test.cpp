module;

#include "tw_includes.h"
#include "tw_test.h"
#ifndef USE_STD_MODULE
    #include <memory>
#endif

export module iterator_test;
import base;
import tw_iterator;
import compute_tool;
import static_space;
import metric_space;

std::unique_ptr<StaticSpace> get_space(tw::Int layers) {
    return std::make_unique<StaticSpace>(StaticSpace(
            tw::node5{1,8,4,4,1},
            tw::vec4(1,1,1,1),
            std_packing,
            tw::node4{0,layers,layers,layers}
        )
    );
}

export struct IteratorTest: ComputeTool {
    IteratorTest(const std::string& name,MetricSpace *ms,Task *tsk): ComputeTool(name,ms,tsk) {}
	virtual void RegisterTests() {
        REGISTER(IteratorTest,TestCellRef1Layer);
        REGISTER(IteratorTest,TestStripRef1Layer);
        REGISTER(IteratorTest, TestVecStripRef1Layer);
        REGISTER(IteratorTest,TestCellRef2Layer);
        REGISTER(IteratorTest,TestStripRef2Layer);
        REGISTER(IteratorTest, TestVecStripRef2Layer);
        REGISTER(IteratorTest,TestCellRange);
        REGISTER(IteratorTest,TestStripRange);
    }
    void TestCellRef1Layer() {
        auto ss = get_space(1);
        auto cell = tw::cell(*ss,1,1,1,1);
        auto idx = cell.Index(0);
        ASSERT_EQ(idx,(1-0)*0 + (1-0)*(4+2)*(4+2) + (1-0)*(4+2) + (1-0));
        ASSERT_EQ(cell.dcd1(),1);
        ASSERT_EQ(cell.dcd2(),1);
        ASSERT_EQ(cell.dcd3(),1);

        // test out of range but ignorable coord, and ghost cell on axis 3
        cell = tw::cell(*ss,5,3,3,0);
        idx = cell.Index(0);
        ASSERT_EQ(idx,(5-0)*0 + (3-0)*(4+2)*(4+2) + (3-0)*(4+2) + (0-0));
    }
    void TestStripRef1Layer() {
        auto ss = get_space(1);
        auto strip = tw::strip(*ss,3,0,tw::node4{1,2,3,1});
        auto idx = strip.Index(4,0);
        ASSERT_EQ(idx,(1-0)*0 + (2-0)*(4+2)*(4+2) + (3-0)*(4+2) + (4-0));
        ASSERT_EQ(strip.dcd1(4),2);
        ASSERT_EQ(strip.dcd2(4),3);
        ASSERT_EQ(strip.dcd3(4),4);
    }
    void TestVecStripRef1Layer() {
        auto ss = get_space(1);
        auto strip = tw::xstrip<3>(*ss,0,tw::node4{1,2,3,1});
        auto idx = strip.Index(4,0);
        ASSERT_EQ(idx,(1-0)*0 + (2-0)*(4+2)*(4+2) + (3-0)*(4+2) + (4-0));
        ASSERT_EQ(strip.dcd1(4),2);
        ASSERT_EQ(strip.dcd2(4),3);
        ASSERT_EQ(strip.dcd3(4),4);
    }

    void TestCellRef2Layer() {
        auto ss = get_space(2);
        auto cell = tw::cell(*ss,1,1,1,1);
        auto idx = cell.Index(0);
        ASSERT_EQ(idx,(1-(-1))*0 + (1-(-1))*(4+4)*(4+4) + (1-(-1))*(4+4) + (1-(-1)));
        ASSERT_EQ(cell.dcd1(),1);
        ASSERT_EQ(cell.dcd2(),1);
        ASSERT_EQ(cell.dcd3(),1);

        // test out of range but ignorable coord, and ghost cell on axis 3
        cell = tw::cell(*ss,5,3,3,0);
        idx = cell.Index(0);
        ASSERT_EQ(idx,(5-(-1))*0 + (3-(-1))*(4+4)*(4+4) + (3-(-1))*(4+4) + (0-(-1)));
    }
    void TestStripRef2Layer() {
        auto ss = get_space(2);
        auto strip = tw::strip(*ss,3,0,tw::node4{1,2,3,1});
        auto idx = strip.Index(4,0);
        ASSERT_EQ(idx,(1-(-1))*0 + (2-(-1))*(4+4)*(4+4) + (3-(-1))*(4+4) + (4-(-1)));
        ASSERT_EQ(strip.dcd1(4),2);
        ASSERT_EQ(strip.dcd2(4),3);
        ASSERT_EQ(strip.dcd3(4),4);
    }
    void TestVecStripRef2Layer() {
        auto ss = get_space(2);
        auto strip = tw::xstrip<3>(*ss,0,tw::node4{1,2,3,1});
        auto idx = strip.Index(4,0);
        ASSERT_EQ(idx,(1-(-1))*0 + (2-(-1))*(4+4)*(4+4) + (3-(-1))*(4+4) + (4-(-1)));
        ASSERT_EQ(strip.dcd1(4),2);
        ASSERT_EQ(strip.dcd2(4),3);
        ASSERT_EQ(strip.dcd3(4),4);
    }

    void TestCellRange() {
        auto ss = get_space(2);
        auto rng = InteriorCellRange(*ss,0);
        ASSERT_EQ(rng.GetReference(0).Index(0),tw::cell(*ss,1,1,1,1).Index(0));
        ASSERT_EQ(rng.GetReference(8).Index(0),tw::cell(*ss,1,1,3,1).Index(0));
        ASSERT_EQ(rng.GetReference(17).Index(0),tw::cell(*ss,1,2,1,2).Index(0));
    }
    void TestStripRange() {
        // For these tests we need to recall the access order of strips is not guaranteed to be stable
        // from one version to the next; as of this writing for z-strips x is the fast varying axis.
        auto ss = get_space(2);
        auto rng = StripRange(*ss,3,0,1,strongbool::no);
        ASSERT_EQ(rng.GetReference(0).Index(3,0),
            tw::strip(*ss,3,0,tw::node4{1,1,1,1})
            .Index(3,0));
        ASSERT_EQ(rng.GetReference(8).Index(3,0),
            tw::strip(*ss,3,0,tw::node4{1,1,2,1})
            .Index(3,0));
        ASSERT_EQ(rng.GetReference(17).Index(3,0),
            tw::strip(*ss,3,0,tw::node4{1,2,3,1})
            .Index(3,0));
    }
};