// macros for testing must be in a header file.
// this has to be included after base is imported.

#define ASSERT_EQ(actual,expected) assertEqualInt(actual,expected,__FILE__,__LINE__,__func__,testName)
#define ASSERT_NEAR(actual,expected,tol) assertClose(actual,expected,tol,__FILE__,__LINE__,__func__,testName)
#define ASSERT_GTREQ(actual,expected) assertGtrEq(actual,expected,__FILE__,__LINE__,__func__,testName)
#define ASSERT_LESSEQ(actual,expected) assertLessEq(actual,expected,__FILE__,__LINE__,__func__,testName)
#define REGISTER_TEST() testName = __func__

