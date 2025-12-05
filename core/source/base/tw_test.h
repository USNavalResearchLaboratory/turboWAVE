// Macros for testing.
// When expanded these macros call into the `base/base` module.

#ifndef USE_STD_MODULE
    #include <functional>
#endif

#define ASSERT_TRUE(actual) assertTrue(actual,__FILE__,__LINE__,__func__,curr_test_name)
#define ASSERT_FALSE(actual) assertTrue(!actual,__FILE__,__LINE__,__func__,curr_test_name)
#define ASSERT_EQ(actual,expected) assertEqualInt(actual,expected,__FILE__,__LINE__,__func__,curr_test_name)
#define ASSERT_NEAR(actual,expected,tol) assertClose(actual,expected,tol,__FILE__,__LINE__,__func__,curr_test_name)
#define ASSERT_GTREQ(actual,expected) assertGtrEq(actual,expected,__FILE__,__LINE__,__func__,curr_test_name)
#define ASSERT_LESSEQ(actual,expected) assertLessEq(actual,expected,__FILE__,__LINE__,__func__,curr_test_name)

#define REGISTER(class_name,func_name) Register(std::bind(&class_name::func_name,this),#func_name)
