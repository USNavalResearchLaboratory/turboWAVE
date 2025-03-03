// The root testing components have to remain as a header because
// of the reliance on macros.  This header must follow `tw_includes.h`.

struct Testable {
	std::string testName;
};

#define ASSERT_EQ(actual,expected) assertEqualInt(actual,expected,__FILE__,__LINE__,__func__,testName)
#define ASSERT_NEAR(actual,expected,tol) assertClose(actual,expected,tol,__FILE__,__LINE__,__func__,testName)
#define ASSERT_GTREQ(actual,expected) assertGtrEq(actual,expected,__FILE__,__LINE__,__func__,testName)
#define ASSERT_LESSEQ(actual,expected) assertLessEq(actual,expected,__FILE__,__LINE__,__func__,testName)
#define REGISTER_TEST() testName = __func__

inline void assertFailed(tw::Float actual,tw::Float expected,const std::string& expr,const std::string& file,int line,const std::string& func)
{
	std::ostringstream mess;
	mess << std::endl << term::err << " function " << term::red << func << term::reset_color << std::endl;
	mess << "  Assertion " << actual << " (actual) " << expr << " " << expected << " (expected) failed." << std::endl;
	mess << "  File: " << file <<  " , Line: " << line << std::endl;
	throw tw::FatalError(mess.str());
}

inline void assertClose(tw::Float actual,tw::Float expected,tw::Float tol,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (fabs(actual-expected)>tol)
		assertFailed(actual,expected,"~",file,line,func);
}

inline void assertEqualInt(tw::Int actual,tw::Int expected,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (actual!=expected)
		assertFailed(actual,expected,"==",file,line,func);
}

inline void assertGtrEq(tw::Int actual,tw::Int expected,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (actual<expected)
		assertFailed(actual,expected,">=",file,line,func);
}

inline void assertLessEq(tw::Int actual,tw::Int expected,const std::string& file,int line,const std::string& func,std::string& testName)
{
	testName = func;
	if (actual>expected)
		assertFailed(actual,expected,"<=",file,line,func);
}

