Writing Unit Tests
==================

The unit testing framework is custom, but is influenced by `Google Test <https://google.github.io/googletest/>`_.

Other Testing
-------------

This document concerns creating unit tests.  Other forms of testing do not require any programming.  For integration tests, one merely adds a special comment to certain input files, see :doc:`run-testing` for discussion.  One should also carry out spot checks using the `Valgrind Memcheck tool <https://valgrind.org/info/tools.html>`_.  TurboWAVE is supposed to be *valgrind clean*, i.e., Memcheck should report *zero* memory errors or leaks originating from turboWAVE itself (there may be errors in external libraries we cannot control).

Unit Test Runner
----------------

The unit test runner is a method of ``Simulation``.  This is what is invoked when you type ``tw3d -n 2 --unit-test --all``.  This same shell command is issued by the python runner ``twtest`` in the case where unit tests are requested.  Ordinarily, there is no need to modify the runner.

Test Suites
-----------

*Test suites* are created by overriding the ``RegisterTests`` method of ``ComputeTool`` or ``Module``.  As a corollary, there is at most one test suite for each derivative of either object.  It is important to note that nothing prevents you from deriving a ``ComputeTool`` for the sole purpose of testing things.

.. highlight:: c++

A typical ``RegisterTest`` function might look like this::

    #include "tw_test.h"
    virtual void RegisterTests() {
        REGISTER(MyTestTool,Test1);
        REGISTER(MyTestTool,Test2);
    }

The name `MyTestTool` is the actual class name of the ``ComputeTool`` or ``Module``.  The functions ``Test1`` and ``Test2`` are methods thereof wherein you program your tests. It is important to understand that a new instance of the object is created for every test case. The function names will show up in the test outputs, so it is helpful if they are descriptive.

There is no requirement on the filename or path wherein the test functions appear, these are simply methods of the ``ComputeTool`` and ``Module`` objects.  However, it is conventional to put all the test functions in a file with postfix ``_test.cpp``.  You are free to define any initialization, cleanup, or other functions for use or re-use by the test cases.

Test Cases
----------

Upon entry, test case functions have access to the grid and to MPI (see :ref:`testenvironment`) through the usual references.  On the other hand, no post-construction initialization is performed on the tested object, i.e., only the object's constructor is invoked.  Therefore one has to understand the object being tested well enough to put in place any dependencies it may have on other objects.  One way to learn how to do this is to look at existing test cases.  These can be discovered by looking in files with the ``_Test.cpp`` postfix.

The content of test cases is varied, but all test cases should use the test assertion macros so that the runner can properly analyze the test results.  As of this writing the available macros are::

    REGISTER_TEST(); // optional, helps runner to find the name of the test
    ASSERT_EQ(actual,expected); // integers are equal
    ASSERT_NEAR(actual,expected,tol); // floating point difference is within tolerance
    ASSERT_GTREQ(actual,expected); // actual integer >= expected integer
    ASSERT_LESSEQ(actual,expected); // actual integer <= expected integer

Parallel Testing
;;;;;;;;;;;;;;;;

Test cases are parallel programs, so generally any assertions will be run in parallel on all the distributed compute nodes.  In some cases it may be useful to only make assertions on a particular node.  This presents no problems, as long as one does not forget to run the ``REGISTER`` macro on *all* nodes.

Example
-------

Suppose we have a ``ComputeTool`` that adds two numbers::

    tw::Float AddTool::Add(tw::Float x,tw::Float y)
    {
        return x + y;
    }

We want to test commutativity and associativity of addition.  The registration is::

    virtual void RegisterTests() {
        REGISTER(AddTool,CommutativityTest);
        REGISTER(AddTool,AssociativityTest);
    }

And the test cases are::

    void AddTool::CommutativityTest()
    {
        tw::Float x=1.0,y=2.0;
        ASSERT_NEAR(Add(x,y),Add(y,x),1e-6);
    }
    void AddTool::AssociativityTest()
    {
        tw::Float x=1.0,y=2.0,z=3.0;
        ASSERT_NEAR(Add(Add(x,y),z),Add(x,Add(y,z)),1e-6);
    }

Polymorphism in Tests
---------------------

Making ``Test`` a virtual function defined on framework objects has both benefits and pitfalls.  The benefits are

    * There is a meaningful and easily computed metric of test coverage, i.e., the number of tests performed by each object.
    * Test outputs provide an idea of tests that are missing, every time you run the tests.
    * Much of the initialization of tests can be handled directly by the test runner.

The pitfalls are

    * If you create a test for an object with child types, the child types will run the same test, unless you explicitly override it.

.. _testenvironment:

Test Environment
----------------

The ``Test`` function is called from within a full turboWAVE simulation environment, i.e., a grid and domain decomposition are already in place by the time ``Test`` is called.  Objects are allowed to throw an error if they are incompatible with the environment that creates them.  In this case the test runner will catch the error and issue a warning that the test could not be carried out.

.. note::

    The test environment is something like a universal "fixture" in `Google Test <https://google.github.io/googletest/>`_.  Locally applied fixtures would correspond to functions defined on specific ``Module`` or ``ComputeTool`` subclasses that are re-used by the various test cases.

Optional Grid Control
;;;;;;;;;;;;;;;;;;;;;

As of this writing, the domain decomposition for all tests is fixed as :math:`1\times 1\times 2`, but the set of grids used for the testing can be controlled for each test suite.  The test grid is controlled by a static function of either ``ComputeTool`` or ``Module``::

	static bool SetTestEnvironment(tw::tool_type theType,tw::Int enviro,MetricSpace *ms,Task *tsk);
	static bool SetTestEnvironment(tw::module_type theType,tw::Int enviro,Simulation* sim);

These functions switch on the first argument, and create a grid that may depend on ``enviro``.  In order to control the grids that are used with a given test suite, cases in the switch must be modified.  The test runner will always start with ``enviro=1``, incrementing by 1.  Note that every test case in the suite will be called with every grid variant.  The individual test cases are free to do whatever they wish with a given grid, including nothing.