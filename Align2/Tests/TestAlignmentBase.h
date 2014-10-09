/*
 * TestAlignmentBase.cpp
 *
 *  Created on: Oct 6th, 2014
 *      Author: Manuel Giollo
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <AlignmentBase.h>
#include <Alignment.h>

using namespace std;
using namespace Biopool;

class TestAlignmentBase : public CppUnit::TestFixture {
private:
    AlignmentBase *testAlignmentBase;
    
    
public:

    TestAlignmentBase() : testAlignmentBase(NULL) {
    }

    virtual ~TestAlignmentBase() {
        delete testAlignmentBase;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAlignmentBase");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentBase>("Test1 - distance greater than zero.",
                &TestAlignmentBase::testAlignmentBase_A));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentBase>("Test2 - distance greater than zero.",
                &TestAlignmentBase::testAlignmentBase_B));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentBase>("Test3 - distance greater than zero.",
                &TestAlignmentBase::testAlignmentBase_C));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testAlignmentBase_A() {
        
        

        CPPUNIT_ASSERT(true);
    }

    void testAlignmentBase_B() {


        CPPUNIT_ASSERT(true);
    }

    void testAlignmentBase_C() {


        CPPUNIT_ASSERT(true);
    }

};
