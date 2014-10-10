/*
 * TestAlign.cpp
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
#include <Align.h>
using namespace std;
using namespace Victor;

class TestAlign : public CppUnit::TestFixture {
private:
    Align *testAlign;
    
    
public:

    TestAlign() : testAlign(NULL) {
    }

    virtual ~TestAlign() {
        delete testAlign;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAlign");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlign>("Test1 - distance greater than zero.",
                &TestAlign::testAlign_A));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlign>("Test2 - distance greater than zero.",
                &TestAlign::testAlign_B));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlign>("Test3 - distance greater than zero.",
                &TestAlign::testAlign_C));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testAlign_A() {
        
      

        CPPUNIT_ASSERT(true);
    }

    void testAlign_B() {


        CPPUNIT_ASSERT(true);
    }

    void testAlign_C() {


        CPPUNIT_ASSERT(true);
    }

};
