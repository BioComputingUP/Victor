/*
 * TestVectorTransformation.cpp
 *
 *  Created on: Oct 6th, 2014
 *      Author: Layla Hirsh 
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <VectorTransformation.h>
#include <PdbLoader.h>
#include <Spacer.h>

using namespace std;
using namespace Victor::Lobo;

class TestVectorTransformation : public CppUnit::TestFixture {
private:
    VectorTransformation *testVectorTransformation;
public:

    TestVectorTransformation() : testVectorTransformation(NULL) {
    }

    virtual ~TestVectorTransformation() {
        delete testVectorTransformation;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestVectorTransformation");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestVectorTransformation>("Test1 - Vector initialized with the identity matrix.",
                &TestVectorTransformation::testTestVectorTransformation_A));

  


        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }
    vg_ieee64 degreesToRadian(const vg_ieee64& deg) {
        return (deg / 180.0) * M_PI;
    }
protected:

    void testTestVectorTransformation_A() {
        // Vector initialized with the identity matrix
        VectorTransformation testVector;
        
        CPPUNIT_ASSERT((testVector.getRot(0).x.x==1 )&&
        (testVector.getRot(0).x.y==0 )&&
        (testVector.getRot(0).x.z==0 )&&
        (testVector.getRot(0).y.x==0 )&&
        (testVector.getRot(0).y.y ==1)&&
        (testVector.getRot(0).y.z ==0 )&&
        (testVector.getRot(0).z.x ==0 )&&
        (testVector.getRot(0).z.y==0 )&&
        (testVector.getRot(0).z.z ==1)) ;
    }

   
};
