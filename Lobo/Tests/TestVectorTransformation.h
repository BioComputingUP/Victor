/*
 * TestRapdfPotential.cpp
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

#include <VectorTransformation.h>

using namespace std;
using namespace Victor;

class TestVectorTransformation : public CppUnit::TestFixture {
private:
	VectorTransformation *testVectorTransformation;
public:
	TestVectorTransformation() : testVectorTransformation(NULL) {}
	virtual ~TestVectorTransformation() {
		delete testVectorTransformation;
	}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestVectorTransformation");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestVectorTransformation>("Test1 - distance greater than zero.",
				&TestVectorTransformation::testVectorTransformation_A ));

		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
	void testVectorTransformation_A() {
		
		CPPUNIT_ASSERT( true );
	}

	
};
