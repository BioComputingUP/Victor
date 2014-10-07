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

using namespace std;
using namespace Biopool;

class TestAlignmentBase : public CppUnit::TestFixture {
private:
	AlignmentBase *testAlignmentBase;
public:
	TestAlignmentBase() : testAlignmentBase(NULL) {}
	virtual ~TestAlignmentBase() {
		delete testAlignmentBase;
	}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAlignmentBase");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentBase>("Test1 - distance greater than zero.",
				&TestAlignmentBase::testAlignmentBase_A ));

		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
	void testAlignmentBase_A() {
		
		CPPUNIT_ASSERT( true );
	}

	
};
