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

#include <RapdfPotential.h>

using namespace std;
using namespace Victor;

class TestRapdfPotential : public CppUnit::TestFixture {
private:
	RapdfPotential *testRapdfPotential;
public:
	TestRapdfPotential() : testRapdfPotential(NULL) {}
	virtual ~TestRapdfPotential() {
		delete testRapdfPotential;
	}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestRapdfPotential");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestRapdfPotential>("Test1 - distance greater than zero.",
				&TestRapdfPotential::testRapdfPotential_A ));

		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
	void testRapdfPotential_A() {
		
		CPPUNIT_ASSERT( true );
	}

	
};
