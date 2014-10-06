/*
 * TestAtom.cpp
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

#include <Atom.h>

using namespace std;
using namespace Biopool;

class TestAtom : public CppUnit::TestFixture {
private:
	Atom *testAtom;
public:
	TestAtom() : testAtom(NULL) {}
	virtual ~TestAtom() {
		delete testAtom;
	}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAtom");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestAtom>("Test1 - distance greater than zero.",
				&TestAtom::testTestAtom_A ));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestAtom>("Test2 - zero distance.",
				&TestAtom::testTestAtom_B ));
                
                suiteOfTests->addTest(new CppUnit::TestCaller<TestAtom>("Test3 - rotation.",
				&TestAtom::testTestAtom_C ));

		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
	void testTestAtom_A() {
		// Calculate the distance between the two atoms
                Atom atom0;
                atom0.setCoords(0,0,0);
                Atom atom1;
                atom1.setCoords(1,1,0);
                //this is the expected value
                double distance = sqrt(2);
                
		cout << endl << "Solution is: x=" << atom0.distance(atom1)
                        << ", y=" << distance << endl;
		CPPUNIT_ASSERT( atom0.distance(atom1) == distance );
	}

	void testTestAtom_B() {
                // Calculate the distance between the two atoms
		Atom atom0;
                atom0.setCoords(0,0,0);
                Atom atom1;
                atom1.setCoords(0,0,0);
                //this is the expected value
                double distance = sqrt(0);
                
		cout << endl << "Solution is: x=" << atom0.distance(atom1)
                        << ", y=" << distance << endl;
		CPPUNIT_ASSERT( atom0.distance(atom1) == distance );
	}
        
        void testTestAtom_C() {
                // Calculate the distance between the two atoms
		Atom atom0;
                atom0.setCoords(0,0,0);
                Atom atom1;
                atom1.setCoords(0,0,0);
                //this is the expected value
                double distance = sqrt(0);
                
		cout << endl << "Solution is: x=" << atom0.distance(atom1)
                        << ", y=" << distance << endl;
		CPPUNIT_ASSERT( atom0.distance(atom1) == distance );
	}


};
