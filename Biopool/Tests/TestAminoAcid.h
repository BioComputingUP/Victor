/*
 * TestAminoAcid.cpp
 *
 *  Created on: Oct 7th, 2014
 *      Author: Layla Hirsh
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>
#include <AminoAcid.h>
#include <Atom.h>

using namespace std;
using namespace Biopool;

class TestAminoAcid : public CppUnit::TestFixture {
private:
	AminoAcid *testAminoAcid;
public:
	TestAminoAcid() : testAminoAcid(NULL) {}
	virtual ~TestAminoAcid() {
		delete testAminoAcid;
	}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAminoAcid");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestAminoAcid>("Test1 - Initialize the object.",
				&TestAminoAcid::testTestAminoAcid_A ));

	/*	suiteOfTests->addTest(new CppUnit::TestCaller<TestAminoAcid>("Test2 - zero distance.",
				&TestAminoAcid::testTestAminoAcid_B ));
                
                suiteOfTests->addTest(new CppUnit::TestCaller<TestAtom>("Test3 - rotation.",
				&TestAminoAcid::testTestAminoAcid_C ));
*/
		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
	void testTestAminoAcid_A() {
		// Calculate the distance between the two atoms
                
                CPPUNIT_ASSERT(  );
	}

	/*void testTestAtom_B() {
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
	}*/


};
