/*
 * TestLin2dSolver.cpp
 *
 *  Created on: Feb 16, 2013
 *      Author: Arvind A de Menezes Pereira
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include "Lin2dSolver.h"

using namespace std;

class TestLin2dSolver : public CppUnit::TestFixture {
private:
	Lin2dSolver *testLin2dSolver;
public:
	TestLin2dSolver() : testLin2dSolver(NULL) {}
	virtual ~TestLin2dSolver() {
		delete testLin2dSolver;
	}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("Test2dLinearEquationSolver");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestLin2dSolver>("Test1 - Unique Solution.",
				&TestLin2dSolver::test2dLinearEqSolver_UniqueSolution ));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestLin2dSolver>("Test2 - No Solution.",
				&TestLin2dSolver::test2dLinEqSolver_NoSolution ));

		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
	void test2dLinearEqSolver_UniqueSolution() {
		// Should work since we can solve this easily by hand.
		double a = 2, b = 8, e = 16, c = 4, d = -8, f = 2; // 2x + 8y = 16; 4x -8y = 2
		Lin2dSolution expSol; expSol.x = 3; expSol.y = 1.25; expSol.isValid = true;
		Lin2dSolution l2dsol = Lin2dSolver::Solve(a,b,e,c, d, f);
		cout<<endl<<"Solution is: x="<<l2dsol.x<<", y="<<l2dsol.y<<endl;
		CPPUNIT_ASSERT( l2dsol==expSol );
	}

	void test2dLinEqSolver_NoSolution() {
		/// Uncomment to see an exception being thrown.
		/// Should throw an exception since this has an error.
		double a2 = 2, b2 = 8, e2=16, c2=4, d2=16, f2=20; // 2x + 8y = 16; 4x + 16y = 20
		Lin2dSolution expSol2; expSol2.x = 0; expSol2.y = 0; expSol2.isValid = false;
		Lin2dSolution l2dsol2 = Lin2dSolver::Solve(a2,b2,e2,c2,d2,f2);
		cout<<endl<<"Solution is: x="<<l2dsol2.x<<", y="<<l2dsol2.y<<endl;
		CPPUNIT_ASSERT( l2dsol2 == expSol2 );
	}


};
