/*
 * TestEnergy.cc
 *
 *  Created on: Oct 7, 2014
 *      Author: Manuel Giollo
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include <TestRapdfPotential.h>
#include <TestSolvationPotential.h>
using namespace std;
using namespace Victor;
using namespace Victor::Biopool;
using namespace Victor::Energy;

int main() {
	CppUnit::TextUi::TestRunner runner;

	cout << "Creating Test Suites for Energy:" << endl;
        runner.addTest(TestRapdfPotential::suite());
        runner.addTest(TestSolvationPotential::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();

	return 0;
}
