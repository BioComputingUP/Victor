/*
 * TestPhylogen.cc
 *
 *  Created on: 02 07 2015
 *      Author: Matteo Del Pioluogo
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include <TestPhylo.h>

using namespace std;

int main(){
	CppUnit::TextUi::TestRunner runner;
	cout << "Creating Test Suites: OK" << endl;
        runner.addTest(TestPhylo::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();

	return 0;
}
























































