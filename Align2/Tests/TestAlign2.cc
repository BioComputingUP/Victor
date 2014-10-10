/*
 * TestAlign2.cc
 *
 *  Created on: Oct 7, 2014
 *      Author: Manuel Giollo
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>
#include <TestAlignmentData.h>
#include <TestAlign.h>
#include <TestAlignment.h>
using namespace std;


int main() {
	CppUnit::TextUi::TestRunner runner;

	cout << "Creating Test Suites:" << endl;
        runner.addTest(TestAlignmentData::suite());
        runner.addTest(TestAlignment::suite());
        runner.addTest(TestAlign::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();

	return 0;
}
