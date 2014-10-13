/*
 * TestBiopool.cc
 *
 *  Created on: Oct 6, 2014
 *      Author: Layla Hirsh
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>
#include <TestVectorTransformation.h>
#include <TestLoopModel.h>
using namespace std;


int main() {
	CppUnit::TextUi::TestRunner runner;

	cout << "Creating Test Suites:" << endl;
        runner.addTest(TestLoopModel::suite());
         runner.addTest(TestVectorTransformation::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();

	return 0;
}
