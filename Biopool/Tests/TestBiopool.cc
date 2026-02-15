/*
 * TestBiopool.cc
 *
 *  Created on: Oct 6, 2014
 *      Author: Manuel Giollo
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>
#include <TestAtom.h>
#include <TestGroup.h>
#include <TestAminoAcid.h>
#include <TestSpacer.h>

using namespace std;


int main() {

    cout << "Main" << endl;

	
	CppUnit::TextUi::TestRunner runner;
	cout << "Creating Test Suites:" << endl;
        runner.addTest(TestAtom::suite());
        runner.addTest(TestGroup::suite());
        runner.addTest(TestAminoAcid::suite());
        runner.addTest(TestSpacer::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();
	
	return 0;
}
