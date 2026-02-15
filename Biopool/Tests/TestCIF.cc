/*
 * TestCIF.cc
 *
 *  Created on: 2016/01/13
 *      Author: Francesco Menniti
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include <TestCIFLoader.h>
#include <TestCIFFileManager.h>
#include <TestCIFLoop.h>

using namespace std;


int main() {

    cout << "Main" << endl;
    
	CppUnit::TextUi::TestRunner runner;
	cout << "Creating Test Suites:" << endl;
        runner.addTest(TestCIFLoop::suite());
        runner.addTest(TestCIFFileManager::suite());
        runner.addTest(TestCIFLoader::suite());
	cout<< "Running the unit tests."<<endl;
	runner.run();

	return 0;
}
