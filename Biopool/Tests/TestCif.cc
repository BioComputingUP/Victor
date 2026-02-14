/*
 * File:   TestCif.cc
 * Author: marco
 *
 * Created on 9-giu-2015, 10.45.40
 */

#include <iostream>
#include <cppunit/ui/text/TestRunner.h>

#include <TestCifLoader.h>
#include <TestCifStructure.h>

using std::cout;
using std::endl;

int main() {
    
    CppUnit::TextUi::TestRunner runner;

    cout << "Creating Test Suites:" << endl;
    
    runner.addTest(TestCifLoader::suite());
    runner.addTest(TestCifStructure::suite());
    
    cout << "Running the unit tests." << endl;
    
    runner.run();

    return 0;
    
}
