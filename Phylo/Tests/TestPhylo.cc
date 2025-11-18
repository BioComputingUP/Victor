/*  This file is part of Victor.
    
    Victor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    Victor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with Victor.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>
#include <TestTree.h>
#include <TestDistanceMatrix.h>
#include <TestClustering.h>
#include <TestMSA.h>

using namespace::std;
using namespace::Victor;

int main(int argc, char *argv[]) {
    CppUnit::TextUi::TestRunner runner;

    cout << "Creating Test Suite:" << endl;
    runner.addTest(TestTree::suite());
    runner.addTest(TestDistanceMatrix::suite());
    runner.addTest(TestClustering::suite());
    runner.addTest(TestMSA::suite());
    cout << "Running the unit tests:" << endl;
    runner.run();

    (void) argc;
    (void) argv;

    return 0;
}