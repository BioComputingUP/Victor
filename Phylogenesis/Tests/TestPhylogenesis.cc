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
#include <TestDistanceMatrix.h>
#include <TestClusterDistanceMatrix.h>
#include <TestUPGMA.h>
#include <TestNJ.h>
#include <TestNewickPhylogeneticTreeExporter.h>
#include <TestClustalw.h>

using namespace std;
using namespace Victor;

int main() {
	CppUnit::TextUi::TestRunner runner;

	cout << "Creating Test Suites:" << endl;
        runner.addTest(TestDistanceMatrix::suite());
        runner.addTest(TestClusterDistanceMatrix::suite());
        runner.addTest(TestUPGMA::suite());
        runner.addTest(TestNJ::suite());
        runner.addTest(TestNewickPhylogeneticTreeExporter::suite());
        runner.addTest(TestClustalw::suite());

	cout<< "Running the unit tests."<<endl;
	runner.run();

	return 0;
}
