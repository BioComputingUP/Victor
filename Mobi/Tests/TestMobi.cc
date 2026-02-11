/*
 @file    TestMobi.cc
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 1.0
 */

/* This file is part of Victor.
 Victor is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 Victor is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with Victor. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include <TestProteinModels.h>
#include <TestStandardDeviation.h>
#include <TestSecondaryStructure.h>
#include <TestSaver.h>


using namespace Victor::Mobi;
using namespace Victor::Biopool;
using namespace std;

int main() {
	CppUnit::TextUi::TestRunner runner;

	cout << "Creating Test Suites:" << endl;
	runner.addTest(TestProteinModels::suite());
	runner.addTest(TestStandardDeviation::suite());
	runner.addTest(TestSecondaryStructure::suite());
	runner.addTest(TestSaver::suite());

	cout << "Running the unit tests. " << endl;
	runner.run();
	return 0;
}
