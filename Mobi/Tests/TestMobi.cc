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
/*!
 *  \author    Luca Demo
 *  \date      2015
 */

#include <iostream>
#include <cppunit/TestSuite.h>
#include <cppunit/ui/text/TestRunner.h>

#include <TestVectorCollection.h>
#include <TestTM.h>
#include <TestMobiMethods.h>
#include <MadTests.h>
#include "TestMobiProtein.h"
using namespace std;


int main() {
	CppUnit::TextUi::TestRunner runner;

	cout << "Creating Test Suites:" << endl;

	runner.addTest(TestVectorCollection::suite());
	runner.addTest(TestMobiProtein::suite());
	runner.addTest(TestTM::suite());
	runner.addTest(TestMobiMethods::suite());
//	runner.addTest(MadTests::suite());
	cout<< "Running the unit tests. " << endl;

	runner.run();
	return 0;
}
