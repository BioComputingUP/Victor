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
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>
#include <DistanceMatrix.h>

using namespace std;
using namespace Victor;
using namespace Victor::Phylogenesis;

class TestDistanceMatrix : public CppUnit::TestFixture {

	private:
		string dataPath;
		DistanceMatrix* sut;

	public:

		TestDistanceMatrix(){
			string path = getenv("VICTOR_ROOT");
			dataPath = path + "Phylogenesis/Tests/data/";
			sut = NULL;
		}

	    virtual ~TestDistanceMatrix() {
	    }

	    void setUp() {
	    }

	    void tearDown() {
	    }

	    void test_shouldLoadDistanceMatrixFromFile() {
	    	sut = new DistanceMatrix(dataPath + "distanceMatrix.dat");

			CPPUNIT_ASSERT(sut->at(0,1) == sut->at(1,0));
			CPPUNIT_ASSERT(sut->at(0,1) == 0.10);
			CPPUNIT_ASSERT(sut->at(0,2) == sut->at(2,0));
			CPPUNIT_ASSERT(sut->at(0,2) == 0.20);
			CPPUNIT_ASSERT(sut->at(1,2) == sut->at(2,1));
			CPPUNIT_ASSERT(sut->at(1,2) == 0.30);
	    }

		static CppUnit::Test *suite() {

			CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestDistanceMatrix");

			suiteOfTests->addTest(new CppUnit::TestCaller<TestDistanceMatrix>("Test1 - ShouldLoadDistanceMatrixFromFile.",
					&TestDistanceMatrix::test_shouldLoadDistanceMatrixFromFile));

			return suiteOfTests;
		}
};
