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
#include <ClusterDistanceMatrix.h>
#include <Cluster.h>

//To avoid "Debug.h" definition
#undef exception
#include <stdexcept>

using namespace std;
using namespace Victor;
using namespace Victor::Phylogenesis;

class TestClusterDistanceMatrix :  public CppUnit::TestFixture {
private:
	ClusterDistanceMatrix* sut;
public:

    virtual ~TestClusterDistanceMatrix() {
    }

    void setUp() {
    }

    void tearDown() {
    }

    void test_shouldInsertAndRetrieve() {
    	Cluster* c1 = new Cluster();
    	Cluster* c2 = new Cluster();

    	sut = new ClusterDistanceMatrix();
    	sut->insert(c1, c2, 0.10);

        CPPUNIT_ASSERT(sut->at(c1,c2) == sut->at(c2,c1));
        CPPUNIT_ASSERT(sut->at(c1,c2) == 0.10);
    }

    void test_shouldEraseCorrectly() {
    	Cluster* c1 = new Cluster();
		Cluster* c2 = new Cluster();
		Cluster* c3 = new Cluster();
		Cluster* c4 = new Cluster();

		sut = new ClusterDistanceMatrix();
		sut->insert(c1, c2, 0.10);
		sut->insert(c1, c3, 0.20);
		sut->insert(c2, c3, 0.30);
		sut->insert(c2, c4, 0.40);
		sut->insert(c3, c4, 0.50);

		sut->erase(c1, c4);

		CPPUNIT_ASSERT(sut->at(c2,c3) == 0.30);
		CPPUNIT_ASSERT_THROW(sut->at(c1,c2), std::out_of_range);
		CPPUNIT_ASSERT_THROW(sut->at(c1,c3), std::out_of_range);
		CPPUNIT_ASSERT_THROW(sut->at(c2,c4), std::out_of_range);
		CPPUNIT_ASSERT_THROW(sut->at(c3,c4), std::out_of_range);
    }

    static CppUnit::Test *suite() {

		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestClusterDistanceMatrix");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestClusterDistanceMatrix>("Test1 - ShouldInsertAndRetrieve.",
				&TestClusterDistanceMatrix::test_shouldInsertAndRetrieve));
		suiteOfTests->addTest(new CppUnit::TestCaller<TestClusterDistanceMatrix>("Test2 - ShouldEraseCorrectly.",
				&TestClusterDistanceMatrix::test_shouldEraseCorrectly));

		return suiteOfTests;
	}


};
