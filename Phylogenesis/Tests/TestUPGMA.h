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
#include <UPGMA.h>
#include <DistanceMatrix.h>
#include <Cluster.h>
#include <vector>
#include "../Sources/Sequence.h"

using namespace std;
using namespace Victor;
using namespace Victor::Phylogenesis;

class TestUPGMA : public CppUnit::TestFixture {
private:
	UPGMA* sut;

public:

	virtual ~TestUPGMA() {
	}

	void setUp() {
	}

	void tearDown() {
	}

	void test_ShouldClusteringOneSpecies() {
		vector<Sequence*> species;
		species.push_back(new Sequence("s1","AAA", 0));

		sut = new UPGMA(species, NULL);
		Cluster* result = sut->performClustering();

		CPPUNIT_ASSERT(result->branches.size() == 0);
		CPPUNIT_ASSERT(result->distance == 0);
		CPPUNIT_ASSERT(result->isSingle() == true);
	}

	void test_ShouldClusteringTwoSpecies() {

		DistanceMatrix* dm = new DistanceMatrix();
		dm->appendDistance(0,100);

		vector<Sequence*> species;
		species.push_back(new Sequence("s1","AAA", 0));
		species.push_back(new Sequence("s2","CCC", 1));

		sut = new UPGMA(species, dm);
		Cluster* result = sut->performClustering();

		CPPUNIT_ASSERT(result->branches.size() == 2);
		CPPUNIT_ASSERT(result->distance == 50);
		CPPUNIT_ASSERT(result->isSingle() == false);
		CPPUNIT_ASSERT(result->branches[0]->isSingle() == true);
		CPPUNIT_ASSERT(result->branches[1]->isSingle() == true);
		CPPUNIT_ASSERT(result->lengths[0] == 50);
		CPPUNIT_ASSERT(result->lengths[1] == 50);
	}

	void test_ShouldClusteringThreeSpecies() {

		DistanceMatrix* dm = new DistanceMatrix();
		dm->appendDistance(0,100);
		dm->appendDistance(0,200);
		dm->appendDistance(1,300);
		//0 1 100
		//0 2 200
		//1 2 300

		vector<Sequence*> species;
		species.push_back(new Sequence("s1","AAA", 0));
		species.push_back(new Sequence("s2","CCC", 1));
		species.push_back(new Sequence("s3","DDD", 2));

		sut = new UPGMA(species, dm);
		Cluster* result = sut->performClustering();

		CPPUNIT_ASSERT(result->branches.size() == 2);
		CPPUNIT_ASSERT(result->distance == 125);
		CPPUNIT_ASSERT(result->lengths[1] == 125);
		CPPUNIT_ASSERT(result->lengths[0] == 75);

		Cluster* left = result->branches[1];
		Cluster* right = result->branches[0];

		CPPUNIT_ASSERT(left->isSingle());
		CPPUNIT_ASSERT(right->branches.size() == 2);
		CPPUNIT_ASSERT(right->distance == 50);
		CPPUNIT_ASSERT(right->isSingle() == false);
		CPPUNIT_ASSERT(right->branches[0]->isSingle() == true);
		CPPUNIT_ASSERT(right->branches[1]->isSingle() == true);
		CPPUNIT_ASSERT(right->lengths[0] == 50);
		CPPUNIT_ASSERT(right->lengths[1] == 50);
	}

	static CppUnit::Test *suite() {

		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestUPGMA");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestUPGMA>("Test1 - ShouldClusteringOneSpecies.",
						&TestUPGMA::test_ShouldClusteringOneSpecies));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestUPGMA>("Test2 - ShouldClusteringTwoSpecies.",
				&TestUPGMA::test_ShouldClusteringTwoSpecies));

		suiteOfTests->addTest(new CppUnit::TestCaller<TestUPGMA>("Test3 - ShouldClusteringThreeSpecies.",
						&TestUPGMA::test_ShouldClusteringThreeSpecies));

		return suiteOfTests;
	}

};
