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

#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <ClustalW.h>
#include <Sequence.h>
#include <MultipleAlignment.h>
#include <UPGMA.h>

#include <vector>
#include <iostream>

using namespace std;
using namespace Victor;
using namespace Victor::Phylogenesis;

class TestClustalw :  public CppUnit::TestFixture {
private:
	DistanceMatrix* dm;
	ClustalW* sut;
public:

	TestClustalw(){
		string path = getenv("VICTOR_ROOT");
		string distanceMatrix = path + "Phylogenesis/Tests/data/distanceMatrix.dat";
		dm = new DistanceMatrix(distanceMatrix);
	}

    virtual ~TestClustalw(){
    }

    void setUp(){
    }

    void tearDown() {
    }

    void test_msa(){

    	Sequence* s1 = new Sequence("A","PEEKSAVTALWGKVNVDEYGG", 0);
    	Sequence* s2 = new Sequence("B","GEEKAAVLALWDKVNEEEYGG", 1);
    	Sequence* s3 = new Sequence("C","PADKTNVKAAWGKVGAHAGEYGA", 2);
    	vector<Sequence*> sequences;
    	sequences.push_back(s1);
    	sequences.push_back(s2);
    	sequences.push_back(s3);

    	UPGMA* builder = new UPGMA(sequences, dm);
    	ClustalW* sut = new ClustalW(dm, builder, sequences, 0.4, 10, 0.2, false, false);

    	MultipleAlignment* result = sut->performMultipleSequenceAlignment();

    	CPPUNIT_ASSERT(result->sequences[0]->sequence == "PEEKSAVTALWGKV--NVDEYGG");
    	CPPUNIT_ASSERT(result->sequences[1]->sequence == "GEEKAAVLALWDKV--NEEEYGG");
    	CPPUNIT_ASSERT(result->sequences[2]->sequence == "PADKTNVKAAWGKVGAHAGEYGA");
    }

	static CppUnit::Test *suite() {

		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestClustalw");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestClustalw>("Test1 - MSA.",
				&TestClustalw::test_msa));

		return suiteOfTests;
	}

};
