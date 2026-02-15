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
#include <Sequence.h>
#include <Cluster.h>
#include <SingleCluster.h>
#include <NewickPhylogeneticTreeExporter.h>
#include <iostream>

using namespace std;
using namespace Victor;
using namespace Victor::Phylogenesis;

class TestNewickPhylogeneticTreeExporter : public CppUnit::TestFixture {
private:
	NewickPhylogeneticTreeExporter* sut;

public:

	TestNewickPhylogeneticTreeExporter(){
		sut = new NewickPhylogeneticTreeExporter();
	}

	virtual ~TestNewickPhylogeneticTreeExporter() {
	}

	void test_ShouldPrintNewickTree1(){

		Sequence* s = new Sequence("A","xxx",0);
		Cluster* tree = new SingleCluster(s);

		string expectation = "A:0";
		std::ostringstream stream;
		sut->save(tree, stream);

		CPPUNIT_ASSERT(expectation == stream.str());
	}

	void test_ShouldPrintNewickTree2(){

		Sequence* s1 = new Sequence("A","xxx",0);
		Sequence* s2 = new Sequence("B","xxx",1);
		Cluster* leaf1 = new SingleCluster(s1);
		Cluster* leaf2 = new SingleCluster(s2);

		Cluster* tree = new Cluster();
		tree->branches.push_back(leaf1);
		tree->branches.push_back(leaf2);
		tree->lengths.push_back(10);
		tree->lengths.push_back(20);

		string expectation = "(A:10,B:20);";
		std::ostringstream stream;
		sut->save(tree, stream);

		CPPUNIT_ASSERT(expectation == stream.str());
	}

	void test_ShouldPrintNewickTree3(){
		Sequence* s1 = new Sequence("A","xxx",0);
		Sequence* s2 = new Sequence("B","xxx",1);
		Cluster* leaf1 = new SingleCluster(s1);
		Cluster* leaf2 = new SingleCluster(s2);
		Cluster* lleaf1 = new Cluster();
		lleaf1->branches.push_back(leaf1);
		lleaf1->branches.push_back(leaf2);
		lleaf1->lengths.push_back(10);
		lleaf1->lengths.push_back(20);
		Sequence* s3 = new Sequence("C","xxx",2);
		Cluster* lleaf2 = new SingleCluster(s3);

		Cluster* tree = new Cluster();
		tree->branches.push_back(lleaf1);
		tree->branches.push_back(lleaf2);
		tree->lengths.push_back(30);
		tree->lengths.push_back(40);

		string expectation = "((A:10,B:20):30,C:40);";
		std::ostringstream stream;
		sut->save(tree, stream);

		CPPUNIT_ASSERT(expectation == stream.str());
	}

	static CppUnit::Test *suite() {

		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestNewickPhylogeneticTreeExporter");

			suiteOfTests->addTest(new CppUnit::TestCaller<TestNewickPhylogeneticTreeExporter>("Test1 - ShouldPrintNewickTree1.",
					&TestNewickPhylogeneticTreeExporter::test_ShouldPrintNewickTree1));
			suiteOfTests->addTest(new CppUnit::TestCaller<TestNewickPhylogeneticTreeExporter>("Test2 - ShouldPrintNewickTree2.",
					&TestNewickPhylogeneticTreeExporter::test_ShouldPrintNewickTree2));
			suiteOfTests->addTest(new CppUnit::TestCaller<TestNewickPhylogeneticTreeExporter>("Test3 - ShouldPrintNewickTree3.",
					&TestNewickPhylogeneticTreeExporter::test_ShouldPrintNewickTree3));

		return suiteOfTests;
	}

	/// Setup method

	void setUp() {
	}

	void tearDown() {
	}

};
