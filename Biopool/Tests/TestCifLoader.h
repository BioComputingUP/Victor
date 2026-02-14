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

#ifndef TESTCIFLOADER_H
#define	TESTCIFLOADER_H

#include <iostream>
#include <fstream>

#include <cppunit/TestFixture.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>

#include <CifLoader.h>

using namespace std;
using namespace Victor::Biopool;
using namespace CppUnit;

class TestCifLoader : public TestFixture {
public:

    TestCifLoader() {
    }

    virtual ~TestCifLoader() {
    }

    static Test* suite() {
	TestSuite* suiteOfTests = new TestSuite("TestCifLoader");

	suiteOfTests->addTest(new TestCaller<TestCifLoader>("Test 1 - Get max numebers of models",
		&TestCifLoader::testGetMaxModels));

	suiteOfTests->addTest(new TestCaller<TestCifLoader>("Test 2 - Get all chain ids",
		&TestCifLoader::testGetAllChains));

	return suiteOfTests;
    }

    void setUp() {
	// inizialize CifLoader
	string path = getenv("VICTOR_ROOT");
	string input = path + "Biopool/Tests/data/modelTest.cif";
	inFile = new ifstream(input.c_str());
	cl = new CifLoader(*inFile);
	
	// initialize test parameters
	maxModel = 5;
	chainIds.push_back('A');
	chainIds.push_back('B');
	chainIds.push_back('C');
    }

    void tearDown() {
	delete cl;
	delete inFile;
    }

private:

    void testGetMaxModels() {
	int max = cl->getMaxModels();
	CPPUNIT_ASSERT_EQUAL(maxModel, max);
    }

    void testGetAllChains() {
	vector<char> testChain = cl->getAllChains();
	CPPUNIT_ASSERT(chainIds[0] == testChain[0] &&
		chainIds[1] == testChain[1] &&
		chainIds[2] == testChain[2]);
    }

    int maxModel;
    vector<char> chainIds;
    
    ifstream* inFile;
    CifLoader* cl;
};
#endif	/* TESTCIFLOADER_H */

