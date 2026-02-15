/*
 * TestCIFLoader.h
 *
 *  Created on: 2016/01/13
 *      Author: Francesco Menniti
 */

#include <iostream>
#include <fstream>
#include <stdlib.h> 

#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <CIFFileManager.h>
#include <CIFLoader.h>


using namespace std;
using namespace Victor::Biopool;

class TestCIFLoader : public CppUnit::TestFixture {

private:

public:

	TestCIFLoader() {;}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestCIFLoader");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestCIFLoader>("TEST TestCIFLoader.",
				&TestCIFLoader::test_A ));

		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
	void test_A() {

	        cout << "TEST TestCIFLoader" << endl;

		ifstream CIFFile("data/CIFTest.cif");
		CIFLoader cl(CIFFile);

		vector<char> chains;
		chains.push_back('A');
		chains.push_back('B');
		chains.push_back('C');
		chains.push_back('D');

		CPPUNIT_ASSERT(cl.getMaxModels() == 2);
		CPPUNIT_ASSERT(cl.getAllChains().size() == 4);
		CPPUNIT_ASSERT(cl.getAllChains()[0] == chains[0]);
		CPPUNIT_ASSERT(cl.getAllChains()[1] == chains[1]);
		CPPUNIT_ASSERT(cl.getAllChains()[2] == chains[2]);
		CPPUNIT_ASSERT(cl.getAllChains()[3] == chains[3]);

		Protein prot;
		cl.setNoVerbose();
		cl.setNoHAtoms();
		prot.load(cl);

	}
};
