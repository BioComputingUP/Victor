/*
 * TestCIFFileManager.h
 *
 *  Created on: Oct 6th, 2014
 *      Author: Manuel Giollo
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

class TestCIFFileManager : public CppUnit::TestFixture {

private:
	CIFFileManager *fm;

public:

	TestCIFFileManager() {

		ifstream CIFFile("data/CIFTest.cif");
		fm = new CIFFileManager(CIFFile);
		fm->Init();
	}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestCIFFileManager");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestCIFFileManager>("TEST CIFFileManager.",
				&TestCIFFileManager::testCFM_A ));

		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
	void testCFM_A() {

	        cout << "TEST CIFFileManager" << endl;

		CPPUNIT_ASSERT(fm->GetInLineValue("entry.id") == "25C8");
		CPPUNIT_ASSERT(fm->AtomLoop.RowsCount() == 72);
		CPPUNIT_ASSERT(fm->AtomLoop.ColumnsCount() == 26);
		CPPUNIT_ASSERT(fm->AtomLoop.GetCategoryGroup() == "atom_site");
		CPPUNIT_ASSERT(fm->AtomLoop.GetCellValue("Cartn_x",0) == "20.906 ");
	}
};
