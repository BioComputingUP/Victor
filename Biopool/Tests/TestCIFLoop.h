/*
 * TestCIFFileManager.h
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
#include <CIFLoop.h>


using namespace std;
using namespace Victor::Biopool;

class TestCIFLoop : public CppUnit::TestFixture {

private:

public:

	TestCIFLoop() {;}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestCIFLoop");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestCIFLoop>("TEST TestCIFLoop.",
				&TestCIFLoop::test_A ));

		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
	void test_A() {

	        cout << "TEST TestCIFLoop" << endl;
        
		CIFLoop loop("atom_site");
		loop.AddColumn("_atom_site.group_PDB ");
		loop.AddColumn("_atom_site.id ");
		loop.AddColumn("_atom_site.type_symbol ");
		loop.AddColumn("_atom_site.label_atom_id ");

		loop.AddRow("ATOM   1    N N   ");
		loop.AddRow("ATOM   2    C CA  ");
		loop.AddRow("ATOM   3    C C   ");

		CPPUNIT_ASSERT(CIFLoop::DecodeCategory("_atom_site.group_PDB ") == "atom_site") ;        
		CPPUNIT_ASSERT(loop.ColumnsCount() == 4);
		CPPUNIT_ASSERT(loop.RowsCount() == 3);
		CPPUNIT_ASSERT(loop.GetCellValue("id", 2) == "3 ");

	}
};
