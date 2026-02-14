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

#ifndef TESTCIFSTRUCTURE_H
#define	TESTCIFSTRUCTURE_H

#include <iostream>
#include <fstream>

#include <cppunit/TestFixture.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>

#include <IoTools.h>
#include <String2Number.h>
#include <CifStructure.h>

using namespace std;
using namespace Victor::Biopool;
using namespace CppUnit;

class TestCifStructure : public TestFixture {
public:

    TestCifStructure() {
    }

    virtual ~TestCifStructure() {
    }

    static Test* suite() {
	TestSuite* suiteOfTests = new TestSuite("TestCifLoader");

	suiteOfTests->addTest(new TestCaller<TestCifStructure>("Get group column number",
		&TestCifStructure::testGetGroupColumnNumber));

	suiteOfTests->addTest(new TestCaller<TestCifStructure>("Get group field",
		&TestCifStructure::testGetGroupField));
	
	suiteOfTests->addTest(new TestCaller<TestCifStructure>("Get inline field",
		&TestCifStructure::testGetInlineField));

	return suiteOfTests;
    }

    void setUp() {
	// inizialize CifLoader
	string path = getenv("VICTOR_ROOT");
	string input = path + "Biopool/Tests/data/modelTest.cif";
	inFile = new ifstream(input.c_str());
	cs = new CifStructure(*inFile);
	
	// initialize test parameters
	idColumn = 1;
	atomId = 1;
	headerLine = "_entry.id   25C8 ";
	headerId = "25C8";
    }

    void tearDown() {
	delete cs;
	delete inFile;
    }

private:

    void testGetGroupField() {
	string line = readLine(cs->getInput());
	cs->parseGroup("atom", line);
	int field = stoiDEF(cs->getGroupField("atom", line, 
		cs->getGroupColumnNumber("atom", "atom id")));
	
	CPPUNIT_ASSERT_EQUAL(atomId, field);
    }

    void testGetGroupColumnNumber() {
	string line = readLine(cs->getInput());
	cs->parseGroup("atom", line);
	int col = cs->getGroupColumnNumber("atom", "atom id");
	
	CPPUNIT_ASSERT_EQUAL(idColumn, col);
    }
    
    void testGetInlineField() {
	string field = cs->getInlineField(headerLine);
	
	CPPUNIT_ASSERT_EQUAL(headerId, field);
    }
    
    int idColumn;
    int atomId;
    string headerLine;
    string headerId;
    
    ifstream* inFile;
    CifStructure* cs;
    
};

#endif	/* TESTCIFSTRUCTURE_H */

