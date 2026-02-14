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
/*!
 *  \author    Luca Demo
 *  \date      2015
 */
#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <Spacer.h>
#include <PdbLoader.h>
#include <MobiProtein.h>

using namespace std;
using namespace Victor::Mobi;

class TestMobiProtein : public CppUnit::TestFixture {

public:

	TestMobiProtein(){}

	virtual ~TestMobiProtein(){}

	static CppUnit::Test *suite() {
	        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestMobiProtein");

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestMobiProtein>("Test1 - Load and Get models",
	                &TestMobiProtein::testLoadAndGet));

	        return suiteOfTests;
	    }

	void setUp(){}

	void tearDown(){}

protected:

	void testLoadAndGet(){
		cout << endl << ">>>\tTestMobiProtein >>> test Load and Get models" << endl;
		string path = getenv("VICTOR_ROOT");
		string inputFile = path + "Mobi/Tests/data/1AB2_input.pdb";
		ifstream inFile(inputFile.c_str());
		if (!inFile)
			ERROR("Input file not found.", exception);

		PdbLoader pl(inFile);
		//pl.setNoVerbose();
	    MobiProtein prot = MobiProtein();
	    vector<unsigned int> models = vector<unsigned int>();
	    models.push_back(1);
	    models.push_back(4);
	    models.push_back(8);
	    prot.load(pl,models);
	    CPPUNIT_ASSERT(prot.size() == 3);
	    Spacer* m1 = prot.getModel(0);
	    Spacer* m2 = prot.getModel(1);
	    Spacer* m3 = prot.getModel(2);
	    CPPUNIT_ASSERT(m1->getAmino(0).getType() == "GLY");
	    CPPUNIT_ASSERT(m2->getAmino(1).getType() == "SER");
	    CPPUNIT_ASSERT(m3->getAmino(2).getType() == "GLY");
	    Spacer* m4 = prot.getModelByPdbID(8);
	    CPPUNIT_ASSERT(m4 == m3);
	}
};
