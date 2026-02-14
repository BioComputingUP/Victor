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
#include <PdbSaver.h>
#include <TMScoreBin.h>
#include <AminoAcid.h>

using namespace std;
using namespace Victor::Mobi;

class TestTM : public CppUnit::TestFixture {

public:

	TestTM(){}

	virtual ~TestTM(){}

	static CppUnit::Test *suite() {
	        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestTM");

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestTM>("Test1 - TM Score",
	                &TestTM::testLoadAndGet));

	        return suiteOfTests;
	    }

	void setUp(){}

	void tearDown(){}

protected:

	void testLoadAndGet(){
		cout << endl << ">>>\tTestTM >>> test TM binary call" << endl;
		string path = getenv("VICTOR_ROOT");
		string inputFile = path + "Mobi/Tests/data/1AB2_input.pdb";
		string TMDir = path + "Mobi/Tests/data/";
		ifstream inFile(inputFile.c_str());
		if (!inFile)
			ERROR("Input file not found.", exception);
		//models to load
		vector<unsigned int> models = vector<unsigned int>();
	    models.push_back(3);
	    models.push_back(4);

	    //loading models
		PdbLoader pl(inFile);
		pl.setNoVerbose();
		MobiProtein prot;
	    prot.load(pl,models);

	    //calling TMScore bin
	    TMScoreBin tmsb(TMDir+"TMScore", TMDir, false);
	    MobiProtein* superImposed = tmsb.TMScore(prot,0,1);


	    //save superimposed model to file
	    std::ofstream fout;
		PdbSaver ps(fout);
		fout.open((TMDir + "superimposed_TM-Test.pdb").c_str());
		ps.saveSpacer(*(superImposed->getModel(0)));
		ps.endFile();
		fout.close();

		CPPUNIT_ASSERT(superImposed->getModel(0)->sizeAmino() == 109);
		CPPUNIT_ASSERT(superImposed->getModel(0)->getAmino(0).getType() == "GLY");
	}
};


