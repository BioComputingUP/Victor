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
#include <MobiMethods.h>
#include <MobiUtils.h>
#include <AtomCode.h>
#include <VectorCollection.h>

using namespace std;
using namespace Victor::Mobi;

class TestMobiMethods : public CppUnit::TestFixture {

public:

	TestMobiMethods(){}

	virtual ~TestMobiMethods(){}

	static CppUnit::Test *suite() {
	        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestMobiMethods");

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestMobiMethods>("Test1 - Test Ordering models",
	        	    &TestMobiMethods::TestOrderModel));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestMobiMethods>("Test2 - Scaled Distance",
					&TestMobiMethods::testDistances));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestMobiMethods>("Test3 - Filters",
					&TestMobiMethods::testFilters));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestMobiMethods>("Test4 - Complete Mobi test",
	                &TestMobiMethods::testMobiMobility));

	        return suiteOfTests;

	    }

	void setUp(){}

	void tearDown(){}

protected:

	void testDistances(){
		cout << ">>>\tTestMobiMethods >>> Test Scaled distance" << endl;
		string path = getenv("VICTOR_ROOT");
		string inputFile = path + "Mobi/Tests/data/1AB2_input.pdb";
		string TMDir = path + "Mobi/Tests/data/";
		ifstream inFile(inputFile.c_str());
		if (!inFile)
			ERROR("Input file not found.", exception);

		TMScoreBinder *tm = new TMScoreBin(TMDir+"TMScore", TMDir, false);

		PdbLoader pl(inFile);
		pl.setNoVerbose();
		MobiProtein* prot = new MobiProtein();
		vector<unsigned int> models = vector<unsigned int>();
		models.push_back(1);
		models.push_back(2);
		models.push_back(3);
		models.push_back(4);
		prot->load(pl,models);
		CPPUNIT_ASSERT(prot->size() == 4);

		VectorCollection<double> sDist;
		VectorCollection<double> dist;
		MobiMethods mm;

		MobiMethods::distances(prot,tm,sDist,dist,mm);
		vector<double> sd = sDist.mean();
		vector<double> sdsd = sDist.stdDev();
		CPPUNIT_ASSERT(sd[0] > 0.106 && sd[0] < 0.107);

//		cout << "isDone()" << mm.isDone() << endl;
//		vector<int> mobilityTrack = mm.mobiMobility(prot,tm);
//		cout << "isDone()" << mm.isDone() << endl;
//		VectorCollection<double> scDist = mm.getScaledDistances();
//		vector<double> scDistVals = scDist.mean();
//		vector<int> track0 = mm.getPsiMobility();
//		for (unsigned int i = 0; i < track0.size(); i++)
//			cout << track0[i];
//		cout << endl;
//		VectorCollection<double> psiAngles;
//		MobiMethods::psis(prot,psiAngles);
//		vector<double> psiDevs = psiAngles.stdDev();
//		for (unsigned int i = 0; i < psiDevs.size(); i++)
//			cout << (psiDevs[i] > 20 ? 1 : 0);
//		cout << endl;
	}



	void testMobiMobility(){
		cout << endl << ">>>\tTestMobiMethods >>> Full Mobi Mobility" << endl;
		string path = getenv("VICTOR_ROOT");
		string inputFile = path + "Mobi/Tests/data/1AB2_input.pdb";
		string TMDir = path + "Mobi/Tests/data/";
		ifstream inFile(inputFile.c_str());
		if (!inFile)
			ERROR("Input file not found.", exception);

		//loading models
		cout << "loading models..." << endl;
		PdbLoader pl(inFile);
		MobiProtein* prot = new MobiProtein();
		vector<unsigned int> models = prot->load(pl);
		TMScoreBin tm(TMDir+"TMScore", TMDir, false);
		MobiMethods mm;
		mm.verbosity(1);
		cout << "Call to mobiMobility..." << endl;
		vector<int> mobility = mm.mobiMobility(prot,&tm);

		//this is the expected mobility as for original Mobi
		string expected = "MMMMMMM.................MMMMMM.......MMMMMMM......................................................MMMMMMMMMMM";
		CPPUNIT_ASSERT(mobility.size() == expected.size());
		double match = .0;
		for (unsigned int i = 0; i < mobility.size(); i++)
			if ((mobility[i] == 1 && expected[i] == 'M') ||(mobility[i] == 0 && expected[i] == '.'))
				match += 1.0;
		match /= mobility.size();
		cout << "Obtained mobility is " << match << "% compatible (>90%) with the expected one\nMismatch on " << (int)(mobility.size()-(match*mobility.size())) << " residues over "<< mobility.size() << endl;
		CPPUNIT_ASSERT(match > .9);
		delete prot;
	}

	void testFilters(){
		int   sd[40] = {1,0,0,1,0,0,0,1,1,0,1,1,0,1,1,0,0,1,1,0,1,1,1,1,1,0,0,0,1,0,1,0,1,0,1,1,0,0,1,0};
		int  dev[40] = {1,1,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,1,0,0,1,0,0,1,1};
		int  psi[40] = {1,1,0,0,0,0,1,1,1,0,0,0,0,1,0,0,0,1,1,1,0,0,1,1,1,0,0,0,1,0,1,0,1,0,0,1,0,0,1,1};
		int  phi[40] = {0,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,1,1,0,1,1,0,0,1,0,1,1,1};
		int dssp[40] = {2,0,1,1,2,0,1,0,2,0,0,0,1,1,0,0,0,1,1,0,0,1,2,1,0,0,0,1,1,2,0,1,1,0,1,2,1,1,2,2};
		vector<int> sdv(sd,sd + sizeof sd / sizeof sd[0]);
		vector<int> devv(dev,dev + sizeof dev / sizeof dev[0]);
		vector<int> psiv(psi,psi + sizeof psi / sizeof psi[0]);
		vector<int> phiv(phi,phi + sizeof phi / sizeof phi[0]);
		vector<int> dsspv(dssp,dssp + sizeof dssp / sizeof dssp[0]);
		vector<unsigned int> a;
		MobiMethods mm;
		mm.verbosity(0);


		vector<int> w = MobiMethods::SDFilters(sdv,devv,dsspv,phiv,psiv,mm);
 	}

	void TestOrderModel(){
		cout << endl << ">>>\tTestMobiUtils >>> Ordering Models" << endl;
		string path = getenv("VICTOR_ROOT");
		string inputFile = path + "Mobi/Tests/data/1AB2_input.pdb";
		string TMDir = path + "Mobi/Tests/data/";
		ifstream inFile(inputFile.c_str());
		if (!inFile)
			ERROR("Input file not found.", exception);

		//loading models
		PdbLoader pl(inFile);
		MobiProtein* prot = new MobiProtein();
		prot->load(pl);
		TMScoreBin tm(TMDir+"TMScore", TMDir, false);

		vector<int> sorted = MobiUtils::sortModels(prot,&tm);
		CPPUNIT_ASSERT(sorted[0] == 16);
		delete prot;
	}
};


