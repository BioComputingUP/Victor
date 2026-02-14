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
#include <MobiMethods.h>
#include "../Sources/VectorCollection.h"

using namespace std;
using namespace Victor::Mobi;

class TestVectorCollection : public CppUnit::TestFixture {

public:

	TestVectorCollection(){}

	virtual ~TestVectorCollection(){}

	static CppUnit::Test *suite() {
	        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestVectorCollection");

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestVectorCollection>("Test1 - Populate the collection",
	                &TestVectorCollection::testPopulation));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestVectorCollection>("Test2 - Calculate means",
	                &TestVectorCollection::testMeanSD));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestVectorCollection>("Test3 - Calculate standard deviation",
	                &TestVectorCollection::testSD));

	        suiteOfTests->addTest(new CppUnit::TestCaller<TestVectorCollection>("Test4 - Calculate mobility index",
	        	                &TestVectorCollection::testMobilityIndex));

	        return suiteOfTests;
	    }

	void setUp(){}

	void tearDown(){}

protected:

	void testPopulation(){
		VectorCollection<double> sdr = VectorCollection<double>();
		cout << endl << ">>>\tTestVectorCollection >>> test Population..." << endl;
		fillWith5Results10(sdr);
		CPPUNIT_ASSERT(sdr.size() == 5);
		CPPUNIT_ASSERT(sdr.vectorsSize() == 10);
		std::map<int,std::vector<double> >::const_iterator it = sdr.begin();
		CPPUNIT_ASSERT(it->second[5] == 5);
		it++;
		CPPUNIT_ASSERT(it->second[5] == 6);
		it++;
		CPPUNIT_ASSERT(it->second[5] == 7);

		vector<int> ids = sdr.getModels();
		CPPUNIT_ASSERT(ids.size() == 6);

		CPPUNIT_ASSERT(sdr.getValuesByModel(4).size() == 2);
		CPPUNIT_ASSERT(sdr.getValuesByModel(1).size() == 1);
	}

	void testMeanSD(){
		cout << endl << ">>>\tTestVectorCollection >>> test Mean..." << endl;
		VectorCollection<double> sdr = VectorCollection<double>();
		fillWith5Results10(sdr);
		vector<double> mean = sdr.mean();
		for (unsigned int i = 0; i < 10; i++)
			CPPUNIT_ASSERT(mean[i] == i+2);
	}

	void testSD(){
		cout << endl << ">>>\tTestVectorCollection >>> test deviation calculation..." << endl;
		VectorCollection<double> sdr = VectorCollection<double>();
		fillWith5Results10(sdr);
		vector<double> sd = sdr.stdDev();
		for (unsigned int i = 0; i < 10; i++)
			CPPUNIT_ASSERT(sd[i] = sqrt(2));
	}

	void testMobilityIndex(){
		cout << endl << ">>>\tTestVectorCollection >>> test mobility index..." << endl;
		VectorCollection<double> sdr = VectorCollection<double>();
		fillWith5Results10(sdr);
		vector<double> sd = MobiMethods::mobilityIndex(sdr);
		CPPUNIT_ASSERT(sd[0] > 2.44 && sd[0] < 2.45);
	}

private:

	void fillWith5Results10(VectorCollection<double>& sdr){
		unsigned int modelLen = 10;
		vector<double> *v1 = new vector<double>(modelLen);
		vector<double> *v2 = new vector<double>(modelLen);
		vector<double> *v3 = new vector<double>(modelLen);
		vector<double> *v4 = new vector<double>(modelLen);
		vector<double> *v5 = new vector<double>(modelLen);
		for (unsigned int i = 0; i < modelLen; i++){
			(*v1)[i] = i;
			(*v2)[i] = i+1;
			(*v3)[i] = i+2;
			(*v4)[i] = i+3;
			(*v5)[i] = i+4;
		}
		sdr.addValue(1,2,*v1);
		sdr.addValue(2,3,*v2);
		sdr.addValue(3,4,*v3);
		sdr.addValue(4,4,*v4);
		sdr.addValue(5,7,*v5);
	}
};

