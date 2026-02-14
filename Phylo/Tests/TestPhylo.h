/*
 * TestPhylo.h
 *
 *  Created on: 02 07 2015
 *      Author: Matteo Del Pioluogo
 */


#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>
#include <PhyloSupport.h>
using namespace std;
using namespace Victor;


class TestPhylo : public CppUnit::TestFixture {
private:
private:

public:

    TestPhylo(){
    }


    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAlign");
        suiteOfTests->addTest(new CppUnit::TestCaller<TestPhylo>("Test1 - Calculate distance of 2 seq.",
                &TestPhylo::testPhylo_A));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestPhylo>("Test2 - Insert Gap in a seq.",
               &TestPhylo::testPhylo_B));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestPhylo>("Test3 - Calculate Multi Aling with weigth.",
                &TestPhylo::testPhylo_C));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestPhylo>("Test4 - Calculate Multi Aling with weigth and deletetion/insertion.",
                &TestPhylo::testPhylo_D));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testPhylo_A() {
    	string seq1="AAABBCCC";
    	string seq2="AAABBCCC";
    	double distance=Phylo::PhyloSupport::distanceCalcTwoSeq(seq1,seq2);
    	CPPUNIT_ASSERT( distance == 0 );
    	seq2="DDDDDDDD";
    	distance=Phylo::PhyloSupport::distanceCalcTwoSeq(seq1,seq2);
    	CPPUNIT_ASSERT( distance == 1 );
    }

    void testPhylo_B() {

    	string seq="AAABBCCC";
    	int pos=3;
    	seq=Phylo::PhyloSupport::insertGapPosition(seq,pos);
        CPPUNIT_ASSERT(seq.substr(3,1) == "-");
    }

    void testPhylo_C() {
    	string a="AAAAAAAA";
    	string b="AAAAAAAA";
    	string c="AAAAAAAA";
    	string o="AAAAAAAA";
    	string t="AAAAAAAA";
    	vector <string> num(2);
    	vector <string> letter(3);
    	num[0]=o;
    	num[1]=t;
    	letter[0]=a;
    	letter[1]=b;
    	letter[2]=c;

    	vector <double> w1(2);
    	vector <double> w2(3);
    	w1[0]=1;
    	w1[1]=1;
    	w2[2]=1;
    	w2[0]=1;
    	w2[1]=1;
    	vector <string> result=Phylo::PhyloSupport::AlingMultiSvsMultiS(num,letter,w1,w2,false,1);
    	for(unsigned int i=0;i<result.size();i++) {
    		for(unsigned int j=0;j<result[i].size();j++) {
    			CPPUNIT_ASSERT(result[i][j] == 'A');
    		}
    	}
    }

    void testPhylo_D() {
		string a="AACCAA";
		string b="AACCAA";
		string c="AACCAA";
		string o="CC";
		string t="CC";
		vector <string> num(2);
		vector <string> letter(3);
		num[0]=o;
		num[1]=t;
		letter[0]=a;
		letter[1]=b;
		letter[2]=c;

		vector <double> w1(3);
		vector <double> w2(2);
		w1[0]=1;
		w1[1]=1;
		w1[2]=1;
		w2[0]=1;
		w2[1]=1;
		vector <string> result=Phylo::PhyloSupport::AlingMultiSvsMultiS(num,letter,w1,w2,false,2);
		CPPUNIT_ASSERT(result[0][0] == '-');
		CPPUNIT_ASSERT(result[0][2] == 'C');

    }

};
