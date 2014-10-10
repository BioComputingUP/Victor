/*
 * TestAlign.cpp
 *
 *  Created on: Oct 6th, 2014
 *      Author: Layla Hirsh 
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>
#include <AGPFunction.h>
#include <ScoringS2S.h>
#include <SequenceData.h>
#include <SecSequenceData.h>
#include <AlignmentBase.h>
#include <Alignment.h>
#include <NWAlign.h>
#include <Align.h>
using namespace std;
using namespace Victor;
using namespace Victor::Align2;

class TestAlign : public CppUnit::TestFixture {
private:
    Align *testAlign;
    ScoringScheme *ss;
    AlignmentData *ad;

    GapFunction *gf;
    Structure *str;
public:

    TestAlign() : testAlign(NULL) {
        string matrixFileName = "blosum62.dat";
        string matrixStrFileName = "secid.dat";
        double cSeq;
        double openGapPenalty = 12, extensionGapPenalty = 3;
        string seq1Name, seq2Name, seq1, seq2, sec1, sec2;
        string path = getenv("VICTOR_ROOT");
        string dataPath = path + "Align2/Tests/data/";
        matrixFileName = dataPath + matrixFileName;
        ifstream matrixFile(matrixFileName.c_str());
        if (!matrixFile)
            ERROR("Error opening substitution matrix file.", exception);
        string inputFileName = "test.fasta";
        if (inputFileName != "!") {
            inputFileName = dataPath + inputFileName;
            ifstream inputFile(inputFileName.c_str());
            if (!inputFile)
                ERROR("Error opening input FASTA file.", exception);
            Alignment ali;
            ali.loadFasta(inputFile);
            if (ali.size() < 1)
                ERROR("Input FASTA file must contain two sequences.", exception);
            seq1Name = ali.getTargetName();
            seq2Name = ali.getTemplateName();
            seq1 = Alignment::getPureSequence(ali.getTarget());
            seq2 = Alignment::getPureSequence(ali.getTemplate());
        }
        matrixStrFileName = dataPath + matrixStrFileName;
        ifstream matrixStrFile(matrixStrFileName.c_str());
        if (!matrixStrFile)
            ERROR("Error opening structural substitution matrix file.", exception);
        string secFileName = "t0111.sec";
        if (secFileName != "!") {
            secFileName = dataPath + secFileName;
            ifstream secFile(secFileName.c_str());
            if (!secFile)
                ERROR("Error opening secondary structure FASTA file.", exception);
            Alignment aliSec;
            aliSec.loadFasta(secFile);
            if (aliSec.size() < 1)
                ERROR("Secondary structure FASTA file must contain two sequences.", exception);
            sec1 = Alignment::getPureSequence(aliSec.getTarget());
            sec2 = Alignment::getPureSequence(aliSec.getTemplate());

        }
        SubMatrix sub(matrixFile);
        SubMatrix subStr(matrixStrFile);
        Structure *str;
        ScoringScheme *ss;
        ad = new SequenceData(2, seq1, seq2, seq1Name, seq2Name);
        str = 0;
        cSeq = 1.00;
        ss = new ScoringS2S(&sub, ad, str, cSeq);
        gf = new AGPFunction(openGapPenalty, extensionGapPenalty);
        testAlign = new NWAlign(ad, gf, ss);
    }

    virtual ~TestAlign() {
        delete testAlign;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAlign");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlign>("Test1 - comparing sequences lengths.",
                &TestAlign::testAlign_A));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlign>("Test2 - loading the ScoringScheme.",
                &TestAlign::testAlign_B));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlign>("Test3 - setting penalty values.",
                &TestAlign::testAlign_C));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testAlign_A() {



        CPPUNIT_ASSERT(testAlign->m == testAlign->n);
    }

    void testAlign_B() {


        CPPUNIT_ASSERT(testAlign->ss->ad->name1 == ad->name1);
    }

    void testAlign_C() {
        //setting penalty adding of 10 and mul of 14
        cout<<"Penalty mul: "<<testAlign->penaltyMul<<" Penalty add "<<testAlign->penaltyAdd<<"\n"; 
        testAlign->setPenalties( 14, 10); 
        cout<<"Penalty mul: "<<testAlign->penaltyMul<<" Penalty add "<<testAlign->penaltyAdd<<"\n"; 
        CPPUNIT_ASSERT((testAlign->penaltyMul== 14 )&&(testAlign->penaltyAdd== 10 ));
    }

};
