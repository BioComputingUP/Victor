/*
 * TestAlignmentData.cpp
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

#include <AlignmentData.h>
#include <Alignment.h>
#include <SecSequenceData.h>

using namespace std;
using namespace Biopool;

class TestAlignmentData : public CppUnit::TestFixture {
private:
    AlignmentData *testAlignmentData;


public:

    TestAlignmentData() : testAlignmentData(NULL) {
         string seq1Name, seq2Name, seq1, seq2, sec1, sec2;
        string path = getenv("VICTOR_ROOT");
        string examplesPath = path + "Align2/Tests/data/";
         
        string inputFileName = examplesPath + inputFileName;
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

        string secFileName = examplesPath + secFileName;
        ifstream secFile(secFileName.c_str());
        if (!secFile)
            ERROR("Error opening secondary structure FASTA file.", exception);
        Alignment aliSec;
        aliSec.loadFasta(secFile);
        if (aliSec.size() < 1)
            ERROR("Secondary structure FASTA file must contain two sequences.", exception);
        sec1 = Alignment::getPureSequence(aliSec.getTarget());
        sec2 = Alignment::getPureSequence(aliSec.getTemplate());
        testAlignmentData = new SecSequenceData(4, seq1, seq2, sec1, sec2, seq1Name, seq2Name);
    }

    virtual ~TestAlignmentData() {
        delete testAlignmentData;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAlignmentData");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentData>("Test1 - distance greater than zero.",
                &TestAlignmentData::testAlignmentData_A));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentData>("Test2 - distance greater than zero.",
                &TestAlignmentData::testAlignmentData_B));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentData>("Test3 - distance greater than zero.",
                &TestAlignmentData::testAlignmentData_C));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testAlignmentData_A() {
       // AlignmentData *ad;
        cout<<testAlignmentData->name1;
       cout<<testAlignmentData->name2;
       // ad = new SecSequenceData(4, seq1, seq2, sec1, sec2, seq1Name, seq2Name);

        CPPUNIT_ASSERT(true);
    }

    void testAlignmentData_B() {


        CPPUNIT_ASSERT(true);
    }

    void testAlignmentData_C() {


        CPPUNIT_ASSERT(true);
    }

};
