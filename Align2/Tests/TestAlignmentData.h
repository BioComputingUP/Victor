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
using namespace Victor::Align2;

class TestAlignmentData : public CppUnit::TestFixture {
private:
    AlignmentData *testAlignmentData;


public:

    TestAlignmentData() : testAlignmentData(NULL) {
        string seq1Name, seq2Name, seq1, seq2, sec1, sec2;
        string path = getenv("VICTOR_ROOT");
        string examplesPath = path + "Align2/Tests/data/";
        string inputFileName = "test.fasta";
        inputFileName = examplesPath + inputFileName;
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
        string secFileName = "t0111.sec";
        secFileName = examplesPath + secFileName;
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

        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentData>("Test1 - Loads the alignment of two identical sequences .",
                &TestAlignmentData::testAlignmentData_A));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentData>("Test2 -Verifies the matching sequences .",
                &TestAlignmentData::testAlignmentData_B));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentData>("Test3 - evaluating match using same positions.",
                &TestAlignmentData::testAlignmentData_C));
        suiteOfTests->addTest(new CppUnit::TestCaller<TestAlignmentData>("Test4 - evaluating match using different positions.",
                &TestAlignmentData::testAlignmentData_D));

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
        // Loads the alignment of two identical sequences 

        CPPUNIT_ASSERT((testAlignmentData->n == 4)&&(testAlignmentData->name2 == testAlignmentData->name1));
    }

    void testAlignmentData_B() {
        // Verifies the matching sequences 

        CPPUNIT_ASSERT((testAlignmentData->getSequence(1) == testAlignmentData->getSequence(2)));
    }

    void testAlignmentData_C() {
        //evaluating match using the same positions from the same sequence
        string path = getenv("VICTOR_ROOT");
        string examplesPath = path + "Align2/Tests/data/";
        string outputFileName = "output1.fasta";
        string spectedoutputFileName = "SpectedOutput.fasta";
        spectedoutputFileName=examplesPath +spectedoutputFileName;
        outputFileName = examplesPath + outputFileName;
        ofstream outputFile(outputFileName.c_str());
        ifstream specoutputFile(spectedoutputFileName.c_str());
        testAlignmentData->calculateMatch(1, 1, 1, 1);
        testAlignmentData->calculateMatch(1, 1, 2, 2);
        testAlignmentData->calculateMatch(1, 1, 3, 3);
        testAlignmentData->calculateMatch(1, 1, 4, 4);
        testAlignmentData->calculateMatch(1, 1, 5, 5);
        testAlignmentData->calculateMatch(1, 1, 6, 6);
        testAlignmentData->outputMatch(outputFile);
        ifstream routputFile(outputFileName.c_str());
        CPPUNIT_ASSERT((readLine(routputFile)==readLine(specoutputFile)));
    }

    void testAlignmentData_D() {
         //evaluating match using different positions from the same sequence
        string path = getenv("VICTOR_ROOT");
        string examplesPath = path + "Align2/Tests/data/";
        string outputFileName = "output2.fasta";
        string spectedoutputFileName = "SpectedOutput.fasta";
        spectedoutputFileName=examplesPath +spectedoutputFileName;
        outputFileName = examplesPath + outputFileName;
        ofstream outputFile(outputFileName.c_str());
        ifstream specoutputFile(spectedoutputFileName.c_str());
        testAlignmentData->calculateMatch(1, 2, 1, 1);
        testAlignmentData->calculateMatch(1, 2, 2, 2);
        testAlignmentData->calculateMatch(1, 2, 3, 3);
        testAlignmentData->calculateMatch(1, 2, 4, 4);
        testAlignmentData->calculateMatch(1, 2, 5, 5);
        testAlignmentData->calculateMatch(1, 2, 6, 6);
        testAlignmentData->outputMatch(outputFile);
        ifstream routputFile(outputFileName.c_str());
        CPPUNIT_ASSERT((readLine(routputFile)!=readLine(specoutputFile)));
    }
};
