/*
 * TestLoopModel.cpp
 *
 *  Created on: Oct 6th, 2014
 *      Author: Manuel Giollo
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <LoopModel.h>
#include <PdbLoader.h>
#include <Spacer.h>

using namespace std;
using namespace Victor::Lobo;

class TestLoopModel : public CppUnit::TestFixture {
private:
    LoopModel *testLoopModel;
public:

    TestLoopModel() : testLoopModel(NULL) {
    }

    virtual ~TestLoopModel() {
        delete testLoopModel;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestLoopModel");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestLoopModel>("Test1 - calculates rms value.",
                &TestLoopModel::testTestLoopModel_A));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestLoopModel>("Test2 - verifies the initialized values.",
                &TestLoopModel::testTestLoopModel_B));


        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testTestLoopModel_A() {
        // Calculate the RMs
         string path = getenv("VICTOR_ROOT");
        string inputFile = path + "Biopool/Tests/data/test.pdb";

        ifstream inFile(inputFile.c_str());
        Spacer* sp;
        if (!inFile)
            ERROR("File not found.", exception);
        PdbLoader pl(inFile);
        Protein prot;
        pl.setNoVerbose();
        pl.setNoHAtoms();
        prot.load(pl);
        sp = prot.getSpacer('A');
        LoopModel lm;
        unsigned int start=1;
        unsigned int end=2;
        double rms= lm.calculateRms(*sp,start,end,*sp);
        CPPUNIT_ASSERT( rms >= 6.8 );
    }

    void testTestLoopModel_B() {
        // verifies the initialized values
        LoopModel lm;
        string tableFile = getenv("VICTOR_ROOT");
        if (tableFile.length() < 3)
            ERROR("Environment variable VICTOR_ROOT was not found.", exception);
        string path = "data/aa0.lt";
        tableFile += path;
        string tableFileFromConstructor= lm.getTableFileName()[0];
        double EndRmsFromConst = lm.getENDRMS_WEIGHT(0);
        CPPUNIT_ASSERT((EndRmsFromConst==125 )&& (tableFileFromConstructor==tableFile));
    }
};
