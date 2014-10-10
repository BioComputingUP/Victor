/*
 * TestAtom.cpp
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

#include <Spacer.h>

#include <PdbLoader.h>

using namespace std;
using namespace Victor::Biopool;

class TestSpacer : public CppUnit::TestFixture {
private:
    Spacer *testSpacer;
public:

    TestSpacer() : testSpacer(NULL) {
    }

    virtual ~TestSpacer() {
        delete testSpacer;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestSpacer");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSpacer>("Test1 - Loading a chain from pdb.",
                &TestSpacer::testTestSpacer_A));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSpacer>("Test2 - Gaps in pdb.",
                &TestSpacer::testTestSpacer_B));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSpacer>("Test3 - loading amino acids from pdb without chain.",
                &TestSpacer::testTestSpacer_C));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {




    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testTestSpacer_A() {
        
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
        cout << "\n\nNumber of atoms in the first amino acid" << sp->getAmino(0).size() << "\n";
        cout << "Number of atoms in the second amino acid" << sp->getAmino(1).size() << "\n";
        cout << "Number of atoms in the third amino acid" << sp->getAmino(2).size() << "\n";
        CPPUNIT_ASSERT((sp->getAmino(0).size() == 7) && (sp->getAmino(1).size() == 5) && (sp->getAmino(2).size() == 11));
    }

    void testTestSpacer_B() {
         
        string path = getenv("VICTOR_ROOT");
        string inputFile = path + "Biopool/Tests/data/test1.pdb";

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
        cout << "\n\nAmino Acid 3 is gap?" << sp->isGap(3) << "\n";
        cout << "Amino Acid 5 is gap?" << sp->isGap(5) << "\n";
        cout << "Amino Acid 8 is gap?" << sp->isGap(8) << "\n";
        CPPUNIT_ASSERT((sp->isGap(5) == 1)&&(sp->isGap(4) == 1));
    }

    void testTestSpacer_C() {
        string path = getenv("VICTOR_ROOT");
        string inputFile = path + "Biopool/Tests/data/test1.pdb";

        ifstream inFile(inputFile.c_str());
        Spacer* sp;
        if (!inFile)
            ERROR("File not found.", exception);
        PdbLoader pl(inFile);
        Protein prot;
        pl.setNoVerbose();
        pl.setNoHAtoms();
        prot.load(pl);
        unsigned int i = 0;
        sp = prot.getSpacer(i);
        cout << "\n\nNumber of atoms in the first amino acid" << sp->getAmino(0).size() << "\n";
        cout << "Number of atoms in the second amino acid" << sp->getAmino(1).size() << "\n";
        cout << "Number of atoms in the third amino acid" << sp->getAmino(2).size() << "\n";
        CPPUNIT_ASSERT((sp->getAmino(0).size() == 7) && (sp->getAmino(1).size() == 5) && (sp->getAmino(2).size() == 11));
    }


};
