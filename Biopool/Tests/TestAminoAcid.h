/*
 * TestAminoAcid.cpp
 *
 *  Created on: Oct 7th, 2014
 *      Author: Layla Hirsh
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>
#include <AminoAcid.h>
#include <Atom.h>

using namespace std;
using namespace Victor::Biopool;

class TestAminoAcid : public CppUnit::TestFixture {
private:
    AminoAcid *testAminoAcid;
public:

    TestAminoAcid() : testAminoAcid(NULL) {
    }

    virtual ~TestAminoAcid() {
        delete testAminoAcid;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestAminoAcid");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestAminoAcid>("Test1 - Initialize the object.",
                &TestAminoAcid::testTestAminoAcid_A));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestAminoAcid>("Test2 - Load from file.",
                &TestAminoAcid::testTestAminoAcid_B));
        /*        
                suiteOfTests->addTest(new CppUnit::TestCaller<TestAtom>("Test3 - rotation.",
                                &TestAminoAcid::testTestAminoAcid_C ));
         */
        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testTestAminoAcid_A() {
        //Initialize
        AminoAcid* aa0 = new AminoAcid();
        //Setting the constructor values by default
        AminoAcid aa1;
        aa1.setPhi(999);
        aa1.setPsi(999);
        aa1.setOmega(999);
        aa1.setState(COIL);
        aa1.setType("XXX");
        CPPUNIT_ASSERT((aa1.getPsi() == aa0->getPsi()) && (aa1.getPhi() == aa0->getPhi()) &&(aa1.getOmega() == aa0->getOmega()) &&(aa1.getState() == aa0->getState()));
    }

    void testTestAminoAcid_B() {
        // Calculate the distance between the two atoms
        //cout << "Start" << endl;

        AminoAcid *aa=new AminoAcid();
        string path = getenv("VICTOR_ROOT");
        string dataPath = path +  "Biopool/Tests/data/ASP.pdb";
        ifstream inFile(dataPath.c_str());
        if (!inFile)
            ERROR("File not found.", exception);
        XyzLoader il(inFile);
        aa->load(il);
        AminoAcid *aa1=new AminoAcid();
        Atom atom1,atom2;
        atom1.setCoords(55.895,39.622,14.803);
        atom1.setCode(N);
        atom2.setCoords(56.598,38.679,13.941);
        atom2.setCode(CA);
        aa1->addAtom(atom1);
        aa1->addAtom(atom2);
        CPPUNIT_ASSERT((aa->getAtom(0).getCoords()== aa1->getAtom(0).getCoords())&&(aa->getAtom(1).getCoords()== aa1->getAtom(1).getCoords())&&
                (aa->getAtom(0).getCode()== aa1->getAtom(0).getCode())&&(aa->getAtom(1).getCode()== aa1->getAtom(1).getCode()));
    }



};
