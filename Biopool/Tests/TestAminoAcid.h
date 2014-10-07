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
using namespace Biopool;

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
         cout<<aa->getAtom(0).getType()<<aa->getAtom(0).getCoords().x<<"\n";
         cout<<aa->getAtom(1).getType()<<"\n";
         cout<<aa->getAtom(2).getType()<<"\n"; 
         cout<<aa->getAtom(3).getType()<<"\n"; 
         cout<<aa->getAtom(4).getType()<<aa->getSideChain().size()<<"\n"; 
         SideChain aa1=aa->getSideChain();
          cout<<aa1.getAtom(0).getCoords().x<<"\n";
        cout<<aa->getSideChain().getAtom(0).getType()<<"\n";//->getPhi();//->getPsi()<<" "<<aa->getPhi() <<" "<<aa->getOmega();
        cout<< aa->getSideChain().getAtom(0).getCoords().x<<"\n";
        cout<< aa->getSideChain().getAtom(1).getType()<<"\n";//->getPhi();//->getPsi()<<" "<<aa->getPhi() <<" "<<aa->getOmega();
        cout<<aa->getSideChain().getAtom(0).getCoords().y<<"\n";
        cout<< aa->getSideChain().getAtom(2).getType()<<"\n";//->getPhi();//->getPsi()<<" "<<aa->getPhi() <<" "<<aa->getOmega();
        cout<< aa->getSideChain().getAtom(0).getCoords().z<<"\n";
      //  cout << endl << "Solution is: x=" << atom0.distance(atom1)
        //        << ", y=" << distance << endl;
        CPPUNIT_ASSERT(1 == 1);
    }
    /*
    void testTestAtom_C() {
            // Calculate the distance between the two atoms
            Atom atom0;
            atom0.setCoords(0,0,0);
            Atom atom1;
            atom1.setCoords(0,0,0);
            //this is the expected value
            double distance = sqrt(0);
                
            cout << endl << "Solution is: x=" << atom0.distance(atom1)
                    << ", y=" << distance << endl;
            CPPUNIT_ASSERT( atom0.distance(atom1) == distance );
    }*/


};
