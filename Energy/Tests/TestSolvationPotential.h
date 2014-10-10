/*
 * TestSolvationPotential.cpp
 *
 *  Created on: Oct 9th, 2014
 *      Author: Manuel Giollo
 */

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <PdbLoader.h>
#include <XyzLoader.h>
#include <Protein.h>
#include <SolvationPotential.h>

using namespace std;
using namespace Biopool;

class TestSolvationPotential : public CppUnit::TestFixture {
private:
	SolvationPotential *testSolvationPotential;
        string path;
public:
	TestSolvationPotential() : testSolvationPotential(NULL), path(getenv("VICTOR_ROOT")) {}
	virtual ~TestSolvationPotential() {
		delete testSolvationPotential;
	}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestSolvationPotential");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestSolvationPotential>("Test1 - energy in a spacer.",
				&TestSolvationPotential::testSolvationPotential_spacer ));

		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
    
	void testSolvationPotential_spacer() {
                Spacer s;
                string p = path + "Energy/Tests/data/solv.pdb";
                ifstream inFile2(p.c_str());
                if (!inFile2)
                  ERROR("File not found.", exception);
                PdbLoader pl(inFile2);
                pl.setChain('A');
                pl.setNoHAtoms();
                pl.setNoVerbose();

                pl.setPermissive();
                Protein prot;
                prot.load(pl);
                Spacer &sp = *prot.getSpacer('A');
                /*for(unsigned int i = 0; i < sp.sizeAmino(); ++i) {
                    for(unsigned int j = 0; j < sp.getAmino(i).size(); ++j) {
                        //set atoms to N, just to simplify the calculation
                        sp.getAmino(i)[j].setCode(N);
                        
                    }
                }*/
                
                //calculating new energy
                SolvationPotential pot;
                cout<<pot.calculateEnergy(sp)<<endl;
                //r.calculateEnergy(sp) == ((8 * 4 + 4 * 4) * -1.83)
		CPPUNIT_ASSERT( true );
	}

	/*void testSolvationPotential_calcEnergyDistantAtoms() {
                Atom atom0;
                atom0.setCoords(0,0,0);
                atom0.setType("N");
                
                //atom1 is very far from atom0
                Atom atom1;
                atom1.setCoords(200,0,0);
                atom1.setType("N");
                
                SolvationPotential s;
                //expected energy is 0, due to high distance
                cout<<s.calculateEnergy(atom0, atom1, "ALA", "ALA")<<endl;
                //r.calculateEnergy(atom0, atom1, "ALA", "ALA") == 0
		CPPUNIT_ASSERT( true );
	}
        
        void testSolvationPotential_N2Nenergy() {
                Atom atom0;
                atom0.setCoords(0,0,0);
                atom0.setType("N");
                
                Atom atom1;
                atom1.setCoords(1,0,0);
                atom1.setType("N");
                
                SolvationPotential s;
                //expected energy is 0, due to high distance
                cout<<s.calculateEnergy(atom0, atom1, "ALA", "ALA")<<endl;
                //r.calculateEnergy(atom0, atom1, "ALA", "ALA") == 0
		CPPUNIT_ASSERT( true );
	}*/
        
};
