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
				&TestSolvationPotential::testSolvationPotential_calculateSolvationSpacer ));
                
                suiteOfTests->addTest(new CppUnit::TestCaller<TestSolvationPotential>("Test2 - energy of a residue.",
				&TestSolvationPotential::testSolvationPotential_calculateSolvationAA ));
                
                suiteOfTests->addTest(new CppUnit::TestCaller<TestSolvationPotential>("Test3 - energy propensity difference have to be positive.",
				&TestSolvationPotential::testSolvationPotential_increasingPropensity ));
                
		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
    
	void testSolvationPotential_calculateSolvationSpacer() {
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
                
                //calculating new energy
                SolvationPotential pot;
		CPPUNIT_ASSERT( pot.calculateEnergy(sp) - 0.869 < 0.01 );
	}
        
        void testSolvationPotential_calculateSolvationAA() {
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
                
                //calculating new energy for each residue
                SolvationPotential pot;
                double en0 = pot.calculateEnergy(sp.getAmino(0), ALA, sp);
                double en1 = pot.calculateEnergy(sp.getAmino(1), ALA, sp);
		CPPUNIT_ASSERT( en0 == en1 && ((en0 + en1 - 0.869) < 0.01 ) );
	}

	void testSolvationPotential_increasingPropensity() {
                bool increasingEnergy = true;
                SolvationPotential s;
                
                //max propensity should be always greater the the min
                for(int i = 0; i < AminoAcid_CODE_SIZE - 1; ++i) {
                    increasingEnergy = increasingEnergy &
                            (s.pReturnMaxPropensity(static_cast<AminoAcidCode>(i))
                            > s.pReturnMinPropensity(static_cast<AminoAcidCode>(i)) );
                }

		CPPUNIT_ASSERT( increasingEnergy );
	}
        
};
