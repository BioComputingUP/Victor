/*
 * TestRapdfPotential.cpp
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

#include <PdbLoader.h>
#include <XyzLoader.h>
#include <Protein.h>
#include <RapdfPotential.h>

using namespace std;
using namespace Biopool;

class TestRapdfPotential : public CppUnit::TestFixture {
private:
	RapdfPotential *testRapdfPotential;
        string path;
public:
	TestRapdfPotential() : testRapdfPotential(NULL), path(getenv("VICTOR_ROOT")) {}
	virtual ~TestRapdfPotential() {
		delete testRapdfPotential;
	}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestRapdfPotential");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestRapdfPotential>("Test1 - energy in a spacer.",
				&TestRapdfPotential::testRapdfPotential_spacer ));
                
                suiteOfTests->addTest(new CppUnit::TestCaller<TestRapdfPotential>("Test2 - energy between 2 very distant atoms.",
				&TestRapdfPotential::testRapdfPotential_calcEnergyDistantAtoms ));
                
                suiteOfTests->addTest(new CppUnit::TestCaller<TestRapdfPotential>("Test3 - N vs N in ALA for distance < 3.",
				&TestRapdfPotential::testRapdfPotential_N2Nenergy ));
                
                suiteOfTests->addTest(new CppUnit::TestCaller<TestRapdfPotential>("Test4 - energy aa vs aa.",
				&TestRapdfPotential::testRapdfPotential_calcEnergyAA2AA ));

		return suiteOfTests;
	}

	/// Setup method
	void setUp() {}

	/// Teardown method
	void tearDown() {}

protected:
	void testRapdfPotential_spacer() {
                Spacer s;
                string p = path + "Energy/Tests/data/test.pdb";
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
                for(unsigned int i = 0; i < sp.sizeAmino(); ++i) {
                    for(unsigned int j = 0; j < sp.getAmino(i).size(); ++j) {
                        //set atoms to N, just to simplify the calculation
                        sp.getAmino(i)[j].setCode(N);
                        
                    }
                }
                
                //calculating new energy
                RapdfPotential r;
                r.calculateEnergy(sp);
		CPPUNIT_ASSERT( r.calculateEnergy(sp) == ((8 * 4 + 4 * 4) * -1.83) );
	}

	void testRapdfPotential_calcEnergyDistantAtoms() {
                Atom atom0;
                atom0.setCoords(0,0,0);
                atom0.setType("N");
                
                //atom1 is very far from atom0
                Atom atom1;
                atom1.setCoords(200,0,0);
                atom1.setType("N");
                
                RapdfPotential r;
                //expected energy is 0, due to high distance
		CPPUNIT_ASSERT( r.calculateEnergy(atom0, atom1, "ALA", "ALA") == 0 );
	}
        
        void testRapdfPotential_N2Nenergy() {
                Atom atom0;
                atom0.setCoords(0,0,0);
                atom0.setType("N");
                
                Atom atom1;
                atom1.setCoords(1,0,0);
                atom1.setType("N");
                
                RapdfPotential r;
                //expected energy according to the ram.paf file is -1.83
		CPPUNIT_ASSERT( r.calculateEnergy(atom0, atom1, "ALA", "ALA") == -1.83 );
	}
        
        void testRapdfPotential_calcEnergyAA2AA() {
                AminoAcid aa0;
                AminoAcid aa1;
                string dataPath0 = path +  "Energy/Tests/data/ALA0.pdb";
                ifstream inFile0(dataPath0.c_str());
                if (!inFile0)
                    ERROR("File not found.", exception);
                XyzLoader il0(inFile0);
                aa0.load(il0);

                //ALA1 is identical to ALA0, except for a rotation.
                string dataPath1 = path +  "Energy/Tests/data/ALA1.pdb";
                ifstream inFile1(dataPath1.c_str());
                if (!inFile1)
                    ERROR("File not found.", exception);
                XyzLoader il1(inFile1);
                aa1.load(il1);
                RapdfPotential r;
                
                //the expected internal energy is just the energy between 2 atom multiplied by the squared number of residues
		CPPUNIT_ASSERT( r.calculateEnergy(aa0, aa1)
                        == r.calculateEnergy(aa0.getAtom(0), aa1.getAtom(0), "ALA", "ALA") * aa0.size() * aa1.size() );
	}
};
