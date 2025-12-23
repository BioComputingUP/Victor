/*  This file is part of Victor.
    Victor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    Victor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with Victor.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <iostream>
#include <math.h>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>
#include <StructuralSuperimposition.h>

#include <PdbLoader.h>

using namespace std;
using namespace Victor::Biopool;

class TestStructuralSuperimposition : public CppUnit::TestFixture 
{
    private:
        
    public:

        // SETUP METHOD
        void setUp() 
        {}

        /// TEARDOWN METHOD

        void tearDown() 
        {}


        void testScores()
        {
            int L=4;
            double d0=3.5;

            double computedTMscore=0.932606;
            double computedGDTscore=0.767476;
            double computedRMSDscore=1.58475;
            double computedMaxSubscore=0.857753;

            char *pPath = getenv("VICTOR_ROOT");
            if(pPath==NULL)
            {
                ERROR("\"VICTOR_ROOT\" not set.", exception);
            }
            string path(pPath);

            string inputFile1= path + "Biopool/Tests/data/15C8_H_input.pdb"; 
            string inputFile2= path + "Biopool/Tests/data/25C8_H_input.pdb";

            ifstream inFile1(inputFile1.c_str());
            
            ifstream inFile2(inputFile2.c_str());

            if(!inFile1 || !inFile2)
                ERROR("Files not found.", exception);
            
            Protein prot1, prot2;   
            Spacer *sp1, *sp2;

            PdbLoader pl1(inFile1);
            PdbLoader pl2(inFile2);

            pl1.setNoVerbose();
            pl2.setNoVerbose();
            pl1.setNoHAtoms();
            pl2.setNoHAtoms();

            
            prot1.load(pl1);       
            prot2.load(pl2);       
            
        
            vector<char> chain1 = prot1.getAllChains();
            vector<char> chain2 = prot2.getAllChains();
            
            sp1=prot1.getSpacer(chain1[0]);
            sp2=prot2.getSpacer(chain2[0]);
            
            StructuralSuperimposition s1, s2= StructuralSuperimposition();
            
            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A=StructuralSuperimposition::convertSpacerToEigen(sp1);

            Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> B=StructuralSuperimposition::convertSpacerToEigen(sp2);
            
            StructuralSuperimposition ssi= StructuralSuperimposition();

            double rmsdscore= ssi.getRmsdScore(A,B,L,d0);

            double tmscore= ssi.getTMScore(A, B, L);

            double msscore= ssi.getMaxSubScore(A, B, L, d0);
                
            double gdttsscore= ssi.getGDTTSScore(A,B,L);
		
            double threshold=0.0001;
            CPPUNIT_ASSERT( fabs(computedTMscore - tmscore) <= threshold );
            CPPUNIT_ASSERT( fabs(computedGDTscore - gdttsscore) <= threshold );
            CPPUNIT_ASSERT( fabs(computedMaxSubscore - msscore) <= threshold );
            CPPUNIT_ASSERT( fabs(computedRMSDscore - rmsdscore) <= threshold );
        }

        static CppUnit::Test *suite()
        {
            CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite( "TestStructuralSuperimposition" );
            suiteOfTests->addTest( new CppUnit::TestCaller<TestStructuralSuperimposition>("testScores",&TestStructuralSuperimposition::testScores ) );
            return suiteOfTests;
        }
    };
