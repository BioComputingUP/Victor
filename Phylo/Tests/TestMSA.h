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
#include <sstream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>
#include <phylo.h>

using namespace std;
using namespace Victor;
using namespace Victor::Phylo;

/** @brief Test on phylogenetic trees.
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class TestMSA: public CppUnit::TestFixture {
 public:
    TestMSA() {
    }



    /** @brief Builds a suite of tests. */
    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestMSA");

        suiteOfTests->addTest(
            new CppUnit::TestCaller<TestMSA>("Test1 - Feng-Doolittle.",
            &TestMSA::Test_A)
        );

        suiteOfTests->addTest(
            new CppUnit::TestCaller<TestMSA>("Test2 - ClustalW.",
            &TestMSA::Test_B)
        );

        return suiteOfTests;
    }



    /** @brief Setup method. */
    void setUp() {
    }



    /** @brief Teardown method. */
    void tearDown() {
    }



protected:
    /** @brief Test for Feng-Doolittle. */
    void Test_A() {
        stringstream input;

        input << ">Seq1\nMAAAAATLRGAMVGPRGAGLP\n"
              << ">Seq2\nMAAAAASLRGVVLGPRGAGL\n"
              << ">Seq3\nMTEFKAGSAKKGATLFKTRCL\n"
              << ">Seq4\nMAAAAASLRRTVLGPRGVGLPGASAPGLL\n"
              << ">Seq5\nMFSQKLLANGKLLSKLAIVSGVVG\n"
              << "\n";

        Alignment alignment;
        alignment.loadFasta(input);
        FengDoolittle feng_doolittle;
        MultipleAlignment MSA  = feng_doolittle(alignment);

        CPPUNIT_ASSERT(MSA.getConsensus() == "MAAFKAAASLRGAVLGPRGAGLAGASAPGLG");
    }



    /** @brief Test for ClustalW. */
    void Test_B() {
        stringstream input;

        input << ">Seq1\nMAAAAATLRGAMVGPRGAGLP\n"
              << ">Seq2\nMAAAAASLRGVVLGPRGAGL\n"
              << ">Seq3\nMTEFKAGSAKKGATLFKTRCL\n"
              << ">Seq4\nMAAAAASLRRTVLGPRGVGLPGASAPGLL\n"
              << ">Seq5\nMFSQKLLANGKLLSKLAIVSGVVG\n"
              << "\n";

        Alignment alignment;
        alignment.loadFasta(input);

        SubMatrix matrix(SubstitutionMatrix(SubstitutionMatrix::BLOSUM62).asSubMatrix());
        AGPFunction gap_function(10.0, 0.2);
        IdentityPercentage build_matrix(matrix, gap_function);
        NJ build_tree;
        ClustalW clustalw(build_matrix, build_tree, ClustalW::BLOSUM);

        MultipleAlignment MSA  = clustalw(alignment);
        
        CPPUNIT_ASSERT(MSA.getConsensus() == "MAAAAASLRGKVLGPRGAGLPGASAGCLG");
    }
};