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
class TestClustering: public CppUnit::TestFixture {
 public:
    TestClustering() {
    }



    /** @brief Builds a suite of tests. */
    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestClustering");

        suiteOfTests->addTest(
            new CppUnit::TestCaller<TestClustering>("Test1 - UPGMA.",
            &TestClustering::Test_A)
        );

        suiteOfTests->addTest(
            new CppUnit::TestCaller<TestClustering>("Test2 - Nighbor joining.",
            &TestClustering::Test_B)
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
    /** @brief Test for UPGMA. */
    void Test_A() {
        DistanceMatrix D;
        D.addOTU("Human")
         .addOTU("Chimp")
         .addOTU("Gorilla")
         .addOTU("Orang");

        D("Human",   "Chimp",   0.095);
        D("Human",   "Gorilla", 0.113);
        D("Human",   "Orang",   0.183);
        D("Chimp",   "Gorilla", 0.118);
        D("Chimp",   "Orang",   0.201);
        D("Gorilla", "Orang",   0.195);

        UPGMA upgma;
        RootedTree tree = upgma(D).asRootedTree();
        const string result = tree.asNewick();

        CPPUNIT_ASSERT(result == "(((Chimp:0.0475,Human:0.0475):0.01025,Gorilla:0.05775):0.03875,Orang:0.0965)");
    }



    /** @brief Test for neighbor joining. */
    void Test_B() {
        DistanceMatrix D;
        D.addOTU("A")
         .addOTU("B")
         .addOTU("C")
         .addOTU("D")
         .addOTU("E")
         .addOTU("F")
         .addOTU("G")
         .addOTU("H");

        D("A", "B",  7);
        D("A", "C",  8);
        D("A", "D", 11);
        D("A", "E", 13);
        D("A", "F", 16);
        D("A", "G", 13);
        D("A", "H", 17);
        D("B", "C",  5);
        D("B", "D",  8);
        D("B", "E", 10);
        D("B", "F", 13);
        D("B", "G", 10);
        D("B", "H", 14);
        D("C", "D",  5);
        D("C", "E",  7);
        D("C", "F", 10);
        D("C", "G",  7);
        D("C", "H", 11);
        D("D", "E",  8);
        D("D", "F", 11);
        D("D", "G",  8);
        D("D", "H", 12);
        D("E", "F",  5);
        D("E", "G",  6);
        D("E", "H", 10);
        D("F", "G",  9);
        D("F", "H", 13);
        D("G", "H",  8);

        NJ nj;
        UnrootedTree tree = nj(D).asUnrootedTree();
        const string result = tree.asNewick();

        CPPUNIT_ASSERT(result == "(((E:1,F:4):2,(G:2,H:6):1):1.5,(D:3,(C:1,((A:5,B:2):1):1):1):0.5)");
    }
};