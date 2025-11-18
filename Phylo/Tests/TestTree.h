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
class TestTree: public CppUnit::TestFixture {
 public:
    TestTree() {
    }



    /** @brief Builds a suite of tests. */
    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestTree");

        suiteOfTests->addTest(
            new CppUnit::TestCaller<TestTree>("Test1 - Parses a Newick string.",
            &TestTree::TestTree_A)
        );

        suiteOfTests->addTest(
            new CppUnit::TestCaller<TestTree>("Test2 - Count leaves.",
            &TestTree::TestTree_B)
        );

        suiteOfTests->addTest(
            new CppUnit::TestCaller<TestTree>("Test3 - Distances.",
            &TestTree::TestTree_C)
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
    /** @brief Test for parse/write. */
    void TestTree_A() {
        RootedTree tree;
        const string load_newick = "(A:0.5,B:0.25):0.36;";
        tree.parseNewick(load_newick);
        const string save_newick = tree.asNewick() + ";";

        CPPUNIT_ASSERT(load_newick == save_newick);
    }



    /** @brief Test for number of leaves. */
    void TestTree_B() {
        RootedTree tree;
        const string load_newick = "(A:0.5,B:0.25):0.36;";
        tree.parseNewick(load_newick);
        const size_t leaves = tree.getLeavesNumber();

        CPPUNIT_ASSERT(leaves == 2);
    }



    /** @brief Test for distance. */
    void TestTree_C() {
        RootedTree tree;
        const string load_newick = "(A:0.5,B:0.25):0.36;";
        tree.parseNewick(load_newick);
        const double max_distance = tree.getMaxDistance();

        CPPUNIT_ASSERT(max_distance == 0.86);
    }
};