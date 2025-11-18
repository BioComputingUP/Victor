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
class TestDistanceMatrix: public CppUnit::TestFixture {
 public:
    TestDistanceMatrix() {
    }



    /** @brief Builds a suite of tests. */
    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestDistanceMatrix");

        suiteOfTests->addTest(
            new CppUnit::TestCaller<TestDistanceMatrix>("Test1 - Builds a matrix.",
            &TestDistanceMatrix::Test_A)
        );

        suiteOfTests->addTest(
            new CppUnit::TestCaller<TestDistanceMatrix>("Test2 - Get maximum value.",
            &TestDistanceMatrix::Test_B)
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
    /** @brief Test for creation. */
    void Test_A() {
        vector<string> OTUs;
        OTUs.push_back("A");
        OTUs.push_back("B");
        OTUs.push_back("C");

        const double ab = 0.1, ac = 9.5, bc = 0.333;

        DistanceMatrix d;
        d.addOTU(OTUs[0])
         .addOTU(OTUs[1])
         .addOTU(OTUs[2]);

        d("A", "B", ab);
        d("A", "C", ac);
        d("B", "C", bc);

        CPPUNIT_ASSERT(
            d.getSize() == 3 &&
            d("A", "B") == ab &&
            d("B", "C") == bc &&
            d("A", "C") == ac
        );
    }



    /** @brief Test for get values. */
    void Test_B() {
        vector<string> OTUs;
        OTUs.push_back("A");
        OTUs.push_back("B");
        OTUs.push_back("C");

        const double ab = 0.1, ac = 9.5, bc = 0.333;

        DistanceMatrix d;
        d.addOTU(OTUs[0])
         .addOTU(OTUs[1])
         .addOTU(OTUs[2]);

        d("A", "B", ab);
        d("A", "C", ac);
        d("B", "C", bc);

        CPPUNIT_ASSERT(
            d.getMaximum() == 9.5 &&
            d.getMinimum() == 0.1
        );
    }
};