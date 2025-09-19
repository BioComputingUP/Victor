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

/* 
 * Author: Stefano Cecconello
 *
 */


#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include <SuperImpositor.h>

#include <PdbLoader.h>

using namespace std;
using namespace Victor::Biopool;

class TestSuperImpositor : public CppUnit::TestFixture {
private:
    Spacer *testSuperImpositor;
public:

    TestSuperImpositor() : testSuperImpositor(NULL) {
    }

    virtual ~TestSuperImpositor() {
        delete testSuperImpositor;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestSuperImpositor");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSuperImpositor>("Test1 - calculate rmsd.",
                &TestSuperImpositor::testTestSuperImpositor_A));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSuperImpositor>("Test2 - calculate maxsub.",
                &TestSuperImpositor::testTestSuperImpositor_B));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSuperImpositor>("Test3 - calculate maxsub.",
                &TestSuperImpositor::testTestSuperImpositor_C));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSuperImpositor>("Test4 - calculate maxsub.",
                &TestSuperImpositor::testTestSuperImpositor_D));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestSuperImpositor>("Test5 - calculate rotation.",
                &TestSuperImpositor::testTestSuperImpositor_E));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {




    }

    /// Teardown method

    void tearDown() {
    }

protected:

    void testTestSuperImpositor_A() {
        string path = getenv("VICTOR_ROOT");
        string inputFile1 = path + "Biopool/Tests/data/TestAlign1.pdb";
        string inputFile2 = path + "Biopool/Tests/data/TestAlign2.pdb";

        ifstream inFile1(inputFile1.c_str());
        ifstream inFile2(inputFile2.c_str());
        if (!inFile1)
            ERROR("First file not found.", exception);
        if (!inFile2)
            ERROR("Second File not found.", exception);
        PdbLoader pl1(inFile1);
        PdbLoader pl2(inFile2);
        Protein* prot1 = new Protein();
        Protein* prot2 = new Protein();
        pl1.setNoVerbose();
        pl2.setNoVerbose();
        pl1.setNoHAtoms();
        pl2.setNoHAtoms();
        prot1->load(pl1);
        prot2->load(pl2);
        string method = "Kabsch";
        SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, method);
        superImpositor->calculateRMSD();
        double rmsd = superImpositor->getRmsdValue();
        cout << "\n\nThe rmsd value between the two protein is: " << rmsd << "\n";
        CPPUNIT_ASSERT((rmsd - 21.732) < EPSILON && (-(rmsd - 21.732)) < EPSILON);
    }

    void testTestSuperImpositor_B() {
        string path = getenv("VICTOR_ROOT");
        string inputFile1 = path + "Biopool/Tests/data/TestAlign1.pdb";
        string inputFile2 = path + "Biopool/Tests/data/TestAlign2.pdb";

        ifstream inFile1(inputFile1.c_str());
        ifstream inFile2(inputFile2.c_str());
        if (!inFile1)
            ERROR("First file not found.", exception);
        if (!inFile2)
            ERROR("Second File not found.", exception);
        PdbLoader pl1(inFile1);
        PdbLoader pl2(inFile2);
        Protein* prot1 = new Protein();
        Protein* prot2 = new Protein();
        pl1.setNoVerbose();
        pl2.setNoVerbose();
        pl1.setNoHAtoms();
        pl2.setNoHAtoms();
        prot1->load(pl1);
        prot2->load(pl2);
        string method = "Kabsch";

        SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, method);

        std::vector<std::pair<int, int> > vectorSet;

        for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
            vectorSet.push_back(std::make_pair(i, i));
        }

        superImpositor->calculateMaxSub(3.5, vectorSet, 'y');
        double maxSub = superImpositor->getMaxsubValue();
        cout << "\n\nThe maxsub value between the two protein is: " << maxSub << "\n";
        CPPUNIT_ASSERT((maxSub - 0.626291) < EPSILON && (-(maxSub - 0.626291)) < EPSILON);
    }

    void testTestSuperImpositor_C() {
        string path = getenv("VICTOR_ROOT");
        string inputFile1 = path + "Biopool/Tests/data/TestAlign1.pdb";
        string inputFile2 = path + "Biopool/Tests/data/TestAlign2.pdb";

        ifstream inFile1(inputFile1.c_str());
        ifstream inFile2(inputFile2.c_str());
        if (!inFile1)
            ERROR("First file not found.", exception);
        if (!inFile2)
            ERROR("Second File not found.", exception);
        PdbLoader pl1(inFile1);
        PdbLoader pl2(inFile2);
        Protein* prot1 = new Protein();
        Protein* prot2 = new Protein();
        pl1.setNoVerbose();
        pl2.setNoVerbose();
        pl1.setNoHAtoms();
        pl2.setNoHAtoms();
        prot1->load(pl1);
        prot2->load(pl2);
        string method = "Kabsch";

        SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, method);

        std::vector<std::pair<int, int> > vectorSet;

        for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
            vectorSet.push_back(std::make_pair(i, i));
        }

        superImpositor->calculateGdt(vectorSet);
        double gdt = superImpositor->getGdtValue();
        cout << "\n\nThe gdt value between the two protein is: " << gdt << "\n";
        CPPUNIT_ASSERT(((gdt < 0.6725) + 0.1) && (gdt > (0.6725 - 0.1)));
    }

    void testTestSuperImpositor_D() {
        string path = getenv("VICTOR_ROOT");
        string inputFile1 = path + "Biopool/Tests/data/TestAlign1.pdb";
        string inputFile2 = path + "Biopool/Tests/data/TestAlign2.pdb";

        ifstream inFile1(inputFile1.c_str());
        ifstream inFile2(inputFile2.c_str());
        if (!inFile1)
            ERROR("First file not found.", exception);
        if (!inFile2)
            ERROR("Second File not found.", exception);
        PdbLoader pl1(inFile1);
        PdbLoader pl2(inFile2);
        Protein* prot1 = new Protein();
        Protein* prot2 = new Protein();
        pl1.setNoVerbose();
        pl2.setNoVerbose();
        pl1.setNoHAtoms();
        pl2.setNoHAtoms();
        prot1->load(pl1);
        prot2->load(pl2);
        string method = "Kabsch";

        SuperImpositor* superImpositor = new SuperImpositor(prot1, prot2, method);

        std::vector<std::pair<int, int> > vectorSet;

        for (unsigned int i = 0; i < superImpositor->getSet1()->size(); i++) {
            vectorSet.push_back(std::make_pair(i, i));
        }

        superImpositor->calculateTMScore(vectorSet);
        double TMScore = superImpositor->getTMScoreValue();
        cout << "\n\nThe TMScore value between the two protein is: " << TMScore << "\n";
        CPPUNIT_ASSERT(((TMScore < 0.7408) + 0.1) && (TMScore > (0.7408 - 0.1)));
    }

    void testTestSuperImpositor_E() {
        Eigen::Matrix3Xd set1(3, 2);

        set1(0, 0) = 1;
        set1(0, 1) = 2;
        set1(1, 0) = 0;
        set1(1, 1) = 0;
        set1(2, 0) = 0;
        set1(2, 1) = 0;

        Eigen::Affine3d* input = new Eigen::Affine3d();
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
        I(0, 0) = -1;
        input->linear() = I;
        input->translation() = Eigen::Vector3d::Zero();
        SuperImpositor::calculateRotation(set1, input);
        cout << "\n\nThe rototransled coordinates are:\n";
        cout << set1(0, 0) << "," << set1(0, 1) << "\n";
        cout << set1(1, 0) << "," << set1(1, 1) << "\n";
        cout << set1(2, 0) << "," << set1(2, 1) << "\n";
        CPPUNIT_ASSERT(set1(0, 0) == -1 && set1(0, 1) == -2 && set1(1, 0) == 0 && set1(1, 1) == 0 && set1(2, 0) == 0 && set1(2, 1) == 0);
    }
};