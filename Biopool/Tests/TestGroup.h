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
#include <vglMath.h>
#include <Atom.h>
#include <XyzLoader.h>
#include <Group.h>

using namespace std;
using namespace Victor::Biopool;

class TestGroup : public CppUnit::TestFixture {
private:
    Group *testGroup;
public:

    TestGroup() : testGroup(NULL) {
    }

    virtual ~TestGroup() {
        delete testGroup;
    }

    static CppUnit::Test *suite() {
        CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestGroup");

        suiteOfTests->addTest(new CppUnit::TestCaller<TestGroup>("Test1 - initialize one atom in the group.",
                &TestGroup::testTestGroup_A));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestGroup>("Test2 - load atoms in the group.",
                &TestGroup::testTestGroup_B));

        suiteOfTests->addTest(new CppUnit::TestCaller<TestGroup>("Calculate the distance between the two loaded atoms.",
        	&TestGroup::testTestGroup_C));

        return suiteOfTests;
    }

    /// Setup method

    void setUp() {
    }

    /// Teardown method

    void tearDown() {
    }

    vg_ieee64 degreesToRadian(const vg_ieee64& deg) {
        return (deg / 180.0) * M_PI;
    }
protected:

    void testTestGroup_A() {
        // Initialize first atom in the group
        Atom atom0;
        atom0.setCoords(0, 0, 0);
        atom0.setType("N");
        Group group1(0, 0);
        Group group0;
        group0.setType("ALA");
        group0.addAtom(atom0);
        vgVector3<double> _coord(0.000, 0.000, 0.000);
        group0[0].setCoords(_coord);
        vgVector3<double> tmp(0, 0, 1);
        vgMatrix3<double> rotationMatrix = vgMatrix3<double>::createRotationMatrix(tmp, degreesToRadian(180));
        group0[0].addRot(rotationMatrix);
        cout << endl << "Solution for initialized group is: x=" << group1.getTrans().x << ", y=" << group1.getTrans().y << ", z=" << group1.getTrans().z << endl;
        cout << endl << "Solution for hard coded  group is: x=" << group0.getTrans().x << ", y=" << group0.getTrans().y << ", z=" << group0.getTrans().z << endl;
        cout << endl << "Solution for initialized group is:\nx.x=" << group1.getRot().x.x << ", y.x=" << group1.getRot().y.x << ", z.x=" << group1.getRot().z.x << endl;
        cout << endl << "x.y=" << group1.getRot().x.y << ", y.y=" << group1.getRot().y.y << ", z.y=" << group1.getRot().z.y << endl;
        cout << endl << "x.z=" << group1.getRot().x.z << ", y.z=" << group1.getRot().y.z << ", z.z=" << group1.getRot().z.z << endl;
        cout << endl << "Solution for hard coded group is:\nx.x=" << group0.getRot().x.x << ", y.x=" << group0.getRot().y.x << ", z.x=" << group0.getRot().z.x << endl;
        cout << endl << "x.y=" << group0.getRot().x.y << ", y.y=" << group0.getRot().y.y << ", z.y=" << group0.getRot().z.y << endl;
        cout << endl << "x.z=" << group0.getRot().x.z << ", y.z=" << group0.getRot().y.z << ", z.z=" << group0.getRot().z.z << endl;
        cout << endl << "Solution type for initialized group is:" << group1.getType() << endl;
        cout << endl << "Solution type for hard coded group is :" << group0.getType() << endl;
        CPPUNIT_ASSERT((group1.getTrans() == group0.getTrans())&& (group1.getRot() == group0.getRot()) &&(group0.getType() != group1.getType()));
    }

    void testTestGroup_B() {
        // load atoms in the group
        Group group1;
        
        /* Data in the file
            ALA
            1  N   0.000   0.000   0.000
            2  CA  0.000   0.000   1.460
            3  C   1.404   0.000   2.016
            4  O   1.732   0.688   2.968
            5  H   -0.397  -0.779  -0.525
            7  XA  -0.710  -1.266   1.974
         
         */
        string path = getenv("VICTOR_ROOT");
        string dataPath = path +  "Biopool/Tests/data/ALA.pdb";
        ifstream inFile(dataPath.c_str());
        if (!inFile)
            ERROR("File not found.", exception);

        XyzLoader il(inFile);
        group1.load(il);
        cout << endl << "Solution for C atom is: x=" << group1[2].getCoords().x << ", y=" << group1[2].getCoords().y << ", z=" << group1[2].getCoords().z << endl;
        CPPUNIT_ASSERT((group1[0].getType() == "N")&&(group1[1].getType() == "CA")&&(group1[2].getType() == "C") && (group1[3].getType() == "O") && (group1[4].getType() == "H")
                && (group1[5].getType() == "XA"));
    }

      void testTestGroup_C() {
              // Calculate the distance between the two loaded atoms
               Group group1;
        
        /* Data in the file
            ALA
            1  N   0.000   0.000   0.000
            2  CA  0.000   0.000   1.460
            3  C   1.404   0.000   2.016
            4  O   1.732   0.688   2.968
            5  H   -0.397  -0.779  -0.525
            7  XA  -0.710  -1.266   1.974
         
         */
        string path = getenv("VICTOR_ROOT");
        string dataPath = path +  "Biopool/Tests/data/ALA.pdb";
        ifstream inFile(dataPath.c_str());
        if (!inFile)
            ERROR("File not found.", exception);

        XyzLoader il(inFile);
        group1.load(il);
        double distance=group1[0].distance(group1[1]);
        //this is the expected value
        double calculatedDistance = sqrt(group1[1].getCoords().z*group1[1].getCoords().z-group1[0].getCoords().z*group1[0].getCoords().z);
        cout << endl << "Solution is: distance =" << distance << endl;
        cout << endl << "Solution is: distance =" <<  calculatedDistance<< endl;
        CPPUNIT_ASSERT(distance==calculatedDistance);
      }
    

};
