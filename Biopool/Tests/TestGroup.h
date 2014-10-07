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
#include <Group.h>

using namespace std;
using namespace Biopool;

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

        suiteOfTests->addTest(new CppUnit::TestCaller<TestGroup>("Test1 - initialize the group.",
                &TestGroup::testTestGroup_A));

        //suiteOfTests->addTest(new CppUnit::TestCaller<TestGroup>("Test2 - zero distance.",
        //		&TestGroup::testTestGroup_B ));

        //suiteOfTests->addTest(new CppUnit::TestCaller<TestGroup>("Test3 - rotation.",
        //		&TestGroup::testTestGroup_C ));

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
        // Calculate the distance between the two atoms
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
        cout << endl << "Solution type for initialized group is:"<<group1.getType()<< endl;
        cout << endl << "Solution type for hard coded group is :"<<group0.getType()<< endl;
                CPPUNIT_ASSERT((group1.getTrans() == group0.getTrans())&& (group1.getRot() == group0.getRot()) &&(group0.getType() != group1.getType()));
    }

    /*void testTestAtom_B() {
            // Calculate the distance between the two atoms
            Atom atom0;
            atom0.setCoords(0,0,0);
            Atom atom1;
            atom1.setCoords(0,0,0);
            //this is the expected value
            double distance = sqrt(0);
                
            cout << endl << "Solution is: x=" << atom0.distance(atom1)
                    << ", y=" << distance << endl;
            CPPUNIT_ASSERT( atom0.distance(atom1) == distance );
    }
        
    void testTestAtom_C() {
            // Calculate the distance between the two atoms
            Atom atom0;
            atom0.setCoords(0,0,0);
            Atom atom1;
            atom1.setCoords(0,0,0);
            //this is the expected value
            double distance = sqrt(0);
                
            cout << endl << "Solution is: x=" << atom0.distance(atom1)
                    << ", y=" << distance << endl;
            CPPUNIT_ASSERT( atom0.distance(atom1) == distance );
    }
     */

};
