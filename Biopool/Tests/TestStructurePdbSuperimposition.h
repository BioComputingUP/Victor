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
#include <StructurePdbSuperimposition.h>
#include <PdbLoader.h>

using namespace std;
using namespace Victor::Biopool;

class TestStructurePdbSuperimposition : public CppUnit::TestFixture {

public:
    void setUp() {}

    void tearDown() {}

	static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestStructurePdbSuperimposition");

		suiteOfTests->addTest(new CppUnit::TestCaller<TestStructurePdbSuperimposition>("testSuperimposition", &TestStructurePdbSuperimposition::testSuperimposition));
		suiteOfTests->addTest(new CppUnit::TestCaller<TestStructurePdbSuperimposition>("testKabsch", &TestStructurePdbSuperimposition::testKabsch));

		return suiteOfTests;
	}

	void testSuperimposition() {
		std::vector<string> inputFile;
		string path = getenv("VICTOR_ROOT");
		inputFile[0] = path + "Biopool/Tests/data/15C8_H_input.pdb";
		inputFile[1] = path + "Biopool/Tests/data/25C8_H_input.pdb";

		ifstream inFile1(inputFile[0].c_str());
		ifstream inFile2(inputFile[1].c_str());
		
		PdbLoader pl1(inFile1);
		PdbLoader pl2(inFile2);
		
		pl1.setNoHAtoms();
		pl2.setNoHAtoms();

		Protein prot1, prot2;
		prot1.load(pl1);
		prot2.load(pl2);
		
		vector<char> chainsIdProt1, chainsIdProt2;
		chainsIdProt1 = prot1.getAllChains();
		chainsIdProt2 = prot2.getAllChains();
		
		Spacer *sp1, *sp2;
		sp1 = prot1.getSpacer(chainsIdProt1[0]);
		sp2 = prot2.getSpacer(chainsIdProt2[0]);
		
		Eigen::MatrixX3d matrix1, matrix2;
		matrix1 = StructurePdbSuperimposition::fromSpacerToEigenMatrix(sp1);
		matrix2 = StructurePdbSuperimposition::fromSpacerToEigenMatrix(sp2);

		double result[4];

		double rmsd = StructurePdbSuperimposition::RMSD(matrix1, matrix2, *sp1, *sp2);
		double maxSub = StructurePdbSuperimposition::MaxSubScore(matrix1, matrix2, *sp1, *sp2, 4, 3.5, true);
		double gdtTs = StructurePdbSuperimposition::GDTTS(matrix1, matrix2, *sp1, *sp2, 4, result);
		double tm = StructurePdbSuperimposition::TMScore(matrix1, matrix2, *sp1, *sp2, 4);

        CPPUNIT_ASSERT((-rmsd - 1.5675) < EPSILON && (rmsd - 1.5675) < EPSILON);
        CPPUNIT_ASSERT((-maxSub - 0.857753) < EPSILON && (maxSub - 0.857753) < EPSILON);
        CPPUNIT_ASSERT((-gdtTs - 0.770826) < EPSILON && (gdtTs - 0.770826) < EPSILON);
        CPPUNIT_ASSERT((-tm - 0.932606) < EPSILON && (tm - 0.932606) < EPSILON);
	}

	void testKabsch() {
		Eigen::MatrixX3d m1, m2;
		m1.resize(2,3);
		m2.resize(2,3);
		
		m1 << 0.6, 0.5, 0.8, -0.2, 0.5, -0.6;
		m2 << -0.3, -0.4, -0.2, 0.5, 0.1, 0.2;

		std::pair<Eigen::Matrix3d,Eigen::Vector3d> Rt;

		Rt = StructurePdbSuperimposition::KabschAlgorithm(m1, m2);  

		CPPUNIT_ASSERT(Rt.first(0,0) == -0.571275 && Rt.first(1,1) == -0.361761 && Rt.first(2,2) == -0.0564811 && Rt.second(0) == -0.0224084 && Rt.second(1) == 0.0231343 && Rt.second(2) == 0.505016);
	}	
};
