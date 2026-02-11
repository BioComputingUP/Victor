/*
 @file    TestSecondaryStructure.h
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 1.0
 */

/* This file is part of Victor.
 Victor is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 Victor is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with Victor. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef MOBI_TESTS_TESTSECONDARYSTRUCTURE_H_
#define MOBI_TESTS_TESTSECONDARYSTRUCTURE_H_


#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include "ProteinModels.h"
#include "SecondaryStructure.h"

using namespace std;
using namespace Victor::Mobi;


class TestSecondaryStructure : public CppUnit::TestFixture {

public:
	TestSecondaryStructure() {}
	virtual ~TestSecondaryStructure(){}


	 static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite(
				"TestSecondaryStructure");

		suiteOfTests->addTest(
				new CppUnit::TestCaller<TestSecondaryStructure>(
						"Test1 - Test Secondary Structure",
						&TestSecondaryStructure::testStructure));

		return suiteOfTests;
	}

	 /** @brief Setup method. */
	void setUp() {
	}
	/** @brief Teardown method. */
	void tearDown() {
	}

protected:
	 /** @brief Test Secondary Structure */
	void testStructure() {

		// --------------------------------------------------
		// 1. Initialization
		// --------------------------------------------------

		ifstream inputFile("/home/riccardo/mobi/1AB2_input.pdb");
		PdbLoader pl(inputFile);
		pl.setVerbose();
		ProteinModels test;
		test.setVerbose();

		// --------------------------------------------------
		// 2. Load
		// --------------------------------------------------
		test.loadSameModels(pl);

		// --------------------------------------------------
		// 3. Secondary Structure
		// --------------------------------------------------

		SecondaryStructure sstr(test, true);

		vector< vector<char> > sec_structures;
		sec_structures = sstr.getSecStructFromModels();

		vector<char> MOB;
			MOB = sstr.getMobilitySecondaryStructure();


		// ----------------------------------------------------
		// 4. Print vector o secondary structure for amino acid
		// ----------------------------------------------------
			cout << "\nSECONDARY STRUCTURE" << endl;
			for (unsigned int i=0; i<109 ; i++){
				cout << "Print secondary structure amino acid number #" << i << endl;
				for (unsigned int j = 0; j < test.original_models.size(); j++)
					cout << sec_structures[i][j] << endl;

			cout << "Final valutation of amino acid #" << i << " is " << MOB[i] <<endl;

			}
			cout << "\n" << endl;

	}

};
#endif /* MOBI_TESTS_TESTSECONDARYSTRUCTURE_H_ */
