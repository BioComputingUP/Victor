/*
 @file    TestProteinModels.h
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

#ifndef MOBI_TESTS_TESTPROTEINMODELS_H_
#define MOBI_TESTS_TESTPROTEINMODELS_H_

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include "ProteinModels.h"
#include <TmScore.h>

using namespace std;
using namespace Victor::Mobi;


class TestProteinModels : public CppUnit::TestFixture {

public:
	TestProteinModels() {}
	virtual ~TestProteinModels(){}


	 static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite(
				"TestProteinModels");

		suiteOfTests->addTest(
				new CppUnit::TestCaller<TestProteinModels>(
						"Test1 - Populate the ProteinModels",
						&TestProteinModels::testPopulation));
		suiteOfTests->addTest(
						new CppUnit::TestCaller<TestProteinModels>(
								"Test2 - Correct aminoacid sequence",
								&TestProteinModels::aminoSequence));
		suiteOfTests->addTest(
								new CppUnit::TestCaller<TestProteinModels>(
										"Test3 - Correct load protein rototralated",
										&TestProteinModels::rototraslated));
		suiteOfTests->addTest(
						new CppUnit::TestCaller<TestProteinModels>(
								"Test4 - Save protein models rototraslated",
								&TestProteinModels::printModels));

		return suiteOfTests;
	}

	 /** @brief Setup method. */
	void setUp() {
	}
	/** @brief Teardown method. */
	void tearDown() {
	}

protected:
	 /** @brief Test load */
	void testPopulation() {

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
		test.load(pl);

		// --------------------------------------------------
		// 3. Test
		// --------------------------------------------------
		cout << endl << "Models in pdb file is: " << pl.getMaxModels() << "\n"
				<< "Models load is" << test.original_models.size() << endl;
		CPPUNIT_ASSERT(pl.getMaxModels() == test.original_models.size());

	}
	 /** @brief Test aminoSequence */
	void aminoSequence() {

		// --------------------------------------------------
		// 1. Initialization
		// --------------------------------------------------
		unsigned d = 0;
		string aminoaced_sequence =
				"GSGNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSISLRYEGRVYHYRINTASDGKLYVSSESRFNTLAELVHHHSTVADGLITTLHYPAPKRGIHRD";
		string amino;
		ifstream inputFile("/home/riccardo/mobi/1AB2_input.pdb");
		PdbLoader pl(inputFile);
		pl.setVerbose();
		ProteinModels* test = new ProteinModels;
		test->setVerbose();

		// --------------------------------------------------
		// 2. Load
		// --------------------------------------------------
		pl.setModel(1);
		pl.checkModel();
		test->Protein::load(pl);

		for (unsigned int i = 0; i < test->getSpacer(d)->sizeAmino(); i++)
			amino += (test->getSpacer(d)->getAmino(i).getType1L());

		// --------------------------------------------------
		// 3. Test
		// --------------------------------------------------
		cout << endl << "Sequence in pdb file is: " << aminoaced_sequence
				<< "\n" << "Sequence load is" << amino << endl;
		CPPUNIT_ASSERT(aminoaced_sequence == amino);

	}
	/** @brief Test function of TmScore */
	void rototraslated() {

		// --------------------------------------------------
		// 1. Initialization
		// --------------------------------------------------
		unsigned int d = 0;

		ifstream inputFile("/home/riccardo/mobi/1AB2_input.pdb");
		string outputFile = "../Mobi/data/stdout";
		PdbLoader pl(inputFile);
		pl.setVerbose();
		ProteinModels test;
		test.setVerbose();

		// --------------------------------------------------
		// 2. Load
		// --------------------------------------------------
		test.loadSameModels(pl);

		// --------------------------------------------------
		// 2. Save
		// --------------------------------------------------

		test.save(outputFile);

		// --------------------------------------------------
		// 2. Rototranslation
		// --------------------------------------------------
		TmScore tm("../Mobi/data/TMscore", outputFile, true);

		Protein* traslata = new Protein();

		for (unsigned int i = 0; i < test.size() - 1; i++)
			for (unsigned int j = i + 1; j < test.size(); j++) {
				traslata = tm.TmImpose(outputFile + (itosDEF(i)),
						outputFile + (itosDEF(j)));
				if (traslata != NULL) {

					test.addModels(*(traslata->getSpacer(d)));
				} else
					ERROR("Error in the creation of shift models.", exeption);

				test.addModels(*(test.getSpacer(j)));

			}
		// --------------------------------------------------------
		// 2. Test save model vector that contain couple of models,
		//    the first rototranslate, the second is the model
		// --------------------------------------------------------
		cout << endl << "Models load is: " << test.original_models.size()
				<< "\n"
				<< "Models load whit models rototranslated (all models imposed on all models) "
				<< test.models.size() << endl;
		unsigned int test_model_size = (test.original_models.size()
				* (test.original_models.size() - 1));
		CPPUNIT_ASSERT(test_model_size == test.models.size());

	}
	/** @brief Test function printmodels */
	void printModels() {

		// --------------------------------------------------
		// 1. Initialization
		// --------------------------------------------------
		unsigned int d = 0;

		ifstream inputFile("/home/riccardo/mobi/1AB2_input.pdb");
		string outputFile = "../Mobi/data/stdout";
		PdbLoader pl(inputFile);
		pl.setVerbose();
		ProteinModels test;
		test.setVerbose();
		// --------------------------------------------------
		// 2. Load
		// --------------------------------------------------
		test.loadSameModels(pl);
		// --------------------------------------------------
		// 2. Save
		// --------------------------------------------------
		test.save(outputFile);
		// --------------------------------------------------
		// 2. Rototraslation
		// --------------------------------------------------
		TmScore tm("../Mobi/data/TMscore", outputFile, true);

		Protein* traslata = new Protein();

		for (unsigned int i = 0; i < test.size() - 1; i++)
			for (unsigned int j = i + 1; j < test.size(); j++) {
				traslata = tm.TmImpose(outputFile + (itosDEF(i)),
						outputFile + (itosDEF(j)));
				if (traslata != NULL) {

					test.addModels(*(traslata->getSpacer(d)));
				} else
					ERROR("Error in the creation of shift models.", exeption);

				test.addModels(*(test.getSpacer(j)));

			}
		// ------------------------------------------------------
		// 5. save model vector that contain couple of models,
		//    the first rototranslate, the second is the model
		// ------------------------------------------------------

		test.printModels(outputFile);
		cout << endl << "File saved in folder data: " << endl;

	}

};

#endif /* MOBI_TESTS_TESTPROTEINMODELS_H_ */
