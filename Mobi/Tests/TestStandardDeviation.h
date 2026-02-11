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


#ifndef MOBI_TESTS_TESTSTANDARDDEVIATION_H_
#define MOBI_TESTS_TESTSTANDARDDEVIATION_H_


#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include "StandardDeviation.h"

using namespace std;
using namespace Victor::Mobi;


class TestStandardDeviation : public CppUnit::TestFixture {

public:
	TestStandardDeviation() {}
	virtual ~TestStandardDeviation(){}


	 static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite(
				"TestProteinModels");

		suiteOfTests->addTest(
				new CppUnit::TestCaller<TestStandardDeviation>(
						"Test1 - test everage scaled distance",
						&TestStandardDeviation::getEverage));

		suiteOfTests->addTest(
						new CppUnit::TestCaller<TestStandardDeviation>(
								"Test2 - test standard deviation",
								&TestStandardDeviation::getSD));
		suiteOfTests->addTest(
								new CppUnit::TestCaller<TestStandardDeviation>(
										"Test3 - test standard deviation angle PHI",
										&TestStandardDeviation::getanglePHI));
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
	void getEverage() {

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

		// --------------------------------------------------
		// 3. everage scaled distance Test
		// --------------------------------------------------

		StandardDeviation st(test, true);

		vector<double> everage = st.get_everage_distance();
		cout << endl << "Test first 5: " << everage[0] << " " << everage[1]
				<< " " << everage[2] << " " << everage[3] << " " << everage[4]
				<< " " << endl;

	}

	void getSD() {

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
		// --------------------------------------------------
		// 3. StandardDeviation Test
		// --------------------------------------------------
		StandardDeviation st(test, true);

		vector<double> SD = st.get_standard_deviation();
		cout << endl << "Test first 5: " << SD[0] << " " << SD[1] << " "
				<< SD[2] << " " << SD[3] << " " << SD[4] << " " << endl;

	}

	void getanglePHI() {

			// --------------------------------------------------
			// 1. Initialization
			// --------------------------------------------------

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
			// 3. StandardDeviation Test
			// --------------------------------------------------
			StandardDeviation st(test, true);

			vector<double> SD = st.get_StandarDev_angle_PHI();
			cout << endl << "Test first 5: " << SD[0] << " " << SD[1] << " "
					<< SD[2] << " " << SD[3] << " " << SD[4] << " " << endl;

		}

};

#endif /* MOBI_TESTS_TESTSTANDARDDEVIATION_H_ */
