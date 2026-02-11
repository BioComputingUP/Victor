/*
 @file    TestSaver.h
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

#ifndef MOBI_TESTS_TESTSAVER_H_
#define MOBI_TESTS_TESTSAVER_H_

#include <iostream>
#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCase.h>

#include "ProteinModels.h"
#include "StandardDeviation.h"
#include "MobiSaver.h"

using namespace std;
using namespace Victor::Mobi;


class TestSaver : public CppUnit::TestFixture {

public:
	TestSaver() {}
	virtual ~TestSaver(){}


	 static CppUnit::Test *suite() {
		CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite(
				"TestSaver");

		suiteOfTests->addTest(
				new CppUnit::TestCaller<TestSaver>(
						"Test1 - Test mask average filtered with DSSP",
						&TestSaver::testFilterDSSP));

		suiteOfTests->addTest(
						new CppUnit::TestCaller<TestSaver>(
								"Test2 - Test mask scaled distance filtered by logic mask",
								&TestSaver::testFilterMask));
		return suiteOfTests;
	}

	 /** @brief Setup method. */
	void setUp() {
	}
	/** @brief Teardown method. */
	void tearDown() {
	}

protected:
	 /** @brief Test save mobility from average filtered by DSSP */
	void testFilterDSSP() {
		// --------------------------------------------------
		// 1. Initialization
		// --------------------------------------------------
		unsigned int d = 0;
		const string OUTPDB = "ResultMobi.pdb";
		ifstream inputFile("/home/riccardo/mobi/1AB2_input.pdb");
		string outputFile = "../Mobi/data/stdout";
		string out_pdb = outputFile + OUTPDB;
		ofstream ofstream(out_pdb.c_str(), ofstream::app);
		PdbLoader pl(inputFile);
		pl.setVerbose();
		ProteinModels test;
		test.setVerbose();
		// --------------------------------------------------
		// 2. Load
		// --------------------------------------------------
		test.loadSameModels(pl);
		// --------------------------------------------------
		// 3. Save
		// --------------------------------------------------
		test.save(outputFile);
		// --------------------------------------------------
		// 4. Rototraslation
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
		// 5. Calculate average
		// --------------------------------------------------
		MobiSaver saver(test, outputFile, ofstream, true);

		StandardDeviation std(test, true);
		vector<double> ever;
		ever = std.get_everage_distance();

		// --------------------------------------------------
		// 6. save everage
		// --------------------------------------------------
		saver.mob_eveScalD(ever);

		// --------------------------------------------------
		// 7. Calculate secondary structure
		// --------------------------------------------------
		SecondaryStructure sstr(test, true);

		vector<char> MOB;
		MOB = sstr.getMobilitySecondaryStructure();

		// --------------------------------------------------
		// 8. save secondary structure
		// --------------------------------------------------

		saver.mob_SecS(MOB);

		// --------------------------------------------------
		// 9. save mobility filtered
		// --------------------------------------------------

		saver.mob_eveScalD_filtSecS(ever, MOB);

		// --------------------------------------------------
		// 9. remove trash
		// --------------------------------------------------

		test.remove(outputFile);
		remove(out_pdb.c_str());
		string trash = outputFile + "TMScore.pdb_atm";
		remove(trash.c_str());

	}


	 /** @brief Test save mobility from standard deviation filtered whith logic mask */
	void testFilterMask() {
			// --------------------------------------------------
			// 1. Initialization
			// --------------------------------------------------
			unsigned int d = 0;
			const string OUTPDB = "ResultMobi.pdb";
			ifstream inputFile("/home/riccardo/mobi/1AB2_input.pdb");
			string outputFile = "../Mobi/data/stdout";
			string out_pdb = outputFile + OUTPDB;
			ofstream ofstream(out_pdb.c_str(), ofstream::app);
			PdbLoader pl(inputFile);
			pl.setVerbose();
			ProteinModels test;
			test.setVerbose();
			// --------------------------------------------------
			// 2. Load
			// --------------------------------------------------
			test.loadSameModels(pl);
			// --------------------------------------------------
			// 3. Save
			// --------------------------------------------------
			test.save(outputFile);
			// --------------------------------------------------
			// 4. Rototraslation
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
			// 5. Calculate average
			// --------------------------------------------------
			MobiSaver saver(test, outputFile, ofstream, true);

			StandardDeviation std(test, true);

			vector<double> SD;
			SD = std.get_standard_deviation();


			// --------------------------------------------------
			// 6. save standard deviation
			// --------------------------------------------------
			saver.mob_stanD(SD);

		// --------------------------------------------------
		// 7. save mobility filtered with logic mask
		//	0 = fixed
		//  1 = mobile
		//
		//			1011   »   1111
		//			1101   »   1111
		//			10011  »  11111
		//			11001  »  11111
		//			01010  »  00000
		//			00100  »  00000
		//			001100 » 000000
		//
		// ---------------------------------------------------

			saver.mob_stanD_withMask(SD);

			// --------------------------------------------------
			// 9. remove trash
			// --------------------------------------------------

			test.remove(outputFile);
			remove(out_pdb.c_str());
			string trash = outputFile + "TMScore.pdb_atm";
			remove(trash.c_str());

		}



};



#endif /* MOBI_TESTS_TESTSAVER_H_ */
