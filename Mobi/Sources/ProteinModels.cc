/*
 @file    ProteinModels.cc
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 1.0
 */

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

// Includes:
#include <iostream>
#include <ProteinModels.h>
#include <String2Number.h>
#include <PdbSaver.h>

using namespace Victor;
using namespace Victor::Mobi;
using namespace Victor::Biopool;

const string PDB = ".pdb";

// CONSTRUCTORS/DESTRUCTOR:

/**
 *   Default Constructor
 * @param nothing.
 */
ProteinModels::ProteinModels() :
		Protein(), verbose(false) {
}

/**
 *  Constructor. Copy another Protein object.
 *@param Protein
 */
ProteinModels::ProteinModels(const Protein& _orig) :
		Protein(_orig), verbose(false) {
}

/**
 *  DESTRUCTOR
 *@param none
 */
ProteinModels::~ProteinModels() {
}

// PREDICATES:

/**
 Overloading method load, this method load all models in Protein Models Object
 @param  (PdbLoader&) , object PdbLoader
 @return void
 */
void ProteinModels::load(PdbLoader& pl) {

	unsigned int d = 0;

	if (!verbose)
		pl.setNoVerbose();
//Load Models
	for (unsigned int i = 1; i <= pl.getMaxModels(); i++) {
		if (verbose)
			cout << ">>>model#" << i << endl;
		pl.setModel(i);
		pl.checkModel();
		this->Protein::load(pl);
	}

	if (this->size() != 0) {
		for (unsigned int i = 0; i < this->size(); i++)
			original_models.push_back(*(this->getSpacer(i)));

		for (unsigned int i = 0; i < this->getSpacer(d)->sizeAmino(); i++)
			sequence.push_back(this->getSpacer(d)->getAmino(i).getType1L());

		if (verbose)
			cout << "\n" << "Models load in protein" << endl;
	}
}

/**
 Internal method
 Chose how much models do you want loa
 @param  (PdbLoader&) , object PdbLoader
 @return unsigned int
 */
unsigned int ProteinModels::selectModels(PdbLoader& pl) {

	unsigned int modelNum;
	string input;

	cout << "How much do you want load? Insert a value between 2 and ";
	cout << pl.getMaxModels() << endl;

	getline(cin, input);

	modelNum = stouiDEF(input);

	while (modelNum > pl.getMaxModels() || modelNum < 2) {
		cerr << "Number of models selected out of bounds! Reinput please..."
				<< endl;
		getline(cin, input);
		modelNum = stouiDEF(input);

	}

	return modelNum;

}

/**
 Load same models
 @param (PdbLoader&) , object PdbLoader
 @return void
 */
void ProteinModels::loadSameModels(PdbLoader& pl) {

	unsigned int d=0;
	unsigned int modelNum = selectModels(pl);

	for (unsigned int i = 1; i <= modelNum; i++) {
		if (verbose)
			cout << "\t>>>model#" << i << endl;
		pl.setModel(i);	//from pdb loader choose model to load in the protein
		pl.checkModel();
		this->Protein::load(pl);          // creates the Protein object
	}

	if (this->size() != 0) {
		for (unsigned int i = 0; i < this->size(); i++)
			original_models.push_back(*(this->getSpacer(i)));

		for (unsigned int i = 0; i < this->getSpacer(d)->sizeAmino(); i++)
			sequence.push_back(this->getSpacer(d)->getAmino(i).getType1L());

		if (verbose)
			cout << "\n" << "Models load in protein" << endl;
	}
}

/**
 Save vector of models in output file
 @param (string output) , name of output file
 @return void
 */
void ProteinModels::printModels(string outputFile) {

	/**stream in output*/
	ofstream fout;
	string outputFile_1;
	PdbSaver ps(fout);
	const string VECTOR = "Vector";

	unsigned int i = 0;

	if (verbose)
		cout << "\n ###Save models on file ###" << endl;
	vector<Spacer>::iterator walk = models.begin();
	while (walk != models.end()) {

		outputFile_1 = outputFile + VECTOR + (itosDEF(i)) + PDB;
		if (verbose)
			cout << outputFile_1 << endl;

		fout.open(outputFile_1.c_str());
		if (!fout) {
			ERROR("Could not open file for writing.", exception);
		} else {

			ps.saveSpacer(*walk);
			walk++;
			i++;
			ps.endFile();
			fout.close();

			if (verbose)
				cout << "file cloased" << endl;
		}
	}

}

/**
 Overloading method save, save models in output file
 @param (string output) , name of output file
 @return void
 */
void ProteinModels::save(string outputFile) {

	/**stream in output*/
	ofstream fout;
	string outputFile_1;
	PdbSaver ps(fout);

	if (verbose)
		cout << "Save models on file..." << endl;
	// Open the proper output stream (file or stdout)
	for (unsigned int i = 0; i < this->size(); i++) {

		outputFile_1 = outputFile + (itosDEF(i)) + PDB;

		if (verbose)
			cout << outputFile_1 << endl;

		fout.open(outputFile_1.c_str());

		if (!fout) {
			ERROR("Could not open file for writing.", exception);
		} else
			ps.saveSpacer(*(this->getSpacer(i)));

		ps.endFile();
		fout.close();

		if (verbose)
			cout << "file closed." << endl;
	}
	if (verbose)
		cout << " " << endl;

}

/**
 Remove file created by save method
 @param (string output) , name of output file
 @return void
 */
void ProteinModels::remove(string outputFile) {

	string outputFile_1;

	if (verbose)
		cout << "Deleting file..." << endl;
	// Open the proper output stream (file or stdout)
	for (unsigned int i = 0; i < this->size(); i++) {

		outputFile_1 = outputFile + (itosDEF(i)) + PDB;
		if (std::remove(outputFile_1.c_str()) != 0)
			ERROR("File don't deleted! ", exception);

	}
	if (verbose)
		cout << " " << endl;

}

/**
 Add Spacer in object ProteinModels
 @param (Spacer&) , object Spacer
 @return void
 */
void ProteinModels::addModels(Spacer& sp) {
	this->models.push_back(sp);

}

/**
 Get Spacer from Protein Models
 @param (unsigned int) , number of Spacer (range 0/number Spacer)
 @return Spacer
 */
Spacer ProteinModels::getModel(unsigned int u) {
	return this->models[u];
}
