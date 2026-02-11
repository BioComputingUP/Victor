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
#include "SecondaryStructure.h"

using namespace Victor;
using namespace Victor::Mobi;
using namespace Victor::Biopool;

// CONSTRUCTORS/DESTRUCTOR:

/**
 *   Default Constructor
 * 	 @param ( ProteinModels& ), object PreteinModels .
 */
SecondaryStructure::SecondaryStructure(const ProteinModels& modelli,
		bool verbose) :
		verbose(verbose), models(modelli.original_models) {
}

/**
 *   Default Deconstructor
 * 	 @param ( ProteinModels& ), object PreteinModels .
 */
SecondaryStructure::~SecondaryStructure() {
}

/**
 * get the secondary structure. Return is a vector of vector of char where vector[i][j]
 * 'i' is a specific amino acid
 * 'j' is a specific model
 * and char is a letter of secondary structure of amino acid 'i' and models 'j'
 * @return vector< vector<char> >  , return vector of vector of char.
 */
vector< vector<char> > SecondaryStructure::getSecStructFromModels(){

	unsigned int num_amino = 0;
		if (models.size() != 0)
			num_amino = models[0].sizeAmino();
		else
			ERROR("Nothing models in the protein, remember load models.", exception);

		//create a support vector to get secondary structure
		vector<set<char> > vettore_di_supporto(num_amino, set<char>());

		//create a vector for manipulated all amiocid's structure
		vector< vector<char> > sec_Structures(num_amino,
				vector<char>(models.size()));

		if (verbose) {
			cout << "\nStart procedure select secondary structure... " << endl;
			cout << "Model size: " << models.size() << endl;
		}
		//for each models
		for (unsigned int i = 0; i < models.size(); i++) {

			models[i].setDSSP(false);
			vettore_di_supporto = models[i].getDSSP();
			//cout << "support vector size: " << vettore_di_supporto.size() << endl;

			//iteration on vector of set
			for (unsigned int j = 0; j < num_amino; j++) {
				//extract first char value and push it in vector char
				sec_Structures[j][i] = (*(vettore_di_supporto[j].begin()));

			}

		}

return sec_Structures;

}

/**
 * get the mobility analyzing secondary structures
 * @return vector<char> , return vector of mobility from Secondary Structure
 */
vector<char> SecondaryStructure::getMobilitySecondaryStructure() {

	mobility.clear();
	vector< vector<char> > sec_Structures;

	unsigned int num_amino = 0;
	if (models.size() != 0)
		num_amino = models[0].sizeAmino();
	else
		ERROR("Nothing models in the protein, remember load models.", exception);


	sec_Structures = getSecStructFromModels();

	char test;
	bool flag;

	for (unsigned int i = 0; i < num_amino; i++) {

		test = sec_Structures[i][0];
		flag = true;

		//if the first letter is empty or is a S
		if ((test != 'G') && (test != 'H') && (test != 'I') && (test != 'T')
				&& (test != 'E') && (test != 'B')) {
			//verifing that is all S and/or empty space.

			for (unsigned int j = 0; j < models.size(); j++) {

				if ((sec_Structures[i][j] != 'G')
						&& (sec_Structures[i][j] != 'H')
						&& (sec_Structures[i][j] != 'I')
						&& (sec_Structures[i][j] != 'T')
						&& (sec_Structures[i][j] != 'E')
						&& (sec_Structures[i][j] != 'B')) {
				} else {

					mobility.push_back('M');
					flag = false;
					break;

				}

			}
			if (flag) {
				mobility.push_back('c');
			}
			//the first letter is the letter of Secondary Structure
		} else {

			//verifing that all letter is the same
			for (unsigned int j = 0; j < models.size(); j++) {

				if ((sec_Structures[i][j] == test)) {
				} else {

					mobility.push_back('M');
					flag = false;
					break;

				}

			}
			if (flag) {
				mobility.push_back('.');
			}

		}

	}

	return mobility;
}

