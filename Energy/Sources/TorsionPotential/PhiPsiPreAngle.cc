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
#include <PhiPsiPreAngle.h>
using namespace Biopool;
// Global constants, typedefs, etc. (to avoid):
// CONSTRUCTORS/DESTRUCTOR:

/**
 * @Description Constructor that sets the arc steps for the table's size
 * @param arc values for table 1 and 2(int ,int), file name where torsion angles are, usually data/tor.par(string)
 */
PhiPsiPreAngle::PhiPsiPreAngle(int SET_ARC1, int SET_ARC2, string knownledge) :
TOR_PARAM_FILE(knownledge), ARC_STEP(SET_ARC1), ARC_STEP2(SET_ARC2),
propensities(), all_propensities() {
    SIZE_OF_TABLE = 360 / ARC_STEP;
    SIZE_OF_TABLE2 = 360 / ARC_STEP2;
    pConstructData();
}

/**
 * @Description  Method that using the torsion file, creates the structure
 * @param none
 * @return changes are made internally(void)
 */
void PhiPsiPreAngle::pConstructData() {
    char *victor = getenv("VICTOR_ROOT");
    if (victor == NULL)
        ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
    string inputFile;
    if (TOR_PARAM_FILE == "data/tor.par") {
        inputFile = getenv("VICTOR_ROOT");
        inputFile += TOR_PARAM_FILE;
    } else {
        string path = "data/tor.par";
        inputFile += path;
    }
    if (inputFile.length() < 3)
        ERROR("Environment variable VICTOR_ROOT was not found.", exception);

    ifstream input(inputFile.c_str());

    if (!input)
        ERROR("Could not read data file. " + inputFile, exception);


    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        vector<vector<vector<vector<int>* >* >* >* tmpB =
                new vector<vector<vector<vector<int>* >* >* >;
        (*tmpB).reserve(SIZE_OF_TABLE);

        for (int i = 0; i < SIZE_OF_TABLE; i++) {
            vector<vector<vector<int>* >* >* tmpA =
                    new vector<vector<vector<int>* >* >;
            (*tmpA).reserve(SIZE_OF_TABLE);

            for (int i = 0; i < SIZE_OF_TABLE; i++) {
                vector<vector<int>* >* tmpC =
                        new vector<vector<int>* >;
                (*tmpC).reserve(SIZE_OF_TABLE2);

                for (int i = 0; i < SIZE_OF_TABLE2; i++) {
                    vector<int>* tmpD = new vector<int>;
                    (*tmpD).reserve(SIZE_OF_TABLE2);

                    for (int i = 0; i < SIZE_OF_TABLE2; i++) {
                        (*tmpD).push_back(0);
                    }
                    (*tmpC).push_back(tmpD);
                }
                (*tmpA).push_back(tmpC);
            }
            (*tmpB).push_back(tmpA);
        }
        (propensities).push_back(tmpB);
    }


    long numCount;
    double phi, psi, prephi, prepsi;
    string name;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_count[i] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i += 1)
        for (int j = 0; j < SIZE_OF_TABLE; j++)
            for (int k = 0; k < SIZE_OF_TABLE; k++)
                for (int l = 0; l < SIZE_OF_TABLE2; l++)
                    for (int m = 0; m < SIZE_OF_TABLE2; m++)
                        (*(*(*(*propensities[i])[j])[k])[l])[m] = 1;

    input >> numCount;

    for (int i = 0; i < numCount; i++) {
        input >> phi >> psi >> name >> prephi >> prepsi;
        skipToNewLine(input);
        sAddProp(aminoAcidThreeLetterTranslator(name), sGetPropBin(phi),
                sGetPropBin(psi), sGetPropBin2(prephi), sGetPropBin2(prepsi));

    }
    input.close();
    // calculate some constant values:
    total = 1;

    for (int i = 0; i < SIZE_OF_TABLE; i++) {
        vector<vector<vector<int>* >* >* tmpE =
                new vector<vector<vector<int>* >* >;
        (*tmpE).reserve(SIZE_OF_TABLE);

        for (int i = 0; i < SIZE_OF_TABLE; i++) {
            vector<vector<int>* >* tmpF =
                    new vector<vector<int>* >;
            (*tmpF).reserve(SIZE_OF_TABLE2);

            for (int i = 0; i < SIZE_OF_TABLE2; i++) {
                vector<int>* tmpG = new vector<int>;
                (*tmpG).reserve(SIZE_OF_TABLE2);

                for (int i = 0; i < SIZE_OF_TABLE2; i++) {
                    (*tmpG).push_back(0);
                }
                (*tmpF).push_back(tmpG);
            }
            (*tmpE).push_back(tmpF);
        }
        (all_propensities).push_back(tmpE);
    }


    for (int j = 0; j < SIZE_OF_TABLE; j++)
        for (int k = 0; k < SIZE_OF_TABLE; k++)
            for (int l = 0; l < SIZE_OF_TABLE2; l++)
                for (int m = 0; m < SIZE_OF_TABLE2; m++)
                    (*(*(*all_propensities[j])[k])[l])[m] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        total += amino_count[i];

        for (int j = 0; j < SIZE_OF_TABLE; j++)
            for (int k = 0; k < SIZE_OF_TABLE; k++)
                for (int l = 0; l < SIZE_OF_TABLE2; l++)
                    for (int m = 0; m < SIZE_OF_TABLE2; m++)
                        (*(*(*all_propensities[j])[k])[l])[m] += (*(*(*(*propensities[i])[j])[k])[l])[m];
    }
    //Construct the max propensities vector
    pConstructMaxPropensities();
}

/**
 * @Description Resets all the propensities set previously
 * @param none
 * @return changes are made internally(void)
 */
void PhiPsiPreAngle::pResetData() {
    for (unsigned int i = 0; i < propensities.size(); i++) {
        for (unsigned int j = 0; j < (*propensities[i]).size(); j++) {
            for (unsigned int k = 0; k < (*(*propensities[i])[j]).size(); k++) {
                for (unsigned int l = 0; l < (*(*(*propensities[i])[j])[k]).size(); l++) {
                    (*(*(*(*propensities[i])[j])[k])[l]).clear();
                    delete ((*(*(*propensities[i])[j])[k])[l]);
                }
                (*(*(*propensities[i])[j])[k]).clear();
                delete ((*(*propensities[i])[j])[k]);
            }
            (*(*propensities[i])[j]).clear();
            delete ((*propensities[i])[j]);
        }
        (*propensities[i]).clear();
        delete (propensities[i]);
    }
    propensities.clear();


    for (unsigned int i = 0; i < all_propensities.size(); i++) {
        for (unsigned j = 0; j < (*all_propensities[i]).size(); j++) {
            for (unsigned int k = 0; k < (*(*all_propensities[i])[j]).size(); k++) {
                (*(*(*all_propensities[i])[j])[k]).clear();
                delete ((*(*all_propensities[i])[j])[k]);
            }
            (*(*all_propensities[i])[j]).clear();
            delete ((*all_propensities[i])[j]);
        }
        (*all_propensities[i]).clear();
        delete (all_propensities[i]);
    }
    all_propensities.clear();
    for (unsigned int i = 0; i < amino_max_propensities_pre_angle.size(); i++) {
        for (unsigned int j = 0; j < (*amino_max_propensities_pre_angle[i]).size(); j++) {
            (*(*amino_max_propensities_pre_angle[i])[j]).clear();
            delete ((*amino_max_propensities_pre_angle[i])[j]);
        }
        (*amino_max_propensities_pre_angle[i]).clear();
        delete (amino_max_propensities_pre_angle[i]);
    }
    amino_max_propensities_pre_angle.clear();
}


// PREDICATES:

/**
 * @Description calculates energy for a protein 
 * @param reference to the spacer containing the protein (Spacer&)
 * @return energy value for the set of amino acids in the spacer (long double)
 */
long double PhiPsiPreAngle::calculateEnergy(Spacer& sp) {
    long double en = 0.0;

    for (unsigned int i = 2; i < sp.sizeAmino() - 1; i++) {
        en += calculateEnergy(sp.getAmino(i));
    }
    return en;
}

/**
 * @Description calculates energy for part of the protein
 * @param reference to the spacer that contains the protein(Spacer&) , start and end position of the protein part (unsigned int , unsigned int)
 * @return energy value for the set of amino acids in the selected section (long double)
 */
long double PhiPsiPreAngle::calculateEnergy(Spacer& sp, unsigned int index1,
        unsigned int index2) {
    long double en = 0.0;

    index1 = (index1 >= 2) ? index1 : 2;
    index2 = (index2 <= sp.sizeAmino() - 1) ? index2 : sp.sizeAmino() - 1;

    for (unsigned int i = index1; i < index2; i++)
        en += calculateEnergy(sp.getAmino(i));

    return en;
}

/**
 * @Description calculates energy for a specific amino acid
 * @param aminoacid reference (AminoAcid&)
 * @return energy value for given amino acid  (long double)
 */
long double PhiPsiPreAngle::calculateEnergy(AminoAcid& aa) {
    // current amino acid type:
    int table_entry = aminoAcidThreeLetterTranslator(aa.getType());

    int n = aa.sizeInBonds();
    if (n == 0)
        ERROR("Invalid Aminoacid, impossible to calculate pre-angle.", exception);
    // offsets for the array

    int x = sGetPropBin(aa.getPhi(true));
    int y = sGetPropBin(aa.getPsi(true));
    int z = sGetPropBin2(aa.getInBond(0).getPhi());
    int l = sGetPropBin2(aa.getInBond(0).getPsi());

    return -log((static_cast<double> ((*(*(*(*propensities[table_entry])[x])[y])[z])[l])
            / (*(*(*all_propensities[x])[y])[z])[l])
            / (static_cast<double> (amino_count[table_entry])
            / total));
}

/**
 * @Description Adds propensity value for a amino acid type
 * @param corresponding amino acid code(int), corresponding propensity values(int, int, int, int))
 * @return changes are made internally(void)
 */
inline void PhiPsiPreAngle::sAddProp(int code, int x, int y, int z, int l) {
    (*(*(*(*propensities[code])[x])[y])[z])[l]++;
    amino_count[code]++;
}

/**
 * @Description Sets the arc step, needed to define the Size table
 * @param quantity of steps(int)
 * @return changes are made internally(void)
 */
inline void PhiPsiPreAngle::setArcStep(int n) {
    ARC_STEP = n;
    SIZE_OF_TABLE = 360 / ARC_STEP;
    pResetData();
    pConstructData();
}

/**
 * @Description Returns the propensity binding value using table 1(granularity) for a specific angle
 * @param angle in degrees
 * @return corresponding prop value (int)
 */
inline int PhiPsiPreAngle::sGetPropBin(double p) {
    int x = (int) (p + 180) / ARC_STEP;
    if (x == SIZE_OF_TABLE)
        x -= 1; // x is exactly 180 degrees and boosts array
    return x;
}

/**
 * @Description Returns the propensity binding using table 2(granularity for i-1 aa) value for a specific angle
 * @param angle in degrees(double)
 * @return corresponding prop value (int)
 */
inline int PhiPsiPreAngle::sGetPropBin2(double p) {
    int x = (int) (p + 180) / ARC_STEP2;
    if (x == SIZE_OF_TABLE2)
        x -= 1; // x is exactly 180 degrees and boosts array
    return x;
}

/**
 * @Description obtains the maximum propensity value for an amino acid type considering all the granularities
 * @param amino acid code(int)
 * @return corresponding propensity value(double)
 */
double PhiPsiPreAngle::pGetMaxPropensities(int amino) {
    double max = 0.00;

    for (int j = 0; j < SIZE_OF_TABLE; j++) {
        for (int k = 0; k < SIZE_OF_TABLE; k++) {
            for (int r = 0; r < SIZE_OF_TABLE2; r++) {
                for (int s = 0; s < SIZE_OF_TABLE2; s++) {
                    int tmp = (*(*(*all_propensities[j])[k])[r])[s];
                    int tmp2 = (*(*(*(*propensities[amino])[j])[k])[r])[s];
                    double propensities = ((static_cast<double> (tmp2) / tmp)
                            / ((static_cast<double> (amino_count[amino])) / total));
                    if (propensities > max) {
                        max = propensities;
                    }
                }
            }
        }
    }

    return max;
}

/**
 * @Description obtains the maximum propensity value for an amino acid type considering a given phi and psi values
 * considering all the granularities
 * @param amino acid code(int), values for the previous phi and psi angles(int , int)
 * @return corresponding propensity value(double)
 */
double PhiPsiPreAngle::pGetMaxPropensities(int amino, int prephi, int prepsi) {
    double max = 0.00;

    for (int j = 0; j < SIZE_OF_TABLE; j++) {
        for (int k = 0; k < SIZE_OF_TABLE; k++) {
            int tmp = (*(*(*all_propensities[j])[k])[prephi])[prepsi];
            int tmp2 = (*(*(*(*propensities[amino])[j])[k])[prephi])[prepsi];
            double propensities = ((static_cast<double> (tmp2) / tmp)
                    / ((static_cast<double> (amino_count[amino])) / total));
            if (propensities > max) {
                max = propensities;
            }
        }
    }
    return max;
}

/**
 * @Description obtains the maximum propensity value for an amino acid type  
 * @param amino acid code(int)
 * @return corresponding propensity value(double)
 */
double PhiPsiPreAngle::pReturnMaxPropensities(int amino) {
    return amino_max_propensities[amino];
}

/**
 * @Description obtains the maximum propensity value for an amino acid type considering a given phi and psi values
 * @param amino acid code(int), values for the previous phi and psi angles(int , int)
 * @return corresponding propensity value(double)
 */
double PhiPsiPreAngle::pReturnMaxPropensitiesPreAngle(int amino, int prephi, int prepsi) {
    return (*(*amino_max_propensities_pre_angle[amino])[prephi])[prepsi];
}

/**
 * @Description Creates the maximum propensities based on the max propensities based on knowledge
 * @param none
 * @return changes are made internally(void)
 */
void PhiPsiPreAngle::pConstructMaxPropensities() {
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_max_propensities.push_back(0);
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        amino_max_propensities[i] = pGetMaxPropensities(i);
    }
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        vector<vector<double>* >* tmpE = new vector<vector<double>* >;
        (*tmpE).reserve(SIZE_OF_TABLE2);
        for (int i = 0; i < SIZE_OF_TABLE2; i++) {
            vector<double>* tmpD = new vector<double>;
            (*tmpD).reserve(SIZE_OF_TABLE2);
            for (int i = 0; i < SIZE_OF_TABLE2; i++) {
                (*tmpD).push_back(0);
            }
            (*tmpE).push_back(tmpD);
        }
        amino_max_propensities_pre_angle.push_back(tmpE);
    }

    for (int j = 0; j < AminoAcid_CODE_SIZE; j++)
        for (int k = 0; k < SIZE_OF_TABLE2; k++)
            for (int l = 0; l < SIZE_OF_TABLE2; l++)
                (*(*amino_max_propensities_pre_angle[j])[k])[l] = pGetMaxPropensities(j, k, l);

}



