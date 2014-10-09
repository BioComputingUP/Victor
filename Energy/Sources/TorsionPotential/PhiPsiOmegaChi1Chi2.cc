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
#include <PhiPsiOmegaChi1Chi2.h>

using namespace Victor;

// Global constants, typedefs, etc. (to avoid):

// CONSTRUCTORS/DESTRUCTOR:

/**
 * @Description  Basic constructor
 */
PhiPsiOmegaChi1Chi2::PhiPsiOmegaChi1Chi2(int SET_ARC1, string knownledge) :
TOR_PARAM_FILE(knownledge), ARC_STEP(SET_ARC1),
propensities(), all_propensities() {
    SIZE_OF_TABLE = 360 / ARC_STEP;
    CHI_RANGE = 8;
    RANGE_OMEGA = 3;
    pConstructData();
}

/**
 * @Description  Sets all the information need for the class
 * @param   none
 * @return    changes the object internally (void)   
 */
void PhiPsiOmegaChi1Chi2::pConstructData() {
    char *victor = getenv("VICTOR_ROOT");
    if (victor == NULL)
        ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
    string inputFile;
    if (TOR_PARAM_FILE == "data/tor.par") {
        inputFile = getenv("VICTOR_ROOT");
        inputFile += TOR_PARAM_FILE;
    } else {
        string path = "data/tor.par";
        inputFile = getenv("VICTOR_ROOT");
        inputFile += path;
    }

    if (inputFile.length() < 3)
        ERROR("Environment variable VICTOR_ROOT was not found.", exception);

    ifstream input(inputFile.c_str());

    if (!input)
        ERROR("Could not read data file. " + inputFile, exception);
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        vector<vector<vector<vector<vector<int>* >* >* >* >* tmpC =
                new vector<vector<vector<vector<vector<int>* >* >* >* >;
        (*tmpC).reserve(SIZE_OF_TABLE);

        for (int j = 0; j < SIZE_OF_TABLE; j++) {
            vector<vector<vector<vector<int>* >* >* >* tmpB =
                    new vector<vector<vector<vector<int>* >* >* >;
            (*tmpB).reserve(SIZE_OF_TABLE);

            for (int z = 0; z < SIZE_OF_TABLE; z++) {
                vector<vector<vector<int>* >* >* tmpA =
                        new vector<vector<vector<int>* >* >;
                (*tmpA).reserve(CHI_RANGE);

                for (int x = 0; x < CHI_RANGE; x++) {
                    vector<vector<int>* >* tmpAA =
                            new vector<vector<int>* >;
                    (*tmpAA).reserve(CHI_RANGE);
                    for (int y = 0; y < CHI_RANGE; y++) {
                        vector<int>* tmpZ = new vector<int>;
                        (*tmpZ).reserve(RANGE_OMEGA);

                        for (int n = 0; n < RANGE_OMEGA; n++) {
                            (*tmpZ).push_back(0);
                        }
                        (*tmpAA).push_back(tmpZ);
                    }
                    (*tmpA).push_back(tmpAA);
                }
                (*tmpB).push_back(tmpA);
            }
            (*tmpC).push_back(tmpB);
        }
        (propensities).push_back(tmpC);
    }

    long numCount;
    double phi, psi, chi1, chi2, omega, numchi;
    double vartmp1, vartmp2; //variable for not include pre angle of TOR_PARAM_FILE 
    string name;
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_count[i] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i += 1)
        for (int j = 0; j < SIZE_OF_TABLE; j++)
            for (int k = 0; k < SIZE_OF_TABLE; k++)
                for (int l = 0; l < CHI_RANGE; l++)
                    for (int n = 0; n < CHI_RANGE; n++)
                        for (int m = 0; m < RANGE_OMEGA; m++)
                            (*(*(*(*(*propensities[i])[j])[k])[l])[n])[m] = 1;

    input >> numCount;
    for (int i = 0; i < numCount; i++) {
        input >> phi >> psi >> name >> vartmp1 >> vartmp2 >> omega >> numchi;

        if (numchi == 0) {
            chi1 = 0.00;
            chi2 = 0.00;
            skipToNewLine(input);
        } else if (numchi == 1) {
            input >> chi1;
            chi2 = 0.00;
            skipToNewLine(input);
        } else if (numchi >= 2) {
            input >> chi1;
            input >> chi2;
            skipToNewLine(input);
        }

        sAddProp(aminoAcidThreeLetterTranslator(name), sGetPropBin(phi),
                sGetPropBin(psi), sGetPropChiBin(chi1), sGetPropChiBin(chi2), sGetPropOmegaBin(omega));
    }
    input.close();

    // calculate some constant values:

    total = 1;
    for (int i = 0; i < SIZE_OF_TABLE; i++) {
        vector<vector<vector<vector<int>* >* >* >* tmpE =
                new vector<vector<vector<vector<int>* >* >* >;
        (*tmpE).reserve(SIZE_OF_TABLE);

        for (int i = 0; i < SIZE_OF_TABLE; i++) {
            vector<vector<vector<int>* >* >* tmpD =
                    new vector<vector<vector<int>* >* >;
            (*tmpD).reserve(CHI_RANGE);

            for (int i = 0; i < CHI_RANGE; i++) {
                vector<vector<int>* >* tmpDD =
                        new vector<vector<int>* >;
                (*tmpDD).reserve(CHI_RANGE);

                for (int i = 0; i < CHI_RANGE; i++) {
                    vector<int>* tmpY = new vector<int>;
                    (*tmpY).reserve(RANGE_OMEGA);

                    for (int i = 0; i < RANGE_OMEGA; i++) {
                        (*tmpY).push_back(0);
                    }
                    (*tmpDD).push_back(tmpY);
                }

                (*tmpD).push_back(tmpDD);
            }
            (*tmpE).push_back(tmpD);
        }
        all_propensities.push_back(tmpE);
    }

    for (int j = 0; j < SIZE_OF_TABLE; j++)
        for (int k = 0; k < SIZE_OF_TABLE; k++)
            for (int l = 0; l < CHI_RANGE; l++)
                for (int n = 0; n < CHI_RANGE; n++)
                    for (int m = 0; m < RANGE_OMEGA; m++)
                        (*(*(*(*all_propensities[j])[k])[l])[n])[m] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        total += amino_count[i];

        for (int j = 0; j < SIZE_OF_TABLE; j++)
            for (int k = 0; k < SIZE_OF_TABLE; k++)
                for (int l = 0; l < CHI_RANGE; l++)
                    for (int n = 0; n < CHI_RANGE; n++)
                        for (int m = 0; m < RANGE_OMEGA; m++)
                            (*(*(*(*all_propensities[j])[k])[l])[n])[m] +=
                                (*(*(*(*(*propensities[i])[j])[k])[l])[n])[m];
    }
    //Construct the max propensities vector
    pConstructMaxPropensities();
    pConstructMinPropensities();
}

/**
 * @Description  Clears all the set data 
 * @param   none
 * @return    changes the object internally (void)  
 */
void PhiPsiOmegaChi1Chi2::pResetData() {

    for (unsigned int i = 0; i < propensities.size(); i++) {
        for (unsigned int j = 0; j < (*propensities[i]).size(); j++) {
            for (unsigned int k = 0; k < (*(*propensities[i])[j]).size(); k++) {
                for (unsigned int l = 0; l < (*(*(*propensities[i])[j])[k]).size(); l++) {
                    for (unsigned int n = 0; n < (*(*(*(*propensities[i])[j])[k])[l]).size(); n++) {
                        (*(*(*(*(*propensities[i])[j])[k])[l])[n]).clear();
                        delete ((*(*(*(*propensities[i])[j])[k])[l])[n]);
                    }
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
        for (unsigned int j = 0; j < (*all_propensities[i]).size(); j++) {
            for (unsigned int k = 0; k < (*(*all_propensities[i])[j]).size(); k++) {
                for (unsigned int l = 0; l < (*(*(*all_propensities[i])[j])[k]).size(); l++) {
                    (*(*(*(*all_propensities[i])[j])[k])[l]).clear();
                    delete ((*(*(*all_propensities[i])[j])[k])[l]);
                }
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
}

// PREDICATES:

/** 
 * @Description calculates the total energy for the amino acids in the spacer 
 * @param spacer reference(Sapcer(&)
 * @return value of the total energy(long double)
 */
long double PhiPsiOmegaChi1Chi2::calculateEnergy(Spacer& sp) {
    long double en = 0.0;
    for (unsigned int i = 1; i < sp.sizeAmino() - 1; i++)
        en += calculateEnergy(sp.getAmino(i));
    return en;
}

/** 
 * @Description calculates the total energy for the amino acids in a portion of the spacer 
 * @param spacer reference(Sapcer(&), positions for start and end of the Sapcer portion (unsigned int, unsigned int) 
 * @return value of the total energy for the aminoacids in the given spacer portion(long double)
 */
long double PhiPsiOmegaChi1Chi2::calculateEnergy(Spacer& sp, unsigned int index1, unsigned int index2) {
    long double en = 0.0;

    index1 = (index1 >= 1) ? index1 : 1;
    index2 = (index2 <= sp.sizeAmino() - 1) ? index2 : sp.sizeAmino() - 1;

    for (unsigned int i = index1; i < index2; i++)
        en += calculateEnergy(sp.getAmino(i));

    return en;
}

/** 
 * @Description calculates the total energy for the amino acid 
 * @param the amino acid reference (AminoAcir&)
 * @return corresponding energy value(long double)
 */
long double PhiPsiOmegaChi1Chi2::calculateEnergy(AminoAcid& aa) {
    // current amino acid type:
    int table_entry = aminoAcidThreeLetterTranslator(aa.getType());
    // offsets for the array
    int a = aa.getMaxChi();
    if (a == 0) {
        int x = sGetPropBin(aa.getPhi(true));
        int y = sGetPropBin(aa.getPsi(true));
        int z = sGetPropChiBin(0.0);
        int n = sGetPropChiBin(0.0);
        int o = sGetPropOmegaBin(aa.getOmega(true));

        return -log((static_cast<double> ((*(*(*(*(*propensities[table_entry])[x])[y])[z])[n])[o])
                / (*(*(*(*all_propensities[x])[y])[z])[n])[o]) / (static_cast<double> (amino_count[table_entry]) / total));
    } else if (a == 1) {
        int x = sGetPropBin(aa.getPhi(true));
        int y = sGetPropBin(aa.getPsi(true));
        int z = sGetPropChiBin(aa.getChi(0));
        int n = sGetPropChiBin(0.0);
        int o = sGetPropOmegaBin(aa.getOmega(true));

        return -log((static_cast<double> ((*(*(*(*(*propensities[table_entry])[x])[y])[z])[n])[o])
                / (*(*(*(*all_propensities[x])[y])[z])[n])[o]) / (static_cast<double>
                (amino_count[table_entry]) / total));
    } else if (a >= 2) {
        int x = sGetPropBin(aa.getPhi(true));
        int y = sGetPropBin(aa.getPsi(true));
        int z = sGetPropChiBin(aa.getChi(0));
        int n = sGetPropChiBin(aa.getChi(1));
        int o = sGetPropOmegaBin(aa.getOmega(true));

        return -log((static_cast<double> ((*(*(*(*(*propensities[table_entry])[x])[y])[z])[n])[o])
                / (*(*(*(*all_propensities[x])[y])[z])[n])[o]) / (static_cast<double>
                (amino_count[table_entry]) / total));
    }
    if (a < 0)
        ERROR("Chi Propension ERROR. Check Chi angles.", exception);

    return (0);
}

/**
 * @Description  Calculates the chi angles index
 * @param   an angle(double)
 * @return   the index for the angle ( int)
 */
inline int PhiPsiOmegaChi1Chi2::sGetPropChiBin(double p) {
    int x = -1;
    if (p > 150)
        x = 0;
    if ((p <= 150) && (p > 100))
        x = 1;
    if ((p <= 100) && (p > 40))
        x = 2;
    if ((p <= 40) && (p > 0))
        x = 3;
    if ((p < 0) && (p > -40))
        x = 3;
    if ((p <= -40) && (p > -100))
        x = 4;
    if ((p <= -100) && (p > -150))
        x = 5;
    if ((p <= -150))
        x = 6;
    if (p == 0)
        x = 7;
    if (x == -1)
        ERROR("ERROR in calculating chi angles. Check chi angles.", exception);
    return x;
}

/**
 * @Description  Adds a propensity for a specific amino acid type
 * @param   code of the amion acid(int), coordinates (int,int,int), values to set (int, int)
 * @return    changes the object internally (void)  
 */
inline void PhiPsiOmegaChi1Chi2::sAddProp(int code, int x, int y, int z, int m, int n) {
    (*(*(*(*(*propensities[code])[x])[y])[z])[m])[n]++;
    amino_count[code]++;
}

/**
 * @Description Sets the step for the arc  
 * @param   the step value(int)
 * @return    changes the object internally (void)  
 */
inline void PhiPsiOmegaChi1Chi2::setArcStep(int n) {
    ARC_STEP = n;
    SIZE_OF_TABLE = 360 / ARC_STEP;
    pResetData();
    pConstructData();
}

/**
 * @Description  obtains and return the propensity binding value
 * @param   angle value(double )
 * @return  the corresponding value ( int)
 */
inline int PhiPsiOmegaChi1Chi2::sGetPropBin(double p) {
    int x = (int) (p + 180) / ARC_STEP;
    if (x == SIZE_OF_TABLE)
        x -= 1; // x is exactly 180 degrees and boosts array
    return x;
}

/**
 * @Description  obtains and return the propensity  Omega binding value
 * @param   angle value(double )
 * @return  the corresponding value ( int)
 */
inline int PhiPsiOmegaChi1Chi2::sGetPropOmegaBin(double p) {
    int x = -1;
    if (RANGE_OMEGA == 11) {
        if ((p < -20) && (p >= -150))
            x = 0;
        if ((p < 150) && (p >= 20))
            x = 1;
        if ((p < -170))
            x = 2;
        if ((p < -160) && (p >= -170))
            x = 3;
        if ((p < -150) && (p >= -160))
            x = 4;
        if ((p < -10) && (p >= -20))
            x = 5;
        if ((p < 10) && (p >= -10))
            x = 6;
        if ((p < 20) && (p >= 10))
            x = 7;
        if ((p < 160) && (p >= 150))
            x = 8;
        if ((p < 170) && (p >= 160))
            x = 9;
        if (p >= 170)
            x = 10;
        return x;
    }
    if (RANGE_OMEGA == 3) {
        if ((p <= -150))
            x = 0;
        if ((p > -150) && (p <= 150))
            x = 1;
        if ((p > 150))
            x = 2;
        return x;
    }
    if (RANGE_OMEGA == 25) {
        if (p <= -175)
            x = 0;
        if ((p <= -170) && (p > -175))
            x = 1;
        if ((p <= -165) && (p > -170))
            x = 2;
        if ((p <= -160) && (p > -165))
            x = 3;
        if ((p <= -155) && (p > -160))
            x = 4;
        if ((p <= -150) && (p > -155))
            x = 5;
        if ((p <= -20) && (p > -150))
            x = 6;
        if ((p <= -15) && (p > -20))
            x = 7;
        if ((p <= -10) && (p > -15))
            x = 8;
        if ((p <= -5) && (p > -10))
            x = 9;
        if ((p <= 0) && (p > -5))
            x = 10;
        if ((p <= 5) && (p > 0))
            x = 11;
        if ((p <= 10) && (p > 5))
            x = 12;
        if ((p <= 15) && (p > 10))
            x = 13;
        if ((p <= 20) && (p > 15))
            x = 14;
        if ((p <= 150) && (p > 20))
            x = 15;
        if ((p <= 155) && (p > 150))
            x = 16;
        if ((p <= 160) && (p > 155))
            x = 17;
        if ((p <= 165) && (p > 160))
            x = 18;
        if ((p <= 170) && (p > 165))
            x = 19;
        if ((p <= 175) && (p > 170))
            x = 20;
        if ((p <= 176) && (p > 175))
            x = 21;
        if ((p <= 177) && (p > 176))
            x = 22;
        if ((p <= 178) && (p > 177))
            x = 23;
        if ((p > 178))
            x = 24;

        return x;
    }

    if (x == -1)
        ERROR("Omega Range not correctly setted.", exception);
    return (0);
}

/**
 * @Description  Sets the value for the Omega range
 * @param   value for the index range( int)
 * @return   the corresponding value for the Range Omega( int)
 */
inline int PhiPsiOmegaChi1Chi2::setRange_Omega(int n) {
    if (n == 1) {
        RANGE_OMEGA = 11;
    } else if (n == 2) {
        RANGE_OMEGA = 3;
    } else if (n == 3) {
        RANGE_OMEGA = 25;
    } else if ((n != 1) && (n != 2) && (n != 3)) {
        ERROR("please insert a valid parameter for omega range.", exception);
    }
    return RANGE_OMEGA;
}

/**
 * @Description  Calculates the Maximum propensities value for the Chi, and Omega Range of the amino acid
 * @param   amino acid index corresponding to the enum (int)
 * @return  corresponding value (double)
 */
double PhiPsiOmegaChi1Chi2::pGetMaxPropensities(int amino) {
    double max = 0.00;

    for (int j = 0; j < SIZE_OF_TABLE; j++) {
        for (int k = 0; k < SIZE_OF_TABLE; k++) {
            for (int l = 0; l < CHI_RANGE; l++) {
                for (int m = 0; m < CHI_RANGE; m++) {
                    for (int x = 0; x < RANGE_OMEGA; x++) {
                        int tmp = (*(*(*(*all_propensities[j])[k])[l])[m])[x];
                        int tmp2 = (*(*(*(*(*propensities[amino])[j])[k])[l])[m])[x];
                        double propensities = ((static_cast<double> (tmp2) / tmp)
                                / ((static_cast<double> (amino_count[amino])) / total));
                        if (propensities > max) {
                            max = propensities;
                        }
                    }
                }
            }
        }
    }
    return max;
}

/**
 * @Description  Calculates the minimum propensities value for the Chi, and Omega Range of the amino acid
 * @param   amino acid index corresponding to the enum (int)
 * @return  corresponding value (double)
 */
double PhiPsiOmegaChi1Chi2::pGetMinPropensities(int amino) {
    double min = 999.99;

    for (int j = 0; j < SIZE_OF_TABLE; j++) {
        for (int k = 0; k < SIZE_OF_TABLE; k++) {
            for (int l = 0; l < CHI_RANGE; l++) {
                for (int m = 0; m < CHI_RANGE; m++) {
                    for (int x = 0; x < RANGE_OMEGA; x++) {
                        int tmp = (*(*(*(*all_propensities[j])[k])[l])[m])[x];
                        int tmp2 = (*(*(*(*(*propensities[amino])[j])[k])[l])[m])[x];
                        double propensities = ((static_cast<double> (tmp2) / tmp)
                                / ((static_cast<double> (amino_count[amino])) / total));
                        if (propensities < min) {
                            min = propensities;
                        }
                    }
                }
            }
        }
    }
    return min;
}

/**
 * @Description  returns the maximum propensities for the corresponding amino
 * @param   amino index from the enum (int)
 * @return   the corresponding value(double)
 */
double PhiPsiOmegaChi1Chi2::pReturnMaxPropensities(int amino) {
    return amino_max_propensities[amino];
}

/**
 * @Description  returns the minimum propensities for the corresponding amino
 * @param   amino index from the enum (int)
 * @return   the corresponding value(double)
 */
double PhiPsiOmegaChi1Chi2::pReturnMinPropensities(int amino) {
    return amino_min_propensities[amino];
}

/**
 * @Description Defines the values for the maximum propensities for each amino acid type
 * @param   none
 * @return    changes the object internally (void)  
 */
void PhiPsiOmegaChi1Chi2::pConstructMaxPropensities() {
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_max_propensities.push_back(0);
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        amino_max_propensities[i] = pGetMaxPropensities(i);
    }
}

/**
 * @Description Defines the values for the minimum propensities for each amino acid type
 * @param   none
 * @return    changes the object internally (void)  
 */
void PhiPsiOmegaChi1Chi2::pConstructMinPropensities() {
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_min_propensities.push_back(0);
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        amino_min_propensities[i] = pGetMinPropensities(i);
    }
}

