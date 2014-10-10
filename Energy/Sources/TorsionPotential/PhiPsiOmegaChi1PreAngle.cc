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
#include <PhiPsiOmegaChi1PreAngle.h>

using namespace Victor;

// Global constants, typedefs, etc. (to avoid):

// CONSTRUCTORS/DESTRUCTOR:

PhiPsiOmegaChi1PreAngle::PhiPsiOmegaChi1PreAngle(int SET_ARC1, int SET_ARC2, string knownledge) :
TOR_PARAM_FILE(knownledge), ARC_STEP(SET_ARC1), ARC_STEP2(SET_ARC2), propensities(), all_propensities() {
    SIZE_OF_TABLE = 360 / ARC_STEP;
    CHI_RANGE = 8;
    RANGE_OMEGA = 3;
    SIZE_OF_TABLE2 = 360 / ARC_STEP2;
    pConstructData();
}

void PhiPsiOmegaChi1PreAngle::pConstructData() {
    char *victor = getenv("VICTOR_ROOT");
    if (victor == NULL)
        ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
    string inputFile;
    if (TOR_PARAM_FILE == "data/tor.par")
        inputFile = getenv("VICTOR_ROOT");

    inputFile += TOR_PARAM_FILE;

    if (inputFile.length() < 3)
        ERROR("Environment variable VICTOR_ROOT was not found.", exception);

    ifstream input(inputFile.c_str());

    if (!input)
        ERROR("Could not read data file. " + inputFile, exception);

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        vector<vector<vector<vector<vector<vector<int>* >* >* >* >* >* tmpA =
                new vector<vector<vector<vector<vector<vector<int>* >* >* >* > *>;
        (*tmpA).reserve(SIZE_OF_TABLE);

        for (int i = 0; i < SIZE_OF_TABLE; i++) {
            vector<vector<vector<vector<vector<int>* >* >* >* >* tmpB =
                    new vector<vector<vector<vector<vector<int>* >* >* >* >;
            (*tmpB).reserve(SIZE_OF_TABLE);

            for (int i = 0; i < SIZE_OF_TABLE; i++) {
                vector<vector<vector<vector<int>* >* >* >* tmpC =
                        new vector<vector<vector<vector<int>* >* >* >;
                (*tmpC).reserve(SIZE_OF_TABLE2);

                for (int i = 0; i < SIZE_OF_TABLE2; i++) {
                    vector<vector<vector<int>* >* >* tmpD =
                            new vector<vector<vector<int>* >* >;
                    (*tmpD).reserve(SIZE_OF_TABLE2);

                    for (int i = 0; i < SIZE_OF_TABLE2; i++) {
                        vector<vector<int>* >* tmpE =
                                new vector<vector<int>* >;
                        (*tmpE).reserve(CHI_RANGE);

                        for (int i = 0; i < CHI_RANGE; i++) {
                            vector<int>* tmpF = new vector<int>;
                            (*tmpF).reserve(RANGE_OMEGA);

                            for (int i = 0; i < RANGE_OMEGA; i++) {
                                (*tmpF).push_back(0);
                            }
                            (*tmpE).push_back(tmpF);
                        }
                        (*tmpD).push_back(tmpE);
                    }
                    (*tmpC).push_back(tmpD);
                }
                (*tmpB).push_back(tmpC);
            }
            (*tmpA).push_back(tmpB);
        }
        (propensities).push_back(tmpA);
    }

    long numCount;
    double phi, psi, chi1, omega, numchi, prephi, prepsi;
    string name;
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_count[i] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i += 1)
        for (int j = 0; j < SIZE_OF_TABLE; j++)
            for (int k = 0; k < SIZE_OF_TABLE; k++)
                for (int h = 0; h < SIZE_OF_TABLE2; h++)
                    for (int n = 0; n < SIZE_OF_TABLE2; n++)
                        for (int l = 0; l < CHI_RANGE; l++)
                            for (int m = 0; m < RANGE_OMEGA; m++)
                                (*(*(*(*(*(*propensities[i])[j])[k])[h])[n])[l])[m] = 1;

    input >> numCount;

    for (int i = 0; i < numCount; i++) {
        input >> phi >> psi >> name >> prephi >> prepsi >> omega >> numchi;

        if (numchi == 0) {
            chi1 = 0;
        } else {
            input >> chi1;
            skipToNewLine(input);
        }
        sAddProp(aminoAcidThreeLetterTranslator(name), sGetPropBin(phi),
                sGetPropBin(psi), sGetPropBin2(prephi), sGetPropBin2(prepsi),
                sGetPropChiBin(chi1), sGetPropOmegaBin(omega));
    }

    input.close();

    // calculate some constant values:

    total = 1;

    for (int i = 0; i < SIZE_OF_TABLE; i++) {
        vector<vector<vector<vector<vector<int>* >* >* >* >* tmpAA =
                new vector<vector<vector<vector<vector<int>* >* >* >* >;
        (*tmpAA).reserve(SIZE_OF_TABLE);

        for (int i = 0; i < SIZE_OF_TABLE; i++) {
            vector<vector<vector<vector<int>* >* >* >* tmpBB =
                    new vector<vector<vector<vector<int>* >* >* >;
            (*tmpBB).reserve(SIZE_OF_TABLE2);

            for (int i = 0; i < SIZE_OF_TABLE2; i++) {
                vector<vector<vector<int>* >* >* tmpCC =
                        new vector<vector<vector<int>* >* >;
                (*tmpCC).reserve(SIZE_OF_TABLE2);

                for (int i = 0; i < SIZE_OF_TABLE2; i++) {
                    vector<vector<int>* >* tmpDD =
                            new vector<vector<int>* >;
                    (*tmpDD).reserve(CHI_RANGE);

                    for (int i = 0; i < CHI_RANGE; i++) {
                        vector<int>* tmpEE =
                                new vector<int>;
                        (*tmpEE).reserve(RANGE_OMEGA);

                        for (int i = 0; i < RANGE_OMEGA; i++) {
                            (*tmpEE).push_back(0);
                        }
                        (*tmpDD).push_back(tmpEE);
                    }

                    (*tmpCC).push_back(tmpDD);
                }
                (*tmpBB).push_back(tmpCC);
            }
            (*tmpAA).push_back(tmpBB);
        }
        all_propensities.push_back(tmpAA);
    }


    for (int j = 0; j < SIZE_OF_TABLE; j++)
        for (int k = 0; k < SIZE_OF_TABLE; k++)
            for (int n = 0; n < SIZE_OF_TABLE2; n++)
                for (int o = 0; o < SIZE_OF_TABLE2; o++)
                    for (int l = 0; l < CHI_RANGE; l++)
                        for (int m = 0; m < RANGE_OMEGA; m++)
                            (*(*(*(*(*all_propensities[j])[k])[n])[o])[l])[m] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        total += amino_count[i];

        for (int j = 0; j < SIZE_OF_TABLE; j++)
            for (int k = 0; k < SIZE_OF_TABLE; k++)
                for (int n = 0; n < SIZE_OF_TABLE2; n++)
                    for (int o = 0; o < SIZE_OF_TABLE2; o++)
                        for (int l = 0; l < CHI_RANGE; l++)
                            for (int m = 0; m < RANGE_OMEGA; m++)
                                (*(*(*(*(*all_propensities[j])[k])[n])[o])[l])[m] +=
                                    (*(*(*(*(*(*propensities[i])[j])[k])[n])[o])[l])[m];
    }

    //Construct the max propensities vector
    pConstructMaxPropensities();

}

/**
 * @Description resets the data
 * @param  none
 */
void PhiPsiOmegaChi1PreAngle::pResetData() {

    for (unsigned int i = 0; i < propensities.size(); i++) {
        for (unsigned int j = 0; j < (*propensities[i]).size(); j++) {
            for (unsigned int k = 0; k < (*(*propensities[i])[j]).size(); k++) {
                for (unsigned int l = 0; l < (*(*(*propensities[i])[j])[k]).size(); l++) {
                    for (unsigned int m = 0; m < (*(*(*(*propensities[i])[j])[k])[l]).size(); m++) {
                        for (unsigned int n = 0; n < (*(*(*(*(*propensities[i])[j])[k])[l])[m]).size(); n++) {
                            (*(*(*(*(*(*propensities[i])[j])[k])[l])[m])[n]).clear();
                            delete((*(*(*(*(*propensities[i])[j])[k])[l])[m])[n]);
                        }

                        (*(*(*(*(*propensities[i])[j])[k])[l])[m]).clear();
                        delete((*(*(*(*propensities[i])[j])[k])[l])[m]);
                    }

                    (*(*(*(*propensities[i])[j])[k])[l]).clear();
                    delete((*(*(*propensities[i])[j])[k])[l]);
                }

                (*(*(*propensities[i])[j])[k]).clear();
                delete((*(*propensities[i])[j])[k]);

            }

            (*(*propensities[i])[j]).clear();
            delete((*propensities[i])[j]);

        }
        (*propensities[i]).clear();
        delete(propensities[i]);
    }
    propensities.clear();

    for (unsigned int i = 0; i < all_propensities.size(); i++) {
        for (unsigned int j = 0; j < (*all_propensities[i]).size(); j++) {
            for (unsigned int k = 0; k < (*(*all_propensities[i])[j]).size(); k++) {
                for (unsigned int l = 0; l < (*(*(*all_propensities[i])[j])[k]).size(); l++) {
                    for (unsigned int m = 0; m < (*(*(*(*all_propensities[i])[j])[k])[l]).size(); m++) {
                        (*(*(*(*(*all_propensities[i])[j])[k])[l])[m]).clear();
                        delete((*(*(*(*all_propensities[i])[j])[k])[l])[m]);
                    }

                    (*(*(*(*all_propensities[i])[j])[k])[l]).clear();
                    delete((*(*(*all_propensities[i])[j])[k])[l]);
                }

                (*(*(*all_propensities[i])[j])[k]).clear();
                delete((*(*all_propensities[i])[j])[k]);

            }

            (*(*all_propensities[i])[j]).clear();
            delete((*all_propensities[i])[j]);

        }

        (*all_propensities[i]).clear();
        delete(all_propensities[i]);
    }

    all_propensities.clear();

    for (unsigned int i = 0; i < amino_max_propensities_pre_angle.size(); i++) {
        for (unsigned int j = 0; j < (*amino_max_propensities_pre_angle[i]).size(); j++) {
            (*(*amino_max_propensities_pre_angle[i])[j]).clear();
            delete((*amino_max_propensities_pre_angle[i])[j]);
        }
        (*amino_max_propensities_pre_angle[i]).clear();
        delete(amino_max_propensities_pre_angle[i]);
    }
    amino_max_propensities_pre_angle.clear();

}



// PREDICATES:

/**
 * @Description Calculates energy 
 * @param spacer reference
 * @return long double
 */
long double PhiPsiOmegaChi1PreAngle::calculateEnergy(Spacer& sp) {
    long double en = 0.0;

    for (unsigned int i = 2; i < sp.sizeAmino() - 1; i++)
        en += calculateEnergy(sp.getAmino(i));

    return en;
}

/**
 * @Description calculates energy
 * @param spacer reference ,unsigned int 
 */
long double PhiPsiOmegaChi1PreAngle::calculateEnergy(Spacer& sp, unsigned int index1,
        unsigned int index2) {
    long double en = 0.0;

    index1 = (index1 >= 2) ? index1 : 2;
    index2 = (index2 <= sp.sizeAmino() - 1) ? index2 : sp.sizeAmino() - 1;

    for (unsigned int i = index1; i < index2; i++)
        en += calculateEnergy(sp.getAmino(i));

    return en;
}

/**
 * @Description calculates energy
 * @param amino acid reference
 * @return long double
 */
long double PhiPsiOmegaChi1PreAngle::calculateEnergy(AminoAcid& aa) {
    // current amino acid type:
    int table_entry = aminoAcidThreeLetterTranslator(aa.getType());

    int t = aa.sizeInBonds();
    if (t == 0)
        ERROR("Invalid Aminoacid, impossible to calculate pre-angle.", exception);

    // offsets for the array
    int a = aa.getMaxChi();
    if (a == 0) {
        int x = sGetPropBin(aa.getPhi(true));
        int y = sGetPropBin(aa.getPsi(true));
        int z = sGetPropChiBin(0.0);
        int o = sGetPropOmegaBin(aa.getOmega(true));
        int m = sGetPropBin2(aa.getInBond(0).getPhi());
        int n = sGetPropBin2(aa.getInBond(0).getPsi());

        return -log((static_cast<double> ((*(*(*(*(*(*propensities[table_entry])[x])[y])[m])[n])[z])[o])
                / (*(*(*(*(*all_propensities[x])[y])[m])[n])[z])[o])
                / (static_cast<double> (amino_count[table_entry])
                / total));
    } else {
        int x = sGetPropBin(aa.getPhi(true));
        int y = sGetPropBin(aa.getPsi(true));
        int z = sGetPropChiBin(aa.getChi(0));
        int o = sGetPropOmegaBin(aa.getOmega(true));
        int m = sGetPropBin2(aa.getInBond(0).getPhi());
        int n = sGetPropBin2(aa.getInBond(0).getPsi());

        return -log((static_cast<double> ((*(*(*(*(*(*propensities[table_entry])[x])[y])[m])[n])[z])[o])
                / (*(*(*(*(*all_propensities[x])[y])[m])[n])[z])[o])
                / (static_cast<double> (amino_count[table_entry])
                / total));
    }
}

/**
 * @Description Identify the correct chi 
 * @param double
 * @return int
 */
inline int PhiPsiOmegaChi1PreAngle::sGetPropChiBin(double p) {
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
        ERROR("Chi Propension ERROR. Check Chi angles.", exception);

    return x;
}

inline void PhiPsiOmegaChi1PreAngle::sAddProp(int code, int x, int y, int z, int m, int n, int o) {
    (*(*(*(*(*(*propensities[code])[x])[y])[z])[m])[n])[o]++;
    amino_count[code]++;
}

inline void PhiPsiOmegaChi1PreAngle::setArcStep(int n) {
    ARC_STEP = n;
    SIZE_OF_TABLE = 360 / ARC_STEP;
    pResetData();
    pConstructData();
}

inline int PhiPsiOmegaChi1PreAngle::sGetPropBin(double p) {
    int x = (int) (p + 180) / ARC_STEP;
    if (x == SIZE_OF_TABLE)
        x -= 1; // x is exactly 180 degrees and boosts array
    return x;
}

inline int PhiPsiOmegaChi1PreAngle::sGetPropBin2(double p) {
    int x = (int) (p + 180) / ARC_STEP2;
    if (x == SIZE_OF_TABLE2)
        x -= 1; // x is exactly 180 degrees and boosts array
    return x;
}

inline int PhiPsiOmegaChi1PreAngle::sGetPropOmegaBin(double p) {
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

inline int PhiPsiOmegaChi1PreAngle::setRange_Omega(int n) {
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

double PhiPsiOmegaChi1PreAngle::pGetMaxPropensities(int amino) {
    double max = 0.00;

    for (int j = 0; j < SIZE_OF_TABLE; j++) {
        for (int k = 0; k < SIZE_OF_TABLE; k++) {
            for (int r = 0; r < SIZE_OF_TABLE2; r++) {
                for (int s = 0; s < SIZE_OF_TABLE2; s++) {
                    for (int l = 0; l < CHI_RANGE; l++) {
                        for (int x = 0; x < RANGE_OMEGA; x++) {
                            int tmp = (*(*(*(*(*all_propensities[j])[k])[r])[s])[l])[x];
                            int tmp2 = (*(*(*(*(*(*propensities[amino])[j])[k])[r])[s])[l])[x];
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
    }

    return max;
}

double PhiPsiOmegaChi1PreAngle::pGetMaxPropensities(int amino, int prephi, int prepsi) {
    double max = 0.00;

    for (int j = 0; j < SIZE_OF_TABLE; j++) {
        for (int k = 0; k < SIZE_OF_TABLE; k++) {
            for (int l = 0; l < CHI_RANGE; l++) {
                for (int x = 0; x < RANGE_OMEGA; x++) {
                    int tmp = (*(*(*(*(*all_propensities[j])[k])[prephi])[prepsi])[l])[x];
                    int tmp2 = (*(*(*(*(*(*propensities[amino])[j])[k])[prephi])[prepsi])[l])[x];
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

double PhiPsiOmegaChi1PreAngle::pReturnMaxPropensities(int amino) {
    return amino_max_propensities[amino];
}

double PhiPsiOmegaChi1PreAngle::pReturnMaxPropensitiesPreAngle(int amino,
        int prephi, int prepsi) {
    return (*(*amino_max_propensities_pre_angle[amino])[prephi])[prepsi];
}

void PhiPsiOmegaChi1PreAngle::pConstructMaxPropensities() {
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

