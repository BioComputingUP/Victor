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
#include <PhiPsiOmega.h>
#include <float.h>

using namespace Biopool;

// Global constants, typedefs, etc. (to avoid):
vector< vector<Potential::ANGLES> > luca;
vector< vector<Potential::ANGLES> > giovanni;

/**
 * @Description Identifies the sign of a number
 * @param number(double)
 * @return corresponding positive/negative one (int)
 */
int sign(double d) {
    if (d < 0)
        return -1;
    return 1;
}
/**
 * @Description Constructor that sets the arc steps for the table's size
 * @param arc values for table 1 and 2(int ,int), file name where torsion angles are, usually data/tor.par(string)
 */
// CONSTRUCTORS/DESTRUCTOR:

PhiPsiOmega::PhiPsiOmega(int SET_ARC1, string knownledge) :
ARC_STEP(SET_ARC1), TOR_PARAM_FILE(knownledge),
propensities(), all_propensities() {
    SIZE_OF_TABLE = 360 / ARC_STEP;
    RANGE_OMEGA = 3;
    pConstructData();
}

/**
 * @Description  Method that using the torsion file, creates the structure
 * @param none
 * @return changes are made internally(void)
 */
void PhiPsiOmega::pConstructData() {
    //Selecting default knownledge or user knownledge
    char *victor = getenv("VICTOR_ROOT");
    if (victor == NULL)
        ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
    string inputFile;
    if (TOR_PARAM_FILE == "data/tor.par")
        inputFile = getenv("VICTOR_ROOT");

    inputFile += TOR_PARAM_FILE;

    if ((inputFile.length() < 3) && (TOR_PARAM_FILE == "data/tor.par"))
        ERROR("Environment variable VICTOR_ROOT was not found.", exception);

    ifstream input(inputFile.c_str());

    if (!input)
        ERROR("Could not read data file. " + inputFile, exception);

    // Propensities tables construct 

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        vector<vector<vector<int>* >* >* tmpC = new vector<vector<vector<int>* >* >;
        (*tmpC).reserve(SIZE_OF_TABLE);

        for (int j = 0; j < SIZE_OF_TABLE; j++) {
            vector<vector<int>* >* tmpB = new vector<vector<int>* >;
            (*tmpB).reserve(SIZE_OF_TABLE);

            for (int k = 0; k < SIZE_OF_TABLE; k++) {
                vector<int>* tmpA = new vector<int>;
                (*tmpA).reserve(RANGE_OMEGA);

                for (int x = 0; x < RANGE_OMEGA; x++) {
                    (*tmpA).push_back(0);
                }
                (*tmpB).push_back(tmpA);
            }
            (*tmpC).push_back(tmpB);
        }
        (propensities).push_back(tmpC);
    }

    // Propensities table filling

    long numCount;
    double phi, psi, omega;
    double vartmp1, vartmp2; //variable for not include pre angle of TOR_PARAM_FILE   
    string name;
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_count[i] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i += 1)
        for (int j = 0; j < SIZE_OF_TABLE; j++)
            for (int k = 0; k < SIZE_OF_TABLE; k++)
                for (int l = 0; l < RANGE_OMEGA; l++)
                    (*(*(*propensities[i])[j])[k])[l] = 1;

    input >> numCount;

    for (int i = 0; i < numCount; i++) {
        input >> phi >> psi >> name >> vartmp1 >> vartmp2 >> omega;
        skipToNewLine(input);
        sAddProp(aminoAcidThreeLetterTranslator(name), sGetPropBin(phi), sGetPropBin(psi), sGetPropOmegaBin(omega));

    }


    input.close();

    // All propensities table construct

    total = 1;

    for (int i = 0; i < SIZE_OF_TABLE; i++) {
        vector<vector<int>* >* tmpE = new vector<vector<int>* >;
        (*tmpE).reserve(SIZE_OF_TABLE);

        for (int j = 0; j < SIZE_OF_TABLE; j++) {
            vector<int>* tmpD = new vector<int>;
            (*tmpD).reserve(RANGE_OMEGA);

            for (int k = 0; k < RANGE_OMEGA; k++) {
                (*tmpD).push_back(0);
            }
            (*tmpE).push_back(tmpD);
        }
        all_propensities.push_back(tmpE);
    }


    for (int j = 0; j < SIZE_OF_TABLE; j++)
        for (int k = 0; k < SIZE_OF_TABLE; k++)
            for (int l = 0; l < RANGE_OMEGA; l++)
                (*(*all_propensities[j])[k])[l] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        total += amino_count[i];

        for (int j = 0; j < SIZE_OF_TABLE; j++)
            for (int k = 0; k < SIZE_OF_TABLE; k++)
                for (int l = 0; l < RANGE_OMEGA; l++)
                    (*(*all_propensities[j])[k])[l] += (*(*(*propensities[i])[j])[k])[l];
    }

    //Construct the max and min propensity vectors
    pConstructMaxPropensities();
    pConstructMinPropensities();

}

/**
 * @Description Resets all the propensities set previously
 * @param none
 * @return changes are made internally(void)
 */
void PhiPsiOmega::pResetData() {
    for (unsigned int i = 0; i < propensities.size(); i++) {
        for (unsigned int j = 0; j < (*propensities[i]).size(); j++) {
            for (unsigned int k = 0; k < (*(*propensities[i])[j]).size(); k++) {
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
            (*(*all_propensities[i])[j]).clear();
            delete (*all_propensities[i])[j];
        }
        (*all_propensities[i]).clear();
        delete all_propensities[i];
    }
    all_propensities.clear();

}


// PREDICATES:

/**
 * @Description calculates energy for a protein 
 * @param reference to the spacer containing the protein (Spacer&)
 * @return energy value for the set of amino acids in the spacer (long double)
 */
long double PhiPsiOmega::calculateEnergy(Spacer& sp) {
    long double en = 0.0;

    for (unsigned int i = 1; i < sp.sizeAmino() - 1; i++)
        en += calculateEnergy(sp.getAmino(i));


    return en;
}

/**
 * @Description calculates energy for part of the protein
 * @param reference to the spacer that contains the protein(Spacer&) , start and end position of the protein part (unsigned int , unsigned int)
 * @return energy value for the set of amino acids in the selected section (long double)
 */
long double PhiPsiOmega::calculateEnergy(Spacer& sp, unsigned int index1,
        unsigned int index2) {
    long double en = 0.0;

    index1 = (index1 >= 1) ? index1 : 1;
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
long double PhiPsiOmega::calculateEnergy(AminoAcid& aa) {
    // current amino acid type:
    int table_entry = aminoAcidThreeLetterTranslator(aa.getType());

    // offsets for the array
    int x = sGetPropBin(aa.getPhi(true));
    int y = sGetPropBin(aa.getPsi(true));
    int z = sGetPropOmegaBin(aa.getOmega(true));

    if ((x > SIZE_OF_TABLE) || (y > SIZE_OF_TABLE) || (z > RANGE_OMEGA))
        return 0;
    else {

        return -log((static_cast<double> ((*(*(*propensities[table_entry])[x])[y])[z])
                / (*(*all_propensities[x])[y])[z])
                / (static_cast<double> (amino_count[table_entry])
                / total));
    }
}

/**
 *@param diheds: a (phi, psi, omega) triple of dihedral angles.
 * code: an amino acid type.
 *@Return phi-psi-omega torsional energy of the amino acid whose (phi, psi, omega) 
 *    triple is specified by 'diheds' and whose type is specified by code.
 */
long double PhiPsiOmega::calculateEnergy(AminoAcid& diheds, AminoAcidCode code) {

    // current amino acid type:
    int table_entry = code;

    // offsets for the array
    int x = sGetPropBin(diheds.getPhi(true));
    int y = sGetPropBin(diheds.getPsi(true));
    int z = sGetPropOmegaBin(diheds.getOmega(true));

    const int amiProp = (*(*(*propensities[table_entry])[x])[y])[z];
    const int allProp = (*(*all_propensities[x])[y])[z];

    const double propFrac = static_cast<double> (amiProp) / allProp;
    const double countFrac = static_cast<double> (amino_count[table_entry]) / total;

    return -log(propFrac / countFrac);
}

/**
 * @Description
 * @param
 * @return
 */
inline int PhiPsiOmega::sGetPropOmegaBin(double p) {
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
 * @Description Adds propensity value for a amino acid type
 * @param corresponding amino acid code(int), corresponding propensity values(int, int, int, int, int))
 * @return changes are made internally(void)
 */
inline void PhiPsiOmega::sAddProp(int code, int x, int y, int z) {
    (*(*(*propensities[code])[x])[y])[z]++;
    amino_count[code]++;
}

/**
 * @Description
 * @param
 * @return
 */

inline void PhiPsiOmega::setArcStep(int n) {
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
inline int PhiPsiOmega::sGetPropBin(double p) {
    int x = (int) (p + 180) / ARC_STEP;
    if (x == SIZE_OF_TABLE)
        x -= 1; // x is exactly 180 degrees and boosts array
    return x;
}

/**
 * @Description Sets the range for omega
 * @param parameter for omega range, only 1,2, and 3 are posible(int)
 * @return corresponding range for omega (int)
 */
inline int PhiPsiOmega::setRange_Omega(int n) {
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
 * @Description obtains the maximum propensity value for an amino acid type considering all the granularities
 * @param amino acid code(int)
 * @return corresponding propensity value(double)
 */
double PhiPsiOmega::pGetMaxPropensities(int amino) {
    double max = 0.00;

    for (int j = 0; j < SIZE_OF_TABLE; j++) {
        for (int k = 0; k < SIZE_OF_TABLE; k++) {
            for (int x = 0; x < RANGE_OMEGA; x++) {
                int tmp = (*(*all_propensities[j])[k])[x];
                int tmp2 = (*(*(*propensities[amino])[j])[k])[x];
                double propensities = ((static_cast<double> (tmp2) / tmp)
                        / ((static_cast<double> (amino_count[amino])) / total));
                if (propensities > max) {
                    max = propensities;
                }

            }
        }
    }
    return max;
}

/**
 * @Description obtains the maximum propensity value for an amino acid type knowledge based
 * @param amino acid code(int)
 * @return corresponding propensity value(double)
 */
double PhiPsiOmega::pReturnMaxPropensities(int amino) {
    return amino_max_propensities[amino];
}

/**
 * @Description
 * @param
 * @return
 */
void PhiPsiOmega::pConstructMaxPropensities() {
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_max_propensities.push_back(0);

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        amino_max_propensities[i] = pGetMaxPropensities(i);
    }
}

/**
 *@param amino: the code of an amino acid type.
 *@Return value:
 *    minimum (phi, psi, omega) propensity of that type according to kwnowledge.
 */
double PhiPsiOmega::pGetMinPropensities(int amino) {
    const double countFrac = static_cast<double> (amino_count[amino]) / total;
    double min = DBL_MAX;

    for (int j = 0; j < SIZE_OF_TABLE; j++)
        for (int k = 0; k < SIZE_OF_TABLE; k++)
            for (int x = 0; x < RANGE_OMEGA; x++) {
                const int amiProp = (*(*(*propensities[amino])[j])[k])[x];
                const int allProp = (*(*all_propensities[j])[k])[x];

                const double propFrac = static_cast<double> (amiProp) / allProp;

                const double cur = propFrac / countFrac;

                if (cur < min)
                    min = cur;
            }

    return min;
}

/**
 * @Description obtains the minimum propensity value for an amino acid type considering a given phi and psi values knowledge based
 * @param amino acid code(int), values for the previous phi and psi angles(int , int)
 * @return corresponding propensity value(double)
 */
double PhiPsiOmega::pReturnMinPropensities(int amino) {
    return amino_min_propensities[amino];
}

/**
 * @Description computes and stores the minimum propensities based on the max propensities based on knowledge
 * @param none
 * @return changes are made internally(void)
 */
void PhiPsiOmega::pConstructMinPropensities() {
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_min_propensities.push_back(0);

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_min_propensities[i] = pGetMinPropensities(i);
}

/**
 * @Description obtains the Omega angle for a given propensity value and Omega range
 * @param propensity value (int), Omega range(long)
 * @return corresponding omega angle (double)
 */

double PhiPsiOmega::getOmegaAngle(int prop, long RANGE_OMEGA) {
    if (RANGE_OMEGA != 3)
        ERROR("Not defined for OMEGA-RANGE <> 3", exception);

    switch (prop) {
        case 0:
            return -175;
        case 1:
            return 0;
        case 2:
            return 175;
        default:
            ERROR("Wrong propensity for omega angle", exception);
    }
}

/**
 * @Description Calculates the value os the psi , phi angle
 * @param propensity value(int), long(not used),arc step(int)
 * @return
 */
inline double PhiPsiOmega::getPhiPsiAngle(int prop, long SIZE_TABLE, int ARC_STEP) {
    return ARC_STEP * prop - 180;
}

/**
 * @Description Returns asort energy table
 * @param none
 * @return pointer to new ordered table (vector< vector<Potential::ANGLES> >*)
 */
vector< vector<Potential::ANGLES> >* PhiPsiOmega::getOrderedEnergyTable() {
    vector< vector<ANGLES> > *arr = new vector< vector<ANGLES> >;

    for (int aa = 0; aa < AminoAcid_CODE_SIZE; ++aa) {
        // get the propensities for each amino acid:
        vector< vector< vector<int>* >* >* data = propensities[aa];
        vector<ANGLES> elem;
        for (int x = 0; x < SIZE_OF_TABLE; ++x)
            for (int y = 0; y < SIZE_OF_TABLE; ++y)
                for (int z = 0; z < RANGE_OMEGA; ++z) {
                    ANGLES tmp = {getPhiPsiAngle(x, SIZE_OF_TABLE, ARC_STEP),
                        getPhiPsiAngle(y, SIZE_OF_TABLE, ARC_STEP),
                        getOmegaAngle(z, RANGE_OMEGA), (*(*(*data)[x])[y])[z]};
                    elem.push_back(tmp);
                }
        // Let's sort the array
        sort(elem.begin(), elem.end(), EnergyGreater());
        arr->push_back(elem);
    }

    return arr;
}






