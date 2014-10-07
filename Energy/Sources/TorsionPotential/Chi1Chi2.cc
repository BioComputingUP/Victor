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
/**
 *@Class               Chi1Chi2
 *@Project        Victor
 *@Description 
 *    This class implements a simple torsion potential based on the 
 *    statistical preference of aminoacid types for certain Chi1 and Chi2 angles.
 */

// Includes:
#include <Chi1Chi2.h>

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:

Chi1Chi2::Chi1Chi2(string knownledge) :
TOR_PARAM_FILE(knownledge),
propensities(), all_propensities() {
    CHI_RANGE = 8;
    pConstructData();
}

void Chi1Chi2::pConstructData() {
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
        vector<vector<int >* >* tmpB = new vector<vector<int >* >;
        (*tmpB).reserve(CHI_RANGE);

        for (int j = 0; j < CHI_RANGE; j++) {
            vector<int >* tmpA = new vector<int >;
            (*tmpA).reserve(CHI_RANGE);

            for (int z = 0; z < CHI_RANGE; z++) {
                (*tmpA).push_back(0);
            }
            (*tmpB).push_back(tmpA);
        }
        propensities.push_back(tmpB);
    }

    // Propensities table filling

    long numCount;
    double phi, psi, omega, chi1, chi2;
    double tmp1, tmp2;
    string name;
    int numchi;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_count[i] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i += 1)
        for (int j = 0; j < CHI_RANGE; j++)
            for (int k = 0; k < CHI_RANGE; k++)
                (*(*propensities[i])[j])[k] = 1;

    input >> numCount;

    for (int i = 0; i < numCount; i++) {
        input >> phi >> psi >> name >> tmp1 >> tmp2 >> omega >> numchi;

        if (numchi != 0) {
            if (numchi == 1) {
                input >> chi1;
                chi2 = 0.0;
            }
            if (numchi == 2) {
                input >> chi1 >> chi2;
            }
            if (numchi >= 3) {
                input >> chi1 >> chi2;
                skipToNewLine(input);
            }
        } else if (numchi == 0) {
            chi1 = 0.0;
            chi2 = 0.0;
        }
        sAddProp(aminoAcidThreeLetterTranslator(name), sGetPropChiBin(chi1),
                sGetPropChiBin(chi2));
    }

    input.close();

    // calculate some constant values:

    total = 1;

    for (int i = 0; i < CHI_RANGE; i++) {
        vector<int>* tmp = new vector<int>;
        (*tmp).reserve(CHI_RANGE);

        for (int j = 0; j < CHI_RANGE; j++) {
            (*tmp).push_back(0);
        }
        all_propensities.push_back(tmp);
    }

    //All propensities tables filling

    for (int j = 0; j < CHI_RANGE; j++)
        for (int k = 0; k < CHI_RANGE; k++)
            (*all_propensities[j])[k] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        total += amino_count[i];

        for (int j = 0; j < CHI_RANGE; j++)
            for (int k = 0; k < CHI_RANGE; k++)
                (*all_propensities[j])[k] += (*(*propensities[i])[j])[k];
    }

    //Construct the max propensities vector
    pConstructMaxPropensities();


}

void
Chi1Chi2::pResetData() {
    //Reset data for propensieties and all propensities tables

    for (unsigned int i = 0; i < propensities.size(); i++) {
        for (unsigned int j = 0; j < (*propensities[i]).size(); j++) {
            (*(*propensities[i])[j]).clear();
            delete ((*propensities[i])[j]);
        }
        (*propensities[i]).clear();
        delete (propensities[i]);
    }
    propensities.clear();


    for (unsigned int i = 0; i < all_propensities.size(); i++) {
        (*all_propensities[i]).clear();
        delete (all_propensities[i]);
    }
    all_propensities.clear();

}



// PREDICATES:

long double
Chi1Chi2::calculateEnergy(Spacer& sp) {
    long double en = 0.0;

    for (unsigned int i = 1; i < sp.sizeAmino() - 1; i++)
        en += calculateEnergy(sp.getAmino(i));

    return en;
}

long double
Chi1Chi2::calculateEnergy(Spacer& sp, unsigned int index1,
        unsigned int index2) {
    long double en = 0.0;

    index1 = (index1 >= 1) ? index1 : 1;
    index2 = (index2 <= sp.sizeAmino() - 1) ? index2 : sp.sizeAmino() - 1;

    for (unsigned int i = index1; i < index2; i++)
        en += calculateEnergy(sp.getAmino(i));

    return en;
}

long double
Chi1Chi2::calculateEnergy(AminoAcid& aa) {
    // current amino acid type:
    int table_entry = aminoAcidThreeLetterTranslator(aa.getType());

    // offsets for the array:
    int a = aa.getMaxChi();
    if (a == 0) {
        int x = sGetPropChiBin(0.00);
        int y = x;
        return -log((static_cast<double> ((*(*propensities[table_entry])[x])[y])
                / (*all_propensities[x])[y])
                / (static_cast<double> (amino_count[table_entry])
                / total));
    }
    if (a == 1) {
        int x = sGetPropChiBin(aa.getChi(0));
        int y = sGetPropChiBin(0.00);
        return -log((static_cast<double> ((*(*propensities[table_entry])[x])[y])
                / (*all_propensities[x])[y])
                / (static_cast<double> (amino_count[table_entry])
                / total));
    }
    if (a >= 2) {
        int x = sGetPropChiBin(aa.getChi(0));
        int y = sGetPropChiBin(aa.getChi(1));
        return -log((static_cast<double> ((*(*propensities[table_entry])[x])[y])
                / (*all_propensities[x])[y])
                / (static_cast<double> (amino_count[table_entry])
                / total));
    }
    if (a < 0)
        ERROR("Chi Propension ERROR. Check Chi angles.", exception);

    return (0);
}

inline int
Chi1Chi2::sGetPropChiBin(double p) {
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
    if (p == 0.00)
        x = 7;

    if (x == -1)
        ERROR("ERROR in calculating chi angles. Check chi angles.", exception);

    return x;
}

inline void
Chi1Chi2::sAddProp(int code, int x, int y) {
    (*(*propensities[code])[x])[y]++;
    amino_count[code]++;
}

double
Chi1Chi2::pGetMaxPropensities(int amino) {
    double max = 0.00;

    for (int j = 0; j < CHI_RANGE; j++) {
        for (int k = 0; k < CHI_RANGE; k++) {
            int tmp = (*all_propensities[j])[k];
            int tmp2 = (*(*propensities[amino])[j])[k];
            double propensities = ((static_cast<double> (tmp2) / tmp)
                    / ((static_cast<double> (amino_count[amino])) / total));
            if (propensities > max) {
                max = propensities;
            }

        }
    }
    return max;
}

double
Chi1Chi2::pReturnMaxPropensities(int amino) {
    return amino_max_propensities[amino];
}

void
Chi1Chi2::pConstructMaxPropensities() {
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_max_propensities.push_back(0);

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        amino_max_propensities[i] = pGetMaxPropensities(i);
    }
}















