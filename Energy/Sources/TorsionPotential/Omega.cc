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
 *@Description:        see .h file
 */
// Includes:
#include <Omega.h>
using namespace Biopool;
// Global constants, typedefs, etc. (to avoid):
// CONSTRUCTORS/DESTRUCTOR:

Omega::Omega(string knownledge) : TOR_PARAM_FILE(knownledge), propensities(), all_propensities() {
    RANGE_OMEGA = 3;
    pConstructData();
}

void Omega::pConstructData() {
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
        vector<int >* tmpA = new vector<int >;
        (*tmpA).reserve(RANGE_OMEGA);

        for (int z = 0; z < RANGE_OMEGA; z++)
            (*tmpA).push_back(0);

        propensities.push_back(tmpA);
    }

    long numCount;
    double phi, psi, omega;
    double tmp1, tmp2;
    string name;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_count[i] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i += 1)
        for (int j = 0; j < RANGE_OMEGA; j++)
            (*propensities[i])[j] = 1;

    input >> numCount;

    for (int i = 0; i < numCount; i++) {
        input >> phi >> psi >> name >> tmp1 >> tmp2 >> omega;
        skipToNewLine(input);
        sAddProp(aminoAcidThreeLetterTranslator(name), sGetPropOmegaBin(omega));
    }
    input.close();
    // calculate some constant values:
    total = 1;
    for (int i = 0; i < RANGE_OMEGA; i++) {
        (all_propensities).push_back(0);
    }

    for (int j = 0; j < RANGE_OMEGA; j++)
        all_propensities[j] = 1;
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        total += amino_count[i];

        for (int j = 0; j < RANGE_OMEGA; j++)
            all_propensities[j] += (*propensities[i])[j];
    }
    //Construct the max propensities vector
    pConstructMaxPropensities();

}

void Omega::pResetData() {
    for (unsigned int i = 0; i < propensities.size(); i++) {
        (*propensities[i]).clear();
        delete (propensities[i]);
    }
    propensities.clear();

    all_propensities.clear();

}

// PREDICATES:

long double Omega::calculateEnergy(Spacer& sp) {
    long double en = 0.0;

    for (unsigned int i = 1; i < sp.sizeAmino() - 1; i++)
        en += calculateEnergy(sp.getAmino(i));

    return en;
}

long double Omega::calculateEnergy(Spacer& sp, unsigned int index1, unsigned int index2) {
    long double en = 0.0;

    index1 = (index1 >= 1) ? index1 : 1;
    index2 = (index2 <= sp.sizeAmino() - 1) ? index2 : sp.sizeAmino() - 1;

    for (unsigned int i = index1; i < index2; i++)
        en += calculateEnergy(sp.getAmino(i));

    return en;
}

long double Omega::calculateEnergy(AminoAcid& aa) {
    // current amino acid type:
    int table_entry = aminoAcidThreeLetterTranslator(aa.getType());
    // offsets for the array
    int x = sGetPropOmegaBin(aa.getOmega(true));

    return -log(((static_cast<double> ((*propensities[table_entry])[x]))
            / (all_propensities[x]))
            / (static_cast<double> (amino_count[table_entry])
            / total));
}

inline int Omega::sGetPropOmegaBin(double p) {
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

inline void Omega::sAddProp(int code, int x) {
    (*propensities[code])[x]++;
    amino_count[code]++;
}

inline void Omega::setRange_Omega(int n) {
    if (n == 1) {
        RANGE_OMEGA = 11;
    } else if (n == 2) {
        RANGE_OMEGA = 3;
    } else if (n == 3) {
        RANGE_OMEGA = 25;
    } else if ((n != 1) && (n != 2) && (n != 3)) {
        ERROR("please insert a valid parameter for omega range.", exception);
    }

    pResetData();
    pConstructData();
}

double Omega::pGetMaxPropensities(int amino) {
    double max = 0.00;

    for (int j = 0; j < RANGE_OMEGA; j++) {
        int tmp = all_propensities[j];
        int tmp2 = (*propensities[amino])[j];
        double propensities = ((static_cast<double> (tmp2) / tmp)
                / ((static_cast<double> (amino_count[amino])) / total));
        if (propensities > max) {
            max = propensities;
        }
    }
    return max;
}

double Omega::pReturnMaxPropensities(int amino) {
    return amino_max_propensities[amino];
}

void Omega::pConstructMaxPropensities() {
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_max_propensities.push_back(0);
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        amino_max_propensities[i] = pGetMaxPropensities(i);
    }
}


