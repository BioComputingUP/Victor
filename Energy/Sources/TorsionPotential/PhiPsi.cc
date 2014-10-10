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
#include <PhiPsi.h>

using namespace Victor::Biopool;
using namespace Victor::Energy;
using namespace Victor;

/**
 * @Description  
 * @param   
 * @return   
 */
PhiPsi::PhiPsi(int SET_ARC1, string knownledge) : TOR_PARAM_FILE(knownledge), ARC_STEP(SET_ARC1),
propensities(), all_propensities() {
    SIZE_OF_TABLE = 360 / ARC_STEP;
    pConstructData();
}

/**
 * @Description  Sets all the information need for the class
 * @param   none
 * @return    changes the object internally (void)   
 */
void PhiPsi::pConstructData() {
    //Selecting default knownledge or user knownledge
    char *victor = getenv("VICTOR_ROOT");
    if (victor == NULL)
        ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
    string inputFile;
    if (TOR_PARAM_FILE == "data/tor.par") {
        inputFile = getenv("VICTOR_ROOT");
        inputFile += TOR_PARAM_FILE;
    } else {
        inputFile = getenv("VICTOR_ROOT");
        string path = "data/tor.par";
        inputFile += path;
    }


    if ((inputFile.length() < 3) && (TOR_PARAM_FILE == "data/tor.par"))
        ERROR("Environment variable VICTOR_ROOT was not found.", exception);

    ifstream input(inputFile.c_str());

    if (!input)
        ERROR("Could not read data file. " + inputFile, exception);

    // Propensities tables construct

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        vector<vector<int >* >* tmpB = new vector<vector<int >* >;
        (*tmpB).reserve(SIZE_OF_TABLE);

        for (int j = 0; j < SIZE_OF_TABLE; j++) {
            vector<int >* tmpA = new vector<int >;
            (*tmpA).reserve(SIZE_OF_TABLE);

            for (int z = 0; z < SIZE_OF_TABLE; z++) {
                (*tmpA).push_back(0);
            }
            (*tmpB).push_back(tmpA);
        }
        propensities.push_back(tmpB);
    }

    // Propensities table filling


    long numCount;
    double phi, psi;
    string name;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_count[i] = 1;

    for (int i = 0; i < AminoAcid_CODE_SIZE; i += 1)
        for (int j = 0; j < SIZE_OF_TABLE; j++)
            for (int k = 0; k < SIZE_OF_TABLE; k++)
                (*(*propensities[i])[j])[k] = 1;

    input >> numCount;

    for (int i = 0; i < numCount; i++) {
        input >> phi >> psi >> name;
        skipToNewLine(input);
        sAddProp(aminoAcidThreeLetterTranslator(name), sGetPropBin(phi),
                sGetPropBin(psi));
    }

    input.close();

    // All propensities table construct

    total = 1;

    for (int i = 0; i < SIZE_OF_TABLE; i++) {
        vector<int>* tmp = new vector<int>;
        (*tmp).reserve(SIZE_OF_TABLE);

        for (int j = 0; j < SIZE_OF_TABLE; j++) {
            (*tmp).push_back(0);
        }
        all_propensities.push_back(tmp);
    }

    for (int i = 0; i < SIZE_OF_TABLE; i++)
        for (int j = 0; j < SIZE_OF_TABLE; j++)
            (*all_propensities[i])[j] = 1;

    //All propensities tables filling

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        total += amino_count[i];

        for (int j = 0; j < SIZE_OF_TABLE; j++)
            for (int k = 0; k < SIZE_OF_TABLE; k++)
                (*all_propensities[j])[k] += (*(*propensities[i])[j])[k];
    }

    //Construct the max/min propensities vector
    pConstructMaxPropensities();
    pConstructMinPropensities();

}

/**
 * @Description  Clears all the set data 
 * @param   none
 * @return    changes the object internally (void)  
 */
void PhiPsi::pResetData() {
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

/** 
 * @Description calculates the total energy for the amino acids in the spacer 
 * @param spacer reference(Sapcer(&)
 * @return value of the total energy(long double)
 */
long double PhiPsi::calculateEnergy(Spacer& sp) {
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
long double PhiPsi::calculateEnergy(Spacer& sp, unsigned int index1, unsigned int index2) {
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
long double PhiPsi::calculateEnergy(AminoAcid& aa) {
    // current amino acid type:
    int table_entry = aminoAcidThreeLetterTranslator(aa.getType());

    // offsets for the array
    int x = sGetPropBin(aa.getPhi(true));
    int y = sGetPropBin(aa.getPsi(true));

    if ((x > SIZE_OF_TABLE) || (y > SIZE_OF_TABLE))
        return 0;
    else {

        return -log((static_cast<double> ((*(*propensities[table_entry])[x])[y])
                / (*all_propensities[x])[y])
                / (static_cast<double> (amino_count[table_entry])
                / total));
    }
}

/** 
 * @Description calculates the total energy for the amino acid 
 * @param the amino acid reference (AminoAcir&)
 * @return corresponding energy value(long double)
 */
long double PhiPsi::calculateEnergySmooth(AminoAcid& aa) {

    // current amino acid type:
    int table_entry = aminoAcidThreeLetterTranslator(aa.getType());

    // offsets for the array
    double diffX, diffY;
    int x = sGetPropBinDiff(aa.getPhi(true), diffX);
    int y = sGetPropBinDiff(aa.getPsi(true), diffY);

    if ((x > SIZE_OF_TABLE) || (y > SIZE_OF_TABLE))
        return 0;

    long double en1 = -log((static_cast<double> (
            (*(*propensities[table_entry])[x])[y])
            / (*all_propensities[x])[y])
            / (static_cast<double> (amino_count[table_entry])
            / total));

    long double enX = 0.0;
    if (diffX >= 0)
        enX = -log((static_cast<double> (
            (*(*propensities[table_entry])[x + 1])[y])
            / (*all_propensities[x + 1])[y])
            / (static_cast<double> (amino_count[table_entry])
            / total));
    else
        enX = -log((static_cast<double> (
            (*(*propensities[table_entry])[x - 1])[y])
            / (*all_propensities[x - 1])[y])
            / (static_cast<double> (amino_count[table_entry])
            / total));


    long double enY = 0.0;
    if (diffY >= 0)
        enY = -log((static_cast<double> (
            (*(*propensities[table_entry])[x])[y + 1])
            / (*all_propensities[x])[y + 1])
            / (static_cast<double> (amino_count[table_entry])
            / total));
    else
        enY = -log((static_cast<double> (
            (*(*propensities[table_entry])[x])[y - 1])
            / (*all_propensities[x])[y - 1])
            / (static_cast<double> (amino_count[table_entry])
            / total));

    long double enXY = 0.0;
    if ((diffX >= 0) && (diffY >= 0))
        enXY = -log((static_cast<double> (
            (*(*propensities[table_entry])[x + 1])[y + 1])
            / (*all_propensities[x + 1])[y + 1])
            / (static_cast<double> (amino_count[table_entry])
            / total));
    else if ((diffX >= 0) && (diffY < 0))
        enXY = -log((static_cast<double> (
            (*(*propensities[table_entry])[x + 1])[y - 1])
            / (*all_propensities[x + 1])[y - 1])
            / (static_cast<double> (amino_count[table_entry])
            / total));
    else if ((diffX < 0) && (diffY >= 0))
        enXY = -log((static_cast<double> (
            (*(*propensities[table_entry])[x - 1])[y + 1])
            / (*all_propensities[x - 1])[y + 1])
            / (static_cast<double> (amino_count[table_entry])
            / total));
    else
        enXY = -log((static_cast<double> (
            (*(*propensities[table_entry])[x - 1])[y - 1])
            / (*all_propensities[x - 1])[y - 1])
            / (static_cast<double> (amino_count[table_entry])
            / total));


    double d = sqrt((diffX * diffX) + (diffY * diffY)) / 1.4;

    return ( ((diffX * en1) + ((1 - diffX) * enX))
            + ((diffY * en1) + ((1 - diffY) * enY))
            + ((d * en1) + ((1 - d) * enXY))
            / 3);
}

/** 
 * @Description calculates the total energy for the amino acid, considering a specific phi an psi value
 * @param the amino acid reference (AminoAcir&), values for phi and psi (double)
 * @return corresponding energy value(long double)
 */
long double PhiPsi::calculateEnergy(AminoAcid& aa, double phi, double psi) {
    if ((phi > 181.00) && (psi > 181.00))
        calculateEnergy(aa);
    else {
        // current amino acid type:
        int table_entry = aminoAcidThreeLetterTranslator(aa.getType());

        // offsets for the array
        int x = sGetPropBin(phi);
        int y = sGetPropBin(psi);

        if ((x > SIZE_OF_TABLE) || (y > SIZE_OF_TABLE))
            return 0;
        else {

            return -log((static_cast<double> ((*(*propensities[table_entry])[x])[y])
                    / (*all_propensities[x])[y])
                    / (static_cast<double> (amino_count[table_entry])
                    / total));
        }
    }
    return 0;
}

/**
 * @Description  Calculates the psi-psi torsional energy of the amino acid whose (phi, psi) pair is 
 *    specified by 'diheds' and whose type is specified by code.
 * @param   pair of dihedral angles(AminoAcir&), the amino acid code(AminoAcidCode)
 * @return   corresponding value(long double)
 */
long double PhiPsi::calculateEnergy(AminoAcid& diheds, AminoAcidCode code) {

    int table_entry = code;

    // locate propensity bin for the (phi, psi) pair
    int x = sGetPropBin(diheds.getPhi(true));
    int y = sGetPropBin(diheds.getPsi(true));

    if ((x > SIZE_OF_TABLE) || (y > SIZE_OF_TABLE)) // why?
        return 0;
    else {
        return -log((static_cast<double> ((*(*propensities[table_entry])[x])[y])
                / (*all_propensities[x])[y])
                / (static_cast<double> (amino_count[table_entry])
                / total));
    }
}

/**
 * @Description  Calculates the propensity binding value for a specific angle
 * @param   angle value (double)
 * @return   the corresponding value( int)
 */
inline int PhiPsi::sGetPropBin(double p) {
    int x = (int) (p + 180) / ARC_STEP;
    if (x == SIZE_OF_TABLE)
        x -= 1; // x is exactly 180 degrees and boosts array
    return x;
}

/**
 * @Description  Calculates the differenciated propensity binding and the 
 * @param   angle value(double), reference to the calculated diff propensity binding(double&)
 * @return   propensity value( int)
 */
inline int PhiPsi::sGetPropBinDiff(double p, double& diff) {
    int x = (int) (p + 180) / ARC_STEP;
    diff = 1 - ((static_cast<int> (p + 180) / ARC_STEP) / 2);

    if (x == SIZE_OF_TABLE)
        x -= 1; // x is exactly 180 degrees and boosts array
    return x;
}

/**
 * @Description  Adds a propensity for a specific amino acid type
 * @param   code of the amion acid(int), coordinates (int,int,int), values to set (int, int)
 * @return    changes the object internally (void)  
 */
inline void PhiPsi::sAddProp(int code, int x, int y) {
    (*(*propensities[code])[x])[y]++;
    amino_count[code]++;
}

/**
 * @Description Sets the step for the arc  
 * @param   the step value(int)
 * @return    changes the object internally (void)  
 */
inline void PhiPsi::setArcStep(int n) {
    ARC_STEP = n;
    SIZE_OF_TABLE = 360 / ARC_STEP;
    pResetData();
    pConstructData();
}

/**
 * @Description  Calculates the Maximum propensities value of the amino acid
 * @param   amino acid index corresponding to the enum (int)
 * @return  corresponding value (double)
 */
double PhiPsi::pGetMaxPropensities(int amino) {
    double max = 0.00;

    for (int j = 0; j < SIZE_OF_TABLE; j++) {
        for (int k = 0; k < SIZE_OF_TABLE; k++) {
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

/**
 * @Description  Calculates the Minimum propensities value of the amino acid
 * @param   amino acid index corresponding to the enum (int)
 * @return  corresponding value (double)
 */
double PhiPsi::pGetMinPropensities(int amino) {
    double min = 999.99;

    for (int j = 0; j < SIZE_OF_TABLE; j++) {
        for (int k = 0; k < SIZE_OF_TABLE; k++) {
            int tmp = (*all_propensities[j])[k];
            int tmp2 = (*(*propensities[amino])[j])[k];
            double propensities = ((static_cast<double> (tmp2) / tmp)
                    / ((static_cast<double> (amino_count[amino])) / total));
            if (propensities < min) {
                min = propensities;
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
double PhiPsi::pReturnMaxPropensities(int amino) {
    return amino_max_propensities[amino];
}

/**
 * @Description  returns the minimum propensities for the corresponding amino
 * @param   amino index from the enum (int)
 * @return   the corresponding value(double)
 */
double PhiPsi::pReturnMinPropensities(int amino) {
    return amino_min_propensities[amino];
}

/**
 * @Description fill an array containing max propensities for each amino
 * @param   none
 * @return    changes the object internally (void)  
 */
void PhiPsi::pConstructMaxPropensities() {
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_max_propensities.push_back(0);

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        amino_max_propensities[i] = pGetMaxPropensities(i);
    }
}

/**
 * @Description fill an array containing min propensities for each amino
 * @param   none
 * @return    changes the object internally (void)  
 */
void PhiPsi::pConstructMinPropensities() {
    for (int i = 0; i < AminoAcid_CODE_SIZE; i++)
        amino_min_propensities.push_back(0);

    for (int i = 0; i < AminoAcid_CODE_SIZE; i++) {
        amino_min_propensities[i] = pGetMinPropensities(i);
    }
}

/**
 * @Description calculates the energy from the PhiPsi angles of an amino acid 
 * @param   code of the amino acid(AminoAcidCode), values for phi and psi(double,double)
 * @return  corresponding value (double)
 */
double PhiPsi::getEnergyFromPhiPsi(AminoAcidCode code, double phi, double psi) {
    int table_entry = code;

    // locate propensity bin for the (phi, psi) pair
    int x = sGetPropBin(phi);
    int y = sGetPropBin(psi);

    if ((x > SIZE_OF_TABLE) || (y > SIZE_OF_TABLE)) // why?
        return 0;
    else {
        return -log((static_cast<double> ((*(*propensities[table_entry])[x])[y])
                / (*all_propensities[x])[y])
                / (static_cast<double> (amino_count[table_entry])
                / total));
    }
}
