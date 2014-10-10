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
#include <PolarSolvationPotential.h>

using namespace Victor;

using namespace Victor::Biopool;

using namespace Victor::Energy;
// Global constants, typedefs, etc. (to avoid):
unsigned int PolarSolvationPotential::MAX_BINS = 15;
unsigned int PolarSolvationPotential::BIN_POLAR = 1;


// CONSTRUCTORS/DESTRUCTOR:

/**
 *@Description  Basic constructor, allocate the information from the solv.par file
 */
PolarSolvationPotential::PolarSolvationPotential() : sumPolar() {
    char *victor = getenv("VICTOR_ROOT");
    if (victor == NULL)
        ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
    string path = "data/solv.par"; //comment this if you are using a no 64bits OC
    string inputFile = getenv("VICTOR_ROOT");
    if (inputFile.length() < 3)
        ERROR("Environment variable VICTOR_ROOT was not found.", exception);
    //inputFile += SOLV_PARAM_FILE;//uncomment this if you are using a no 64bits OC
    inputFile += path; //comment this if you are using a no 64bits OC
    ifstream input(inputFile.c_str());

    if (!input)
        ERROR("Could not read data file.", exception);

    vector<int> tmp;
    for (unsigned int i = 0; i < MAX_BINS + 1; i++)
        tmp.push_back(0);

    vector<vector<int> > sum;
    for (unsigned int i = 0; i < AminoAcid_CODE_SIZE; i++)
        sum.push_back(tmp);

    for (unsigned int i = 0; i < (MAX_BINS / BIN_POLAR); i++)
        sumPolar.push_back(sum);

    while (input) {
        vector<unsigned int> num;
        string type;
        for (unsigned int i = 0; i < (MAX_BINS / BIN_POLAR); i++) {
            unsigned int tmp;
            input >> tmp;
            num.push_back(tmp);
        }
        input >> type;
        AminoAcidCode code = aminoAcidThreeLetterTranslator(type);

        for (unsigned int k = 0; k < (MAX_BINS / BIN_POLAR); k++)
            sumPolar[k][code][MAX_BINS] = num[k];

        for (unsigned int k = 0; k < (MAX_BINS / BIN_POLAR); k++)
            for (unsigned int i = 0; i < MAX_BINS; i++)
                input >> sumPolar[k][code][i];
    }

    input.close();
}

// PREDICATES:

/**
 *@Description Calculates the maximum energy for the amino acids in the spacer
 *@param  spacer reference(Spacer&)
 *@return  the corresponding value( double)
 */
long double PolarSolvationPotential::calculateSolvation(Spacer& sp) {
    long double solv = 0.0;
    for (unsigned int i = 0; i < sp.sizeAmino(); i++)
        solv += calculateSolvation(sp.getAmino(i), sp);
    return solv;
}

/**
 *@Description  Calculates energy the energy for a portion of the amino acids in the spacer. considering the portion from start(index 1) to en(Index2).
 *@param spacer reference(Spacer&), start and end of the Spacer portion (unsigned int, unsigned int)
 *@return  the sum of all the solvation of each amino acid in the  spacer portion(long double) 
 */
long double PolarSolvationPotential::calculateEnergy(Spacer& sp, unsigned int index1,
        unsigned int index2) {
    long double solv = 0.0;

    for (unsigned int i = index1; i < index2; i++)
        solv += calculateSolvation(sp.getAmino(i), sp);

    return solv;
}

/**
 *@Description Calculates solvation for a portion of amino acids in the spacer considering one specific amino acid
 *@param amino acid reference(AminoAcid&), spacer reference(Spacer&), start and end of the Spacer portion(unsigned int, unsigned int)
 *@return corresponding value of solvation (long double)
 */
long double PolarSolvationPotential::calculateSolvation(AminoAcid& aa, Spacer& sp,
        unsigned int start, unsigned int end) {
    if ((aa.getCode() == GLY) || (!aa.getSideChain().isMember(CB)))
        return 0.0;

    start = (start >= 0) ? start : 0;
    end = (end < sp.sizeAmino()) ? end : sp.sizeAmino() - 1;

    unsigned int count = 0;
    unsigned int countPolar = 0;
    for (unsigned int j = start; j < end; j++) {
        if (sp.getAmino(j).getSideChain().isMember(CB))
            if ((aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB])
                    <= SOLVATION_CUTOFF_DISTANCE_POLAR)
                    && (aa.getSideChain()[CB].distance(sp.getAmino(j).getSideChain()[CB]) > 0.0)) {// check not identical

                count++;
                if (isPolar(static_cast<AminoAcidCode> (sp.getAmino(j).getCode())))
                    countPolar++;
            }
    }

    unsigned int fracPolar = countPolar / BIN_POLAR;
    if (fracPolar >= BIN_POLAR)
        fracPolar = BIN_POLAR - 1;

    if (count >= MAX_BINS)
        count = MAX_BINS - 1;

    long double maxCount = 0.0;
    for (unsigned int k = 0; k < (MAX_BINS / BIN_POLAR); k++)
        for (unsigned int j = 0; j < sumPolar[k][aa.getCode()].size() - 1; j++)
            maxCount += static_cast<long double> (sumPolar[k][aa.getCode()][j]);

    long double maxCode = 0.0;
    for (unsigned int i = 0; i < sumPolar[fracPolar].size() - 1; i++)
        maxCode += static_cast<long double> (sumPolar[fracPolar][i][count]);

    long double maxAll;
    maxAll = 0.0;
    for (unsigned int k = 0; k < (MAX_BINS / BIN_POLAR); k++)
        for (unsigned int i = 0; i < sumPolar[k].size() - 1; i++)
            for (unsigned int j = 0; j < sumPolar[k][i].size() - 1; j++)
                maxAll += static_cast<long double> (sumPolar[k][i][j]);



    long double a = static_cast<long double> (
            sumPolar[fracPolar][aa.getCode()][count] + 1) /
            static_cast<long double> (maxCount);
    long double b = static_cast<long double> (maxCode + 1) / maxAll;


    return -0.582 * log(a / b);
}

