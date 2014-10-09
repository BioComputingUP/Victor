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
// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:     Calculate scores for profile to profile alignment using
//                  sum of pairs method. Some explanations can be found in:
//
//                  Guoli Wang, Roland L. Dunbrack jr.
//                  Scoring profile-to-profile sequence alignments.
//                  Institute for Cancer Research, Fox Chase Cancer Center,
//                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
//
// -----------------x-----------------------------------------------------------

/*************************************************
 * sub is BLOSUM62 matrix standard log-odds form *
 *************************************************/

#include <LogAverage.h>

namespace Biopool {

    // CONSTRUCTORS:
    /**
     * @description
     * @param sub
     * @param pro1
     * @param pro2
     */
    LogAverage::LogAverage(SubMatrix *sub, Profile *pro1, Profile *pro2)
    : ScoringFunction(), sub(sub), pro1(pro1), pro2(pro2) {
    }

    LogAverage::LogAverage(const LogAverage &orig) : ScoringFunction(orig) {
        copy(orig);
    }

    LogAverage::~LogAverage() {
    }


    // OPERATORS:
    /**
     * @description
     * @param orig
     * @return 
     */
    LogAverage&
            LogAverage::operator =(const LogAverage &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:
    /**
     * @description
     * @param i
     * @param j
     * @return 
     */
    double
    LogAverage::scoringSeq(int i, int j) {
        double s = 0.00;

        for (AminoAcidCode amino1 = ALA; amino1 <= TYR; amino1++) {
            double tmp = 0.00;

            for (AminoAcidCode amino2 = ALA; amino2 <= TYR; amino2++)
                tmp += exp(sub->score[aminoAcidOneLetterTranslator(amino1)]
                    [aminoAcidOneLetterTranslator(amino2)]) *
                pro2->getAminoFrequencyFromCode(amino2, (j - 1));
            ostringstream convert; // stream used for the conversion
            convert << amino1;

            s += (pro1->getAminoFrequencyFromCode(amino1, (i - 1)) * tmp);
        }


        return log(s);
    }


    // MODIFIERS:
    /**
     * @description
     * @param orig
     */
    void
    LogAverage::copy(const LogAverage &orig) {
        ScoringFunction::copy(orig);
        sub = orig.sub->newCopy();
        pro1 = orig.pro1->newCopy();
        pro2 = orig.pro2->newCopy();
    }
    /**
     * @description
     * @return 
     */
    LogAverage*
    LogAverage::newCopy() {
        LogAverage *tmp = new LogAverage(*this);
        return tmp;
    }

} // namespace
