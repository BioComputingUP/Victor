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
//                  Panchenko method. Some explanations can be found in:
//
//                  Panchenko AR.
//                  Finding weak similarities between proteins by sequence
//                  profile comparison.
//                  Nucleic Acids Res. 2003 Jan 15;31(2):683-9.
//
// -----------------x-----------------------------------------------------------

#include <Panchenko.h>

namespace Biopool {

    // CONSTRUCTORS:

    Panchenko::Panchenko(Profile *pro1, Profile *pro2, PssmInput *pssm1,
            PssmInput *pssm2) : ScoringFunction(), pro1(pro1), pro2(pro2), pssm1(pssm1),
    pssm2(pssm2) {
    }

    Panchenko::Panchenko(const Panchenko &orig) : ScoringFunction(orig) {
        copy(orig);
    }

    Panchenko::~Panchenko() {
    }


    // OPERATORS:

    Panchenko&
            Panchenko::operator =(const Panchenko &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:

    double
    Panchenko::scoringSeq(int i, int j) {
        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

        double freq1 = 0.00;
        double freq2 = 0.00;
        double odds1 = 0.00;
        double odds2 = 0.00;
        double s1 = 0.00;
        double s2 = 0.00;

        int ni = Panchenko::returnAaColumnTarget(i);
        int nj = Panchenko::returnAaColumnTemplate(j);

        for (unsigned int k = 0; k < 20; k++) {
            freq1 = pro1->getAminoFrequency(residue_indices[k], (i - 1));
            freq2 = pro2->getAminoFrequency(residue_indices[k], (j - 1));

            odds1 = pssm1->score((i - 1), k);
            odds2 = pssm2->score((j - 1), k);

            s1 += (freq1 * odds2);
            s2 += (freq2 * odds1);
        }

        return ((ni * s1) + (nj * s2)) / (ni + nj);
    }

    int
    Panchenko::returnAaColumnTarget(int i) {
        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

        int aaColumn = 0;
        double freq1 = 0.00;

        for (unsigned int k = 0; k < 20; k++) {
            freq1 = pro1->getAminoFrequency(residue_indices[k], (i - 1));
            if (freq1 < 0.00001)
                continue;
            aaColumn++;
        }

        return aaColumn;
    }

    int
    Panchenko::returnAaColumnTemplate(int i) {
        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

        int aaColumn = 0;
        double freq2 = 0.00;

        for (unsigned int k = 0; k < 20; k++) {
            freq2 = pro2->getAminoFrequency(residue_indices[k], (i - 1));
            if (freq2 < 0.00001)
                continue;
            aaColumn++;
        }

        return aaColumn;
    }


    // MODIFIERS:

    void
    Panchenko::copy(const Panchenko &orig) {
        ScoringFunction::copy(orig);
        pro1 = orig.pro1->newCopy();
        pro2 = orig.pro2->newCopy();
        pssm1 = orig.pssm1->newCopy();
        pssm2 = orig.pssm2->newCopy();
    }

    Panchenko*
    Panchenko::newCopy() {
        Panchenko *tmp = new Panchenko(*this);
        return tmp;
    }

} // namespace
