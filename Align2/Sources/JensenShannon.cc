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
// Description:     Calculate scores for profile to profile alignment using
//                  Jensen-Shannon method. Some explanations can be found in:
//
//                  Guoli Wang, Roland L. Dunbrack jr.
//                  Scoring profile-to-profile sequence alignments.
//                  Institute for Cancer Research, Fox Chase Cancer Center,
//                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
//
// -----------------x-----------------------------------------------------------

#include <JensenShannon.h>

namespace Biopool {

    // CONSTRUCTORS:

    JensenShannon::JensenShannon(Profile *pro1, Profile *pro2) : ScoringFunction(),
    pro1(pro1), pro2(pro2) {
        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

        for (unsigned int k = 0; k < 20; k++) {
            p1[k] = 0.00001;
            for (unsigned int n = 0; n < pro1->getSequenceLength(); n++)
                p1[k] += (pro1->getAminoFrequency(residue_indices[k], n) / pro1->getSequenceLength());
        }

        for (unsigned int k = 0; k < 20; k++) {
            p2[k] = 0.00001;
            for (unsigned int n = 0; n < pro2->getSequenceLength(); n++)
                p2[k] += (pro2->getAminoFrequency(residue_indices[k], n) / pro2->getSequenceLength());
        }
    }

    JensenShannon::JensenShannon(const JensenShannon &orig) : ScoringFunction(orig) {
        copy(orig);
    }

    JensenShannon::~JensenShannon() {
    }


    // OPERATORS:

    JensenShannon&
            JensenShannon::operator =(const JensenShannon &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:

    double
    JensenShannon::scoringSeq(int i, int j) {
        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

        double freq1 = 0.00;
        double freq2 = 0.00;
        double tmp1 = 0.00;
        double tmp2 = 0.00;
        double s1 = 0.00;
        double s2 = 0.00;
        double s3 = 0.00;
        double s4 = 0.00;

        for (unsigned int k = 0; k < 20; k++) {
            freq1 = pro1->getAminoFrequency(residue_indices[k], (i - 1)) + 0.00001;
            freq2 = pro2->getAminoFrequency(residue_indices[k], (j - 1)) + 0.00001;

            tmp1 = (freq1 + freq2) / 2;
            tmp2 = (tmp1 + ((p1[k] + p2[k]) / 2)) / 2;

            s1 += (freq1 * log2(freq1 / tmp1));
            s2 += (freq2 * log2(freq2 / tmp1));
            s3 += (tmp1 * log2(tmp1 / tmp2));
            s4 += (((p1[k] + p2[k]) / 2) * log2(((p1[k] + p2[k]) / 2) / tmp2));
        }

        double D = (s1 + s2) / 2;
        double S = (s3 + s4) / 2;

        return ((1 - D) * (1 + S)) / 2;
    }


    // MODIFIERS:

    void
    JensenShannon::copy(const JensenShannon &orig) {
        ScoringFunction::copy(orig);
        pro1 = orig.pro1->newCopy();
        pro2 = orig.pro2->newCopy();

        double p1[20];
        for (unsigned int k = 0; k < 20; k++) {
            p1[k] = orig.p1[k];
            cout << p1[k];
        }

        double p2[20];
        for (unsigned int k = 0; k < 20; k++) {
            p2[k] = orig.p2[k];
            cout << p2[k];
        }
    }

    JensenShannon*
    JensenShannon::newCopy() {
        JensenShannon *tmp = new JensenShannon(*this);
        return tmp;
    }

} // namespace
