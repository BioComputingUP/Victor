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
//                  Zhou-Zhou method. Some explanations can be found in:
//
//                  Zhou H., Zhou Y.
//                  Single-body residue-level knowledge-based energy score
//                  combined with sequence-profile and secondary structure
//                  information for fold recognition.
//                  Proteins. 2004 Jun 1;55(4):1005-13.
//                  PMID: 15146497 [PubMed - indexed for MEDLINE]
//
// -----------------x-----------------------------------------------------------

#include <Zhou.h>

namespace Victor { namespace Align2{

    // CONSTRUCTORS:
    /**
     * 
     * @param pro1
     * @param pssm2
     */
    Zhou::Zhou(Profile *pro1, PssmInput *pssm2) : ScoringFunction(), pro1(pro1),
    pssm2(pssm2) {
    }

    Zhou::Zhou(const Zhou &orig) : ScoringFunction(orig) {
        copy(orig);
    }

    Zhou::~Zhou() {
    }


    // OPERATORS:

    Zhou&
            Zhou::operator =(const Zhou &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:
    /**
     * 
     * @param i
     * @param j
     * @return 
     */
    double
    Zhou::scoringSeq(int i, int j) {
        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

        double freq1 = 0.00;
        double odds2 = 0.00;
        double s = 0.00;

        for (unsigned int k = 0; k < 20; k++) {
            freq1 = pro1->getAminoFrequency(residue_indices[k], (i - 1));
            odds2 = pssm2->score((j - 1), k);

            s += (freq1 * odds2);
        }

        return s;
    }


    // MODIFIERS:
    /**
     * 
     * @param orig
     */
    void
    Zhou::copy(const Zhou &orig) {
        ScoringFunction::copy(orig);
        pro1 = orig.pro1->newCopy();
        pssm2 = orig.pssm2->newCopy();
    }
    /**
     * 
     * @return 
     */
    Zhou*
    Zhou::newCopy() {
        Zhou *tmp = new Zhou(*this);
        return tmp;
    }

}} // namespace
