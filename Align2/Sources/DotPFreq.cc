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
//                  dot product method. Some explanations can be found in:
//
// 	                Mittelman D., Sadreyev R., Grishin N.
//                  Probabilistic scoring measures for profile-profile
//                  comparison yield more accurate short seed alignments.
//                  Bioinformatics. 2003 Aug 12;19(12):1531-9.
//                  PMID: 12912834 [PubMed - indexed for MEDLINE]
//
//                  Marti-Renom MA., Madhusudhan MS., Sali A.
//                  Alignment of protein sequences by their profiles.
//                  Protein Sci. 2004 Apr;13(4):1071-87.
//                  PMID: 15044736 [PubMed - indexed for MEDLINE]
//
// -----------------x-----------------------------------------------------------

#include <DotPFreq.h>

namespace Victor { namespace Align2{

    // CONSTRUCTORS:

    /**
     *  
     * @param pro1
     * @param pro2
     */
    DotPFreq::DotPFreq(Profile *pro1, Profile *pro2) : ScoringFunction(),
    pro1(pro1), pro2(pro2) {
    }

    /**
     *  
     * @param orig
     */
    DotPFreq::DotPFreq(const DotPFreq &orig) : ScoringFunction(orig) {
        copy(orig);
    }

    /**
     *  
     */
    DotPFreq::~DotPFreq() {
    }


    // OPERATORS:

    /**
     *  
     * @param orig
     * @return 
     */
    DotPFreq&
            DotPFreq::operator =(const DotPFreq &orig) {
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
    double    DotPFreq::scoringSeq(int i, int j) {
        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

        double freq1 = 0.00;
        double freq2 = 0.00;
        double s = 0.00;

        for (unsigned int k = 0; k < 20; k++) {
            freq1 = pro1->getAminoFrequency(residue_indices[k], (i - 1));
            freq2 = pro2->getAminoFrequency(residue_indices[k], (j - 1));

            s += (freq1 * freq2);
        }

        return s;
    }


    // MODIFIERS:

    /**
     *  
     * @param orig
     */
    void DotPFreq::copy(const DotPFreq &orig) {
        ScoringFunction::copy(orig);
        pro1 = orig.pro1->newCopy();
        pro2 = orig.pro2->newCopy();
    }

    /**
     *  
     * @return 
     */
    DotPFreq*
    DotPFreq::newCopy() {
        DotPFreq *tmp = new DotPFreq(*this);
        return tmp;
    }

}} // namespace
