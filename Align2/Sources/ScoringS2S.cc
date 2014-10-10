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
// Description:     Calculate scores for sequence to sequence alignment.
//
// -----------------x-----------------------------------------------------------

#include <ScoringS2S.h>

namespace Victor { namespace Align2{

    // CONSTRUCTORS:
    /**
     * 
     * @param sub
     * @param ad
     * @param str
     * @param cSeq
     */
    ScoringS2S::ScoringS2S(SubMatrix *sub, AlignmentData *ad, Structure *str,
            double cSeq) : ScoringScheme(sub, ad, str), seq1(ad->getSequence(1)),
    seq2(ad->getSequence(2)), cSeq(cSeq) {
    }
    /**
     * 
     * @param orig
     */
    ScoringS2S::ScoringS2S(const ScoringS2S &orig) : ScoringScheme(orig) {
        copy(orig);
    }

    ScoringS2S::~ScoringS2S() {
    }


    // OPERATORS:

    ScoringS2S&
            ScoringS2S::operator =(const ScoringS2S &orig) {
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
    ScoringS2S::scoring(int i, int j) {
        double s = cSeq * sub->score[seq1[i - 1]][seq2[j - 1]];

        if (str != 0)
            s += str->scoringStr(i, j);

        return s;
    }


    // MODIFIERS:
    /**
     * 
     * @param orig
     */
    void
    ScoringS2S::copy(const ScoringS2S &orig) {
        ScoringScheme::copy(orig);
        seq1 = orig.seq1;
        seq2 = orig.seq2;
        cSeq = orig.cSeq;
    }

    ScoringS2S*
    ScoringS2S::newCopy() {
        ScoringS2S *tmp = new ScoringS2S(*this);
        return tmp;
    }
    /**
     * 
     */
    void
    ScoringS2S::reverse() {
        ScoringScheme::reverse();

        string tmp = "";
        for (unsigned int i = seq2.length(); i > 0; i--)
            tmp.push_back(seq2[i - 1]);
        seq2 = tmp;
    }

}} // namespace
