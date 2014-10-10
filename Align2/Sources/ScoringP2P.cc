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
// Description:     Calculate scores for profile to profile alignment.
//
// -----------------x-----------------------------------------------------------

#include <ScoringP2P.h>

namespace Victor { namespace Align2{

    // CONSTRUCTORS:
    /**
     * @description
     * @param sub
     * @param ad
     * @param str
     * @param pro1
     * @param pro2
     * @param fun
     * @param cSeq
     */
    ScoringP2P::ScoringP2P(SubMatrix *sub, AlignmentData *ad, Structure *str,
            Profile *pro1, Profile *pro2, ScoringFunction *fun, double cSeq)
    : ScoringScheme(sub, ad, str), seq1(ad->getSequence(1)),
    seq2(ad->getSequence(2)), pro1(pro1), pro2(pro2), fun(fun),
    cSeq(cSeq) {
    }

    ScoringP2P::ScoringP2P(const ScoringP2P &orig) : ScoringScheme(orig) {
        copy(orig);
    }

    ScoringP2P::~ScoringP2P() {
    }


    // OPERATORS:
    /**
     * @description
     * @param orig
     * @return 
     */
    ScoringP2P&
            ScoringP2P::operator =(const ScoringP2P &orig) {
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
    ScoringP2P::scoring(int i, int j) {
        double s = cSeq * fun->scoringSeq(i, j);
        if (str != 0)
            s += str->scoringStr(i, j);
        return s;
    }


    // MODIFIERS:
    /**
     * @description
     * @param orig
     */
    void
    ScoringP2P::copy(const ScoringP2P &orig) {
        ScoringScheme::copy(orig);
        seq1 = orig.seq1;
        seq2 = orig.seq2;
        pro1 = orig.pro1->newCopy();
        pro2 = orig.pro2->newCopy();
        fun = orig.fun->newCopy();
        cSeq = orig.cSeq;
    }
    /**
     * @description
     * @return 
     */
    ScoringP2P*
    ScoringP2P::newCopy() {
        ScoringP2P *tmp = new ScoringP2P(*this);
        return tmp;
    }

    void
    ScoringP2P::reverse() {
        ScoringScheme::reverse();

        string tmp = "";
        for (unsigned int i = seq2.length(); i > 0; i--)
            tmp.push_back(seq2[i - 1]);
        seq2 = tmp;

        pro2->reverse();
    }

}} // namespace
