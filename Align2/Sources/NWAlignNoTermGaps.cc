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
// Description:     Implement Needleman-Wunsh global alignment with no penalty
//                  for the terminal hangouts of the sequence.
//
// -----------------x-----------------------------------------------------------

#include <NWAlignNoTermGaps.h>

namespace Biopool {

    // CONSTRUCTORS:
    /**
     * @description
     * @param ad
     * @param gf
     * @param ss
     */
    NWAlignNoTermGaps::NWAlignNoTermGaps(AlignmentData *ad, GapFunction *gf,
            ScoringScheme *ss) : Align(ad, gf, ss) {
        pCalculateMatrix(true);
    }
    /**
     * @description
     * @param ad
     * @param gf
     * @param ss
     * @param v1
     * @param v2
     */
    NWAlignNoTermGaps::NWAlignNoTermGaps(AlignmentData *ad, GapFunction *gf,
            ScoringScheme *ss, const vector<unsigned int> &v1,
            const vector<unsigned int> &v2) : Align(ad, gf, ss) {
        pCalculateMatrix(v1, v2, true);
    }
    /**
     * @description
     * @param orig
     */
    NWAlignNoTermGaps::NWAlignNoTermGaps(const NWAlignNoTermGaps &orig)
    : Align(orig) {
    }

    NWAlignNoTermGaps::~NWAlignNoTermGaps() {
    }


    // OPERATORS:
    /**
     * @description
     * @param orig
     * @return 
     */
    NWAlignNoTermGaps&
            NWAlignNoTermGaps::operator =(const NWAlignNoTermGaps &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:
    /**
     * @description
     */
    void
    NWAlignNoTermGaps::getMultiMatch() {
        Traceback tb = B0;
        int i = tb.i;
        int j = tb.j;
        pModifyMatrix(i, j);

        tb = next(tb);
        while (((tb.i >= 0) || (tb.j >= 0)) && ((i != tb.i) || (j != tb.j))) {
            ad->calculateMatch(i, tb.i, j, tb.j);
            i = tb.i;
            j = tb.j;
            pModifyMatrix(i, j);
            tb = next(tb);
        }

        pCalculateMatrix(false); // recalculate F, B and B0 from modified data
        ad->getMatch();
    }


    // MODIFIERS:
    /**
     * @description
     * @param orig
     */
    void
    NWAlignNoTermGaps::copy(const NWAlignNoTermGaps &orig) {
        Align::copy(orig);
    }

    NWAlignNoTermGaps*
    NWAlignNoTermGaps::newCopy() {
        NWAlignNoTermGaps *tmp = new NWAlignNoTermGaps(*this);
        return tmp;
    }


    // HELPERS:
    /**
     * @description
     * @param update
     */
    void
    NWAlignNoTermGaps::pCalculateMatrix(bool update) {
        if (update)
            F[0][0] = 0;

        for (int i = 1; i <= static_cast<int> (n); i++) {
            if (update)
                F[i][0] = 0;
            B[i][0] = Traceback(i - 1, 0);
        }

        for (int j = 1; j <= static_cast<int> (m); j++) {
            if (update)
                F[0][j] = 0;
            B[0][j] = Traceback(0, j - 1);
        }

        for (int i = 1; i <= static_cast<int> (n); i++)
            for (int j = 1; j <= static_cast<int> (m); j++) {
                double s = ss->scoring(i, j);
                double extI, extJ;

                if ((i != 1) && (j != 1)) {
                    if (B[i - 1][j].j == j)
                        extI = F[i - 1][j] - gf->getExtensionPenalty(j);
                    else
                        if (B[i - 1][j].j == (j - 1))
                        extI = F[i - 1][j] - gf->getOpenPenalty(j);
                } else
                    extI = F[i - 1][j] - gf->getOpenPenalty(j);

                if ((i != 1) && (j != 1)) {
                    if (B[i][j - 1].i == i)
                        extJ = F[i][j - 1] - gf->getExtensionPenalty(j);
                    else
                        if (B[i][j - 1].i == (i - 1))
                        extJ = F[i][j - 1] - gf->getOpenPenalty(j);
                } else
                    extJ = F[i][j - 1] - gf->getOpenPenalty(j);

                double z = F[i - 1][j - 1] + s;
                double val = max(max(z, extI), extJ);

                if (update)
                    F[i][j] = val;

                if (EQUALS(val, z))
                    B[i][j] = Traceback(i - 1, j - 1);
                else
                    if (EQUALS(val, extJ))
                    B[i][j] = Traceback(i, j - 1);
                else
                    if (EQUALS(val, extI))
                    B[i][j] = Traceback(i - 1, j);
                else
                    ERROR("Error in NWAlignNoTermGaps: NW 1", exception);
            }

        B0 = Traceback(n, m);
    }


    // SSEA variant
    /**
     * @description
     * @param v1
     * @param v2
     * @param update
     */
    void
    NWAlignNoTermGaps::pCalculateMatrix(const vector<unsigned int> &v1,
            const vector<unsigned int> &v2, bool update) {
        // start SSEA variant code
        PRECOND((v1.size() == sq1.size()) && (v2.size() == sq2.size()), exception);
        unsigned int minL = 0;
        // end SSEA variant code

        if (update)
            F[0][0] = 0;

        for (int i = 1; i <= static_cast<int> (n); i++) {
            if (update)
                F[i][0] = 0;
            B[i][0] = Traceback(i - 1, 0);
        }

        for (int j = 1; j <= static_cast<int> (m); j++) {
            if (update)
                F[0][j] = 0;
            B[0][j] = Traceback(0, j - 1);
        }

        for (int i = 1; i <= static_cast<int> (n); i++)
            for (int j = 1; j <= static_cast<int> (m); j++) {
                // start SSEA variant code
                if (v1[i - 1] < v2[j - 1])
                    minL = v1[i - 1];
                else
                    minL = v2[j - 1];

                double s = ss->scoring(i, j) * minL;
                // end SSEA variant code

                double extI, extJ;

                if ((i != 1) && (j != 1)) {
                    if (B[i - 1][j].j == j)
                        extI = F[i - 1][j] - gf->getExtensionPenalty(j);
                    else
                        if (B[i - 1][j].j == (j - 1))
                        extI = F[i - 1][j] - gf->getOpenPenalty(j);
                } else
                    extI = F[i - 1][j] - gf->getOpenPenalty(j);

                if ((i != 1) && (j != 1)) {
                    if (B[i][j - 1].i == i)
                        extJ = F[i][j - 1] - gf->getExtensionPenalty(j);
                    else
                        if (B[i][j - 1].i == (i - 1))
                        extJ = F[i][j - 1] - gf->getOpenPenalty(j);
                } else
                    extJ = F[i][j - 1] - gf->getOpenPenalty(j);

                double z = F[i - 1][j - 1] + s;
                double val = max(max(z, extI), extJ);

                if (update)
                    F[i][j] = val;

                if (EQUALS(val, z))
                    B[i][j] = Traceback(i - 1, j - 1);
                else
                    if (EQUALS(val, extJ))
                    B[i][j] = Traceback(i, j - 1);
                else
                    if (EQUALS(val, extI))
                    B[i][j] = Traceback(i - 1, j);
                else
                    ERROR("Error in NWAlignNoTermGaps: NW 1", exception);
            }

        B0 = Traceback(n, m);
    }

} // namespace
