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
// Description:     Implement Smith-Waterman local alignment.
//
// -----------------x-----------------------------------------------------------

#include <SWAlign.h>
#include <limits.h>

namespace Victor { namespace Align2{

    // CONSTRUCTORS:
    /**
     * @description
     * @param ad
     * @param gf
     * @param ss
     */
    SWAlign::SWAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss)
    : Align(ad, gf, ss) {
        pCalculateMatrix(true);
    }

    SWAlign::SWAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss,
            const vector<unsigned int> &v1, const vector<unsigned int> &v2)
    : Align(ad, gf, ss) {
        pCalculateMatrix(v1, v2, true);
    }

    SWAlign::SWAlign(const SWAlign &orig) : Align(orig) {
    }

    SWAlign::~SWAlign() {
    }


    // OPERATORS:

    SWAlign&
            SWAlign::operator =(const SWAlign &orig) {
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
    SWAlign::getMultiMatch() {
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

    void
    SWAlign::copy(const SWAlign &orig) {
        Align::copy(orig);
    }
    /**
     * @description
     * @return 
     */
    SWAlign*
    SWAlign::newCopy() {
        SWAlign *tmp = new SWAlign(*this);
        return tmp;
    }


    // HELPERS:
    /**
     * @description
     * @param update
     */
    void
    SWAlign::pCalculateMatrix(bool update) {
        int maxi = n;
        int maxj = m;
        double maxval = INT_MIN;

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
                double val = max(max(max(z, extI), extJ), 0.00);

                if (update)
                    F[i][j] = val;

                if (EQUALS(val, 0))
                    B[i][j] = Traceback::getInvalidTraceback();
                else
                    if (val > 0) {
                    if (EQUALS(val, z))
                        B[i][j] = Traceback(i - 1, j - 1);
                    else
                        if (EQUALS(val, extJ))
                        B[i][j] = Traceback(i, j - 1);
                    else
                        if (EQUALS(val, extI))
                        B[i][j] = Traceback(i - 1, j);
                    else
                        ERROR("Error in SWAlign: SW 1", exception);
                } else
                    ERROR("Error in SWAlign: SW 2", exception);

                if (val > maxval) {
                    maxval = val;
                    maxi = i;
                    maxj = j;
                }

                B0 = Traceback(maxi, maxj);
            }
    }


    // SSEA variant
    /**
     * @description
     * @param v1
     * @param v2
     * @param update
     */
    void
    SWAlign::pCalculateMatrix(const vector<unsigned int> &v1,
            const vector<unsigned int> &v2, bool update) {
        // start SSEA variant code
        PRECOND((v1.size() == sq1.size()) && (v2.size() == sq2.size()), exception);
        unsigned int minL = 0;
        // end SSEA variant code

        int maxi = n;
        int maxj = m;
        double maxval = INT_MIN;

        for (int i = 1; i <= static_cast<int> (n); i++)
            for (int j = 1; j <= static_cast<int> (m); j++) {
                // start SSEA variant code
                if (v1[i - 1] < v2[j - 1])
                    minL = v1[i - 1];
                else
                    minL = v2[j - 1];

                double s = ss->scoring(i, j) * minL;
                // end SSEA variant

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
                double val = max(max(max(z, extI), extJ), 0.00);

                if (update)
                    F[i][j] = val;

                if (EQUALS(val, 0))
                    B[i][j] = Traceback::getInvalidTraceback();
                else
                    if (val > 0) {
                    if (EQUALS(val, z))
                        B[i][j] = Traceback(i - 1, j - 1);
                    else
                        if (EQUALS(val, extJ))
                        B[i][j] = Traceback(i, j - 1);
                    else
                        if (EQUALS(val, extI))
                        B[i][j] = Traceback(i - 1, j);
                    else
                        ERROR("Error in SWAlign: SW 1", exception);
                } else
                    ERROR("Error in SWAlign: SW 2", exception);

                if (val > maxval) {
                    maxval = val;
                    maxi = i;
                    maxj = j;
                }

                B0 = Traceback(maxi, maxj);
            }
    }

}} // namespace
