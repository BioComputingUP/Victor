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
// $Id: ThreadingSs2.cc,v 1.1 2008-05-21 09:45:27 biocomp Exp $
//
// Class:           ThreadingSs2

//
// Project name:    Victor/Align
//
// Date:            12/2007
//
// Description:     Calculate structural scores with info derived from
//                  threading and PSI-PRED.
//
// -----------------x-----------------------------------------------------------

#include <ThreadingSs2.h>

namespace Victor { namespace Align2{

    // CONSTRUCTORS:
    /**
     * 
     * @param subStr
     * @param ad
     * @param thread
     * @param psipred1
     * @param psipred2
     * @param cThr
     * @param cSs2
     */
    ThreadingSs2::ThreadingSs2(SubMatrix *subStr, AlignmentData *ad,
            ThreadingInput *thread, Ss2Input *psipred1, Ss2Input *psipred2, double cThr,
            double cSs2) : Structure(subStr), seq1(ad->getSequence(1)),
    sec1(ad->getSequence(3)), sec2(ad->getSequence(4)), thread(thread),
    psipred1(psipred1), psipred2(psipred2), cThr(cThr), cSs2(cSs2) {
    }

    ThreadingSs2::ThreadingSs2(const ThreadingSs2 &orig) : Structure(orig) {
        copy(orig);
    }

    ThreadingSs2::~ThreadingSs2() {
    }


    // OPERATORS:

    ThreadingSs2&
            ThreadingSs2::operator =(const ThreadingSs2 &orig) {
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
    ThreadingSs2::scoringStr(int i, int j) {
        //
        // THREADING
        //

        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

        char aaTarget = seq1[i - 1];
        int targetIndex = 0;

        for (int k = 0; k <= 19; k++)
            if (aaTarget == residue_indices[k]) {
                targetIndex = k;
                break;
            }

        double s1 = thread->score(targetIndex, (j - 1));


        //
        // SS2
        //

        const string ss2_indices = "HEC";

        // Target PSI-PRED
        double weigthH = psipred1->score((i - 1), 1);
        double weigthE = psipred1->score((i - 1), 2);
        double weigthC = psipred1->score((i - 1), 0);
        int substiH = subStr->score[sec2[j - 1]][ss2_indices[0]];
        int substiE = subStr->score[sec2[j - 1]][ss2_indices[1]];
        int substiC = subStr->score[sec2[j - 1]][ss2_indices[2]];
        double tmp1 = weigthH * substiH + weigthE * substiE + weigthC * substiC;

        // Template PSI-PRED
        weigthH = psipred2->score((j - 1), 1);
        weigthE = psipred2->score((j - 1), 2);
        weigthC = psipred2->score((j - 1), 0);
        substiH = subStr->score[sec1[i - 1]][ss2_indices[0]];
        substiE = subStr->score[sec1[i - 1]][ss2_indices[1]];
        substiC = subStr->score[sec1[i - 1]][ss2_indices[2]];
        double tmp2 = weigthH * substiH + weigthE * substiE + weigthC * substiC;

        double s2 = (tmp1 + tmp2) / 2;


        return cThr * s1 + cSs2 * s2;
    }


    // MODIFIERS:
    /**
     * 
     * @param orig
     */
    void
    ThreadingSs2::copy(const ThreadingSs2 &orig) {
        Structure::copy(orig);
        seq1 = orig.seq1;
        sec1 = orig.sec1;
        sec2 = orig.sec2;
        thread = orig.thread->newCopy();
        psipred1 = orig.psipred1->newCopy();
        psipred2 = orig.psipred2->newCopy();
        cThr = orig.cThr;
        cSs2 = orig.cSs2;
    }

    ThreadingSs2*
    ThreadingSs2::newCopy() {
        ThreadingSs2 *tmp = new ThreadingSs2(*this);
        return tmp;
    }

}} // namespace
