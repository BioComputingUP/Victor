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
// Description:     Calculate structural scores with info derived from
//                  threading.
//
// -----------------x-----------------------------------------------------------

#include <Threading.h>

namespace Biopool {

    // CONSTRUCTORS:

    Threading::Threading(AlignmentData *ad, ThreadingInput *thread, double cThr)
    : Structure(0), seq1(ad->getSequence(1)), thread(thread), cThr(cThr) {
    }

    Threading::Threading(const Threading &orig) : Structure(orig) {
        copy(orig);
    }

    Threading::~Threading() {
    }


    // OPERATORS:
    /**
     * @description
     * @param orig
     * @return 
     */
    Threading&
            Threading::operator =(const Threading &orig) {
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
    Threading::scoringStr(int i, int j) {
        const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

        char aaTarget = seq1[i - 1];
        int targetIndex = 0;

        for (int k = 0; k < 20; k++)
            if (aaTarget == residue_indices[k]) {
                targetIndex = k;
                break;
            }

        return cThr * thread->score(targetIndex, (j - 1));
    }


    // MODIFIERS:
    /**
     * @description
     * @param orig
     */
    void
    Threading::copy(const Threading &orig) {
        Structure::copy(orig);
        seq1 = orig.seq1;
        thread = orig.thread->newCopy();
        cThr = orig.cThr;
    }

    Threading*
    Threading::newCopy() {
        Threading *tmp = new Threading(*this);
        return tmp;
    }

} // namespace
