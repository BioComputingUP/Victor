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
// Description:     Calculate structural scores with info derived from secondary
//                  structure.
//
// -----------------x-----------------------------------------------------------

#include <Sec.h>

namespace Biopool {

    // CONSTRUCTORS:

    Sec::Sec(SubMatrix *subStr, AlignmentData *ad, double cSec) : Structure(subStr),
    sec1(ad->getSequence(3)), sec2(ad->getSequence(4)), cSec(cSec) {
    }

    Sec::Sec(const Sec &orig) : Structure(orig) {
        copy(orig);
    }

    Sec::~Sec() {
    }


    // OPERATORS:

    Sec&
            Sec::operator =(const Sec &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // PREDICATES:

    double
    Sec::scoringStr(int i, int j) {
        return cSec * subStr->score[sec1[i - 1]][sec2[j - 1]];
    }


    // MODIFIERS:

    void
    Sec::copy(const Sec &orig) {
        Structure::copy(orig);
        sec1 = orig.sec1;
        sec2 = orig.sec2;
        cSec = orig.cSec;
    }

    Sec*
    Sec::newCopy() {
        Sec *tmp = new Sec(*this);
        return tmp;
    }

    void
    Sec::reverse() {
        string tmp = "";
        for (unsigned int i = sec2.length(); i > 0; i--)
            tmp.push_back(sec2[i - 1]);
        sec2 = tmp;
    }

} // namespace
