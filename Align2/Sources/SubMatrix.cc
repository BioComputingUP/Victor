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
// Description:     Implement a standard substitution matrix.
//
// -----------------x-----------------------------------------------------------

#include <SubMatrix.h>

namespace Biopool {

    // CONSTRUCTORS:

    SubMatrix::SubMatrix() : Substitution() {
    }

    SubMatrix::SubMatrix(istream &is) : Substitution() {
        is >> *this;
    }

    SubMatrix::SubMatrix(const SubMatrix &orig) {
        copy(orig);
    }

    SubMatrix::~SubMatrix() {
    }


    // OPERATORS:

    SubMatrix&
            SubMatrix::operator =(const SubMatrix &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }

    ostream&
    operator <<(ostream &os, const SubMatrix &object) {
        os << object.residues << endl;
        Substitution::pWriteDoubleVector(os, object.residuescores);
        return os;
    }

    istream&
    operator >>(istream &is, SubMatrix &object) {
        is >> object.residues;
        Substitution::pReadDoubleVector(is, object.residuescores);
        object.buildscore(object.residues, object.residuescores);
        return is;
    }


    // MODIFIERS:

    void
    SubMatrix::copy(const SubMatrix &orig) {
        Substitution::copy(orig);

        residuescores.clear();
        for (unsigned int n = 0; n < orig.residuescores.size(); n++)
            residuescores.push_back(orig.residuescores[n]);

        residues = orig.residues;
    }

    SubMatrix*
    SubMatrix::newCopy() {
        SubMatrix *tmp = new SubMatrix(*this);
        return tmp;
    }

} // namespace
