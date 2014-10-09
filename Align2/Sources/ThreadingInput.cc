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

// Description:     Implement I/O objects for handling threading files.
//
// -----------------x-----------------------------------------------------------

#include <ThreadingInput.h>

namespace Biopool {

    // CONSTRUCTORS:

    ThreadingInput::ThreadingInput() {
    }

    ThreadingInput::ThreadingInput(istream &is) {
        is >> *this;
    }

    ThreadingInput::ThreadingInput(const ThreadingInput &orig) {
        copy(orig);
    }

    ThreadingInput::~ThreadingInput() {
    }


    // OPERATORS:
    /**
     * @description
     * @param orig
     * @return 
     */
    ThreadingInput&
            ThreadingInput::operator =(const ThreadingInput &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }

    ostream&
    operator <<(ostream &os, const ThreadingInput &object) {
        ThreadingInput::pWriteDoubleVector(os, object.residuescores);
        return os;
    }

    istream&
    operator >>(istream &is, ThreadingInput &object) {
        ThreadingInput::pReadDoubleVector(is, object.residuescores);
        return is;
    }


    // MODIFIERS:
    /**
     * @description
     * @param orig
     */
    void
    ThreadingInput::copy(const ThreadingInput &orig) {
        residuescores.clear();
        residuescores.reserve(orig.residuescores.size());
        for (unsigned int i = 0; i < orig.residuescores.size(); i++) {
            vector<double> tmp;
            tmp.reserve(orig.residuescores[i].size());
            for (unsigned int j = 0; j < orig.residuescores[i].size(); j++)
                tmp.push_back(orig.residuescores[i][j]);
            residuescores.push_back(tmp);
        }
    }
    /**
     * @description
     * @return 
     */
    ThreadingInput*
    ThreadingInput::newCopy() {
        ThreadingInput *tmp = new ThreadingInput(*this);
        return tmp;
    }


    // HELPERS:

    template<class T> void
    ThreadingInput::pWriteDoubleVector(ostream &os, vector< vector<T> > data) {
        os << data.size() << "   " << endl;

        for (unsigned int i = 0; i < data.size(); i++) {
            os << data[i].size() << "   ";
            for (unsigned int j = 0; j < data[j].size(); j++)
                os << data[i][j] << " ";
            os << endl;
        }
        os << endl;
    }
    /**
     * @description
     * @param is
     * @param data
     */
    template<class T> void
    ThreadingInput::pReadDoubleVector(istream &is, vector< vector<T> > &data) {
        unsigned int size1;
        is >> size1;

        for (unsigned int i = 0; i < size1; i++) {
            unsigned int size2;
            is >> size2;
            vector<T> row;
            row.reserve(size2);
            for (unsigned int j = 0; j < size2; j++) {
                T tmp;
                is >> tmp;
                row.push_back(tmp);
            }
            data.push_back(row);
        }
    }


    template void
    ThreadingInput::pReadDoubleVector(istream &is, vector< vector<int> > &data);


    template void
    ThreadingInput::pReadDoubleVector(istream &is, vector< vector<double> > &data);

} // namespace
