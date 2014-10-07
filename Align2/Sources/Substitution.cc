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
// Description:     Base class for deriving substitution matrices.
//
// -----------------x-----------------------------------------------------------

#include <Substitution.h>

namespace Biopool {

    // CONSTRUCTORS:

    Substitution::Substitution() {
    }

    Substitution::Substitution(const Substitution &orig) {
        copy(orig);
    }

    Substitution::~Substitution() {
    }


    // OPERATORS:

    Substitution&
            Substitution::operator =(const Substitution &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }

    ostream&
    operator <<(ostream &os, const Substitution &object) {
        Substitution::pWriteDoubleVector(os, object.score);
        return os;
    }

    istream&
    operator >>(istream &is, Substitution &object) {
        Substitution::pReadDoubleVector(is, object.score);
        return is;
    }


    // MODIFIERS:

    void
    Substitution::copy(const Substitution &orig) {
        score.clear();
        score.reserve(orig.score.size());
        for (unsigned int i = 0; i < orig.score.size(); i++) {
            vector<int> tmp;
            tmp.reserve(orig.score[i].size());
            for (unsigned int j = 0; j < orig.score[i].size(); j++)
                tmp.push_back(orig.score[i][j]);
            score.push_back(tmp);
        }
    }

    void
    Substitution::buildscore(const string &residues,
            const vector< vector<int> > &residuescores) {
        // Allow lowercase and uppercase residues (ASCII code <= 127)
        vector<int> row128(128, 0); // not sure if this should be 127!

        for (unsigned int i = 0; i < 128; ++i)
            score.push_back(row128);

        for (unsigned int i = 0; i < residues.size(); i++) {
            char res1 = residues[i];

            for (unsigned int j = 0; j <= i; j++) {
                char res2 = residues[j];
                score[res1][res2] = score[res2][res1] =
                        score[res1][res2 + 32] = score[res2 + 32][res1] =
                        score[res1 + 32][res2] = score[res2][res1 + 32] =
                        score[res1 + 32][res2 + 32] = score[res2 + 32][res1 + 32] =
                        residuescores[i][j];
            }
        }
    }


    // HELPERS:

    void Substitution::pWriteDoubleVector(ostream &os, vector< vector<int> > data) {
        os << data.size() << endl;

        for (unsigned int i = 0; i < data.size(); i++) {
            os << data[i].size() << "   ";
            for (unsigned int j = 0; j < data[i].size(); j++)
                os << data[i][j] << " ";
            os << "\n";
        }

        os << "\n";
    }

    template<class T> void
    Substitution::pReadDoubleVector(istream &is, vector< vector<T> > &data) {
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
    Substitution::pReadDoubleVector(istream &is, vector< vector<int> > &data);


    template void
    Substitution::pReadDoubleVector(istream &is, vector< vector<double> > &data);

} // namespace
