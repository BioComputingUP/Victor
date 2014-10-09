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
// Description:     Implement I/O objects for handling PHD files.
//
// -----------------x-----------------------------------------------------------

#include <ProfInput.h>

namespace Biopool {

    // CONSTRUCTORS:

    ProfInput::ProfInput() {
    }

    ProfInput::ProfInput(istream &is) {
        is >> *this;
    }

    ProfInput::ProfInput(const ProfInput &orig) {
        copy(orig);
    }

    ProfInput::~ProfInput() {
    }


    // OPERATORS:
    /**
     * @description
     * @param orig
     * @return 
     */
    ProfInput&
            ProfInput::operator =(const ProfInput &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }
    /**
     * @description
     * @param os
     * @param object
     * @return 
     */
    ostream&
    operator <<(ostream &os, const ProfInput &object) {
        ProfInput::pWriteString(os, object.seq, object.profSSPred,
                object.profBEPred);
        return os;
    }
    /**
     * @description
     * @param is
     * @param object
     * @return 
     */
    istream&
    operator >>(istream &is, ProfInput &object) {
        ProfInput::pReadString(is, object.seq, object.profSSPred,
                object.profBEPred);
        return is;
    }


    // PREDICATES:
    /**
     * @description
     * @param ssChar
     * @param beChar
     * @return 
     */
    char
    ProfInput::getProfMixSSBE(char ssChar, char beChar) {
        char correctChar;

        if (ssChar == 'H') // H case
        {
            correctChar = 'H';
            if (beChar == 'e')
                correctChar = 'I';
            return correctChar;
        } else
            if (ssChar == 'E') // E case
        {
            correctChar = 'E';
            if (beChar == 'e')
                correctChar = 'F';
            return correctChar;
        } else // C case
        {
            correctChar = 'C';
            if (beChar == 'e')
                correctChar = 'D';
            return correctChar;
        }
    }


    // MODIFIERS:
    /**
     * @description
     * @param orig
     */
    void
    ProfInput::copy(const ProfInput &orig) {
        seq = orig.seq;
        profSSPred = orig.profSSPred;
        profBEPred = orig.profBEPred;
    }
    /**
     * @description
     * @return 
     */
    ProfInput*
    ProfInput::newCopy() {
        ProfInput *tmp = new ProfInput(*this);
        return tmp;
    }


    // HELPERS:
    /**
     * @description
     * @param os
     * @param data1
     * @param data2
     * @param data3
     */
    void
    ProfInput::pWriteString(ostream &os, string data1, string data2, string data3) {
        os << "     #    AA   Pss   Pbe\n" << endl;

        for (unsigned int i = 0; i < data1.size(); i++) {
            os << setw(6) << i
                    << setw(6) << data1[i]
                    << setw(6) << data2[i]
                    << setw(6) << data3[i] << endl;
        }
        os << endl;
    }
    /**
     * @description
     * @param is
     * @param data1
     * @param data2
     * @param data3
     */
    void
    ProfInput::pReadString(istream &is, string &data1, string &data2, string &data3) {
        string line;
        data1 = "";
        data2 = "";
        data3 = "";
        unsigned int lineType = 0;

        while ((is.good()) && (!is.eof()))
            switch (lineType) {
                case 0:
                    line = readLine(is);
                    lineType++;
                    break;
                case 1:
                    data1 += readLine(is);
                    lineType++;
                    break;
                case 2:
                    data2 += readLine(is);
                    lineType++;
                    break;
                case 3:
                    data3 += readLine(is);
                    lineType = 0;
            }
    }

} // namespace
