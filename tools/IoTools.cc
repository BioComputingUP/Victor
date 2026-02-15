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
/**
 * @Description of auxiliary I/O functions 
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctype.h>
#include <IoTools.h>

/**
 * @Description   Function tries to read and return an int from the 
 *  same input line. (as delimited by '\n')
 *  Precondition:  Stream must be opened.
 * Postcondition:  Returns the int if it existed or 0 otherwise.
 */
template<class T> int readOnSameLine(istream& is, T& result) {
    char c;
    do {
        is.get(c);
        // number found:
        if (isdigit(c) || (c == '-')) {
            is.putback(c);
            is >> result;
            return 1;
        }
    } while ((!is.eof()) && (c != '\n') && ((c == ' ') || (c == '\t')));
    // no number on line:
    return 0;
}

template int readOnSameLine<unsigned int>(istream& is, unsigned int& result);
template int readOnSameLine<double>(istream& is, double& result);

/**
 * @Description:  Function checks if next line is blank. 
 *  (ie. whitespace only)
 *
 *  Precondition:  Stream must be opened.
 *
 *  Postcondition:  Returns 1 if the next line is blank, 0 otherwise.
 *
 */

int checkForBlankLine(istream& is, bool returnbuffer) {
    char c;
    // already at end of file?
    if (is.eof())
        return 0;
    do {
        is.get(c);
        // line is blank:
        if (c == '\n')
            return 1;
    } while ((!is.eof()) && ((c == ' ') || (c == '\t')));
    // line is not blank:
    if (returnbuffer)
        is.putback(c);
    return 0;
}

/**
 * @Description:  checks if keyword is present at current position in stream
 * */
int checkForKeyword(istream& is, string key) {
    eatComment(is);
    char c;
    // already at end of file?
    if (is.eof())
        return 0;
    is.get(c);
    // found key?:
    if (c == key[0]) {
        string tmp;
        is >> tmp;
        if (key == (c + tmp))
            return 1;
    } else
        is.putback(c);
    return 0;
}

/**
 * @Description:  reads until newline is found*/
void skipToNewLine(istream& is) {
    char c;
    do {
        is.get(c);
    } while ((!is.eof()) && (c != '\n'));
}

/**
 * @Description:  skips over whitespace*/
void eatWhite(istream& is) {
    char c;
    while (is.get(c)) {
        if (!isspace(c)) {
            is.putback(c);
            break;
        }
    }
}

/**
 * @Description:  skips over comments (starting with '#')*/
void eatComment(istream& is) {
    char c;
    do {
        eatWhite(is);
        is.get(c);
        if (c != '#') {
            is.putback(c);
            break;
        }
        skipToNewLine(is); // search next new line
    } while (!is.eof());
}

/**
 * @Description:  peek if number follows, if so return it (1) else not (0)*/
template<class T> int readNumber(istream& is, T& result) {
    result = 0;
    eatComment(is);
    char c;
    is.get(c);
    // number found:
    if (isdigit(c) || (c == '-')) {
        is.putback(c);
        is >> result;
        return 1;
    }
    // no number on line:
    is.putback(c);
    return 0;
}

template int readNumber<int>(istream& is, int& result);

/**
 * @Description:  read a whole line into a string*/
string readLine(istream& is) {
    return readLine(is,true);
}

/**
 * @Description:  read a whole line into a string*/
string readLine(istream& is, bool _eatComment){
    string str = "";
    if(_eatComment)
    eatComment(is);
    char c;
    while (is.get(c)) {
        if ((c == '\n') || (c == '\r'))
            break;
        else
            str = str + c;
    }
    return str;
}

/**
 * @Description:  writes a vgVector3 in one statement*/
template<class T> ostream& operator <<(ostream &os, const vgVector3<T> &vec) {
    unsigned old_prec = os.precision();
    ios::fmtflags old_flags = os.flags();
    os.setf(ios::fixed, ios::floatfield);

    os << setw(8) << setprecision(3) << vec[0] << "   "
            << setw(8) << setprecision(3) << vec[1] << "   "
            << setw(8) << setprecision(3) << vec[2];

    os.precision(old_prec);
    os.flags(old_flags);
    return os;
}

template ostream& operator <<<double>(ostream &os, const vgVector3<double> &vec);

/**
 * @Description: reads a vgVector3 in one statement*/
template<class T> istream& operator >>(istream &is, vgVector3<T> &vec) {
    is >> vec[0] >> vec[1] >> vec[2];
    return is;
}

template istream& operator >><double>(istream &is, vgVector3<double> &vec);

/**
 * @Description:  writes a vgVector3 in one statement*/
template<class T> ostream& operator <<(ostream &os, const vgMatrix3<T> &vec) {
    unsigned old_prec = os.precision();
    ios::fmtflags old_flags = os.flags();
    os.setf(ios::fixed, ios::floatfield);

    os << setw(8) << setprecision(3) << vec[0] << "   "
            << setw(8) << setprecision(3) << vec[1] << "   "
            << setw(8) << setprecision(3) << vec[2];
    os << setw(8) << setprecision(3) << vec[3] << "   "
            << setw(8) << setprecision(3) << vec[4] << "   "
            << setw(8) << setprecision(3) << vec[5];
    os << setw(8) << setprecision(3) << vec[6] << "   "
            << setw(8) << setprecision(3) << vec[7] << "   "
            << setw(8) << setprecision(3) << vec[8];

    os.precision(old_prec);
    os.flags(old_flags);
    return os;
}

template ostream& operator <<<double>(ostream &os, const vgMatrix3<double> &vec);

/**
 * @Description:  reads a vgVector3 in one statement*/
template<class T> istream& operator >>(istream &is, vgMatrix3<T> &vec) {
    is >> vec[0] >> vec[1] >> vec[2];
    is >> vec[3] >> vec[4] >> vec[5];
    is >> vec[6] >> vec[7] >> vec[8];
    return is;
}

template istream& operator >><double>(istream &is, vgMatrix3<double> &vec);
