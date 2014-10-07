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
// -*-C++-*-------------------------------------------------------------------
//
//
// descriptipton: changes string to integer, double, float or long
// changes integer, double, float or long to string 
//
// error: if the input is no number (e.g. char or string)  then error,
// if wrong input (e.g. integer expected but float given) then error
//
// ---------------------------------------------------------------------------

#ifndef __STRING_2_NUMBER_H__
#define __STRING_2_NUMBER_H__

#include <iostream>
#include <sstream>
#include <string> // stl string class
#include <vector>
#include <Debug.h>
#include <limits>
#include <cmath>
// changes string into integer:
int stoi(const string&);
// changes string into unsigned integer:
unsigned int stoui(const string&);
// changes string into long:
long stol(const string&);
// changes string into float:
float stof(const string&);
// changes string into double:
double stod(const string&);
// changes string to vector of integer. Format: n1,n2,n3, ...
vector<int> sToVectorOfInt(const string&);
// changes string to vector of unsigned integer. Format: n1,n2,n3, ...
vector<unsigned int> sToVectorOfUInt(const string&);

// changes integer into string:
string itos(const int&);
// changes unsigned integer into string:
string uitos(const unsigned int&);
// changes long into string;
string ltos(const long&);
// changes float into string:
string ftos(const float&);
// changes double into string:
string dtos(const double&);
/** tokenize text line */
vector<string> getTokens(const string& s);

/** replace each occurence of c1 by c2 */
string translate(const string& s, char oldChar, char newChar);

/** return vector with positions, at which character c is found in string */
vector<unsigned int> findPositions(const string& s, char c);

#endif // __STRING_2_NUMBER_H__


