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
// -*- C++ -*---------------------------------------------------------------
//
// Description of "fileNamePath": get the path (if existing!!) of
// the file-name
//
// Precondition: file-name of which you get the path
//
// Postcondition: path of the file-name
//
//
// Description of "fileNameEnding": get the ending (if existing!!!)
// of the file-name
//
// Precondition: file-name of which you get the ending
//
// Postcondition: ending of the file-name
//
//
// Description of "fileNameBase": get the name-base (that means,
// the name without path and ending) of the file-name
//
// Precondition: file-name
//
// Postcondition: file-name without path and ending
//
// Description of "fileName": test, if a file with the given
// file-name already exists, if it exists make new file-name by adding
// an integer at the end (before the ending!!!!)
//
// Precondition: file-name which should be tested
//
// Postcondition: file-name which doesn't exists so far
//
// -------------------------------------------------------------------------


#ifndef __FILE_NAME_H__
#define __FILE_NAME_H__

// #define EG_COMPILER // define egcs compiler
// #define G_COMPILER // define g++ compiler

#include <iostream>
#include <fstream>
#include <string> // STL string class
#include "String2Number.h" // convert number to string and and vice versa

// input: filename:
string fileNamePath(const string&);
// input: filename:
string fileNameEnding(const string&);
// input: filename:
string fileNameBase(const string&);
// input: filename:
string fileName(const string&);

#endif // __FILE_NAME_H__


