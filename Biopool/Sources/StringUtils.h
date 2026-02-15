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

/*
 * File:   CIFFileManager.h
 * Author: Francesco Menniti
 *
 * Created on 14 dicembre 2015
 */

#ifndef STRINGUTILS_H
#define STRINGUTILS_H

#include<iostream>
#include<stdio.h>
#include<string.h>
#include<vector>

using namespace std;

string Trim(string _input);
bool StartsWith(string _input, string _prefix);
void PrintOut(string line);
void PrintOut(int line);
vector<string> SplitOnSpace(string _input);

#endif /* STRINGUTILS_H */

