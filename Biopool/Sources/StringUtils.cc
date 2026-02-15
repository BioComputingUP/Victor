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

#include<StringUtils.h>

using namespace std;

string Trim(string _input)
{
    size_t first = _input.find_first_not_of(' ');
    size_t last = _input.find_last_not_of(' ');
    return _input.substr(first, (last-first+1));
}

bool StartsWith(string _input, string _prefix)
{
    return _input.compare(0, _prefix.length(), _prefix) == 0;
}

void PrintOut(string line)
{
    cout << line << endl;
}

void PrintOut(int line)
{
    cout << line << endl;
}

vector<string> SplitOnSpace(string _input) {
    
    vector<string>* result = new vector<string>();
    _input = Trim(_input);

    while (_input.size() > 0) {

        _input = Trim(_input);
        size_t first = _input.find_first_of(' ');

        if (first > 0 && first <= _input.size() ) {
            result->push_back(_input.substr(0, first + 1));
            _input = _input.substr(first + 1);
        } else {
            result->push_back(_input);
            _input = "";
        }
    }

    return *result;
}