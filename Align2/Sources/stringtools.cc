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
// stringtools.cc

#include <stringtools.h>
#include <cstdlib>
using namespace std;

/**
 * @description
 * @param str
 */
void
strip(string &str) {
    string temp = "";

    for (string::iterator pos = str.begin(); pos != str.end(); pos++)
        if (*pos != ' ' and *pos != '\n')
            temp += *pos;

    str = temp;
}

/**
 * @description
 * @param str
 * @return 
 */
string
strip_and_return_string(string str) {
    string temp = "";

    for (string::iterator pos = str.begin(); pos != str.end(); pos++)
        if (*pos != ' ' and *pos != '\n')
            temp += *pos;

    return temp;
}
/**
 * @description
 * @param str
 * @return 
 */
int
strip_and_return_int(string str) {
    string temp = "";

    for (string::iterator pos = str.begin(); pos != str.end(); pos++)
        if (*pos != ' ' and *pos != '\n')
            temp += *pos;

    int integ = atoi(temp.c_str());
    return integ;
}
/**
 * @description
 * @param str
 * @return 
 */
unsigned int
strip_and_return_unsigned_int(string str) {
    string temp = "";

    for (string::iterator pos = str.begin(); pos != str.end(); pos++)
        if (*pos != ' ' and *pos != '\n')
            temp += *pos;

    unsigned int integ = atoi(temp.c_str());
    return integ;
}
/**
 * @description
 * @param str
 * @return 
 */
float
strip_and_return_float(string str) {
    string temp = "";

    for (string::iterator pos = str.begin(); pos != str.end(); pos++) {
        if (*pos == ',')
            temp += '.';
        else
            if (*pos != ' ' and *pos != '\n')
            temp += *pos;
    }

    float doub = atof(temp.c_str());
    return doub;
}
/**
 * @description
 * @param text_to_split
 * @param delimiter
 * @return 
 */
vector<string>
split(string text_to_split, char delimiter) {
    vector<string> vector_with_splited_text;
    string temp_text = "";

    for (string::iterator pos = text_to_split.begin(); pos != text_to_split.end(); pos++)
        if (*pos != delimiter and *pos != '\n')
            temp_text += *pos;
        else {
            vector_with_splited_text.push_back(temp_text);
            temp_text = "";
        }
    vector_with_splited_text.push_back(temp_text); // last piece

    return vector_with_splited_text;
}
/**
 * @description
 * @param str
 * @return 
 */
int
seq_length(const string str) {
    unsigned int size = 0;

    for (string::const_iterator pos = str.begin(); *pos != ' '; pos++)
        size++;

    return size;
}


/// Change each element of the string to upper case.
/**
 * @description
 * @param strToConvert
 * @return 
 */
string
StringToUpper(string strToConvert) {
    for (unsigned int i = 0; i < strToConvert.length(); i++)
        strToConvert[i] = toupper(strToConvert[i]);
    return strToConvert; // return the converted string
}


/// Change each element of the string to lower case.
/**
 * @description
 * @param strToConvert
 * @return 
 */
string
StringToLower(string strToConvert) {
    for (unsigned int i = 0; i < strToConvert.length(); i++)
        strToConvert[i] = tolower(strToConvert[i]);
    return strToConvert; // return the converted string
}
