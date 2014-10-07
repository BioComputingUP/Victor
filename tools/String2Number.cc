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
#include "String2Number.h"

/**
 * @Description string to integer
 * @param string to convert
 */
int stoi(const string &s) {
    int i;
    double temp;
    const char *pointer;
    pointer = s.c_str();
    istringstream ist(pointer);
    ist >> temp;
    if (ist == NULL) {
        cerr << "\"" << s << "\" is no integer!!" << endl;
        ERROR("Integer expected", exception);
    }
    if ((temp - static_cast<int> (temp)) != 0) {
        cerr << "\"" << s << "\" is no integer!!" << endl;
        ERROR("Wrong format", exception);
    }
    i = static_cast<int> (temp);
    return i;
}

/** 
 * @name of function sToVectorOfUInt
 * @Description   if the input is no number (e.g. char or string)  then error, 
 *  if wrong input (integer expected but e.g. float given) then error
 * changes string into integer:
 * */
unsigned int stoui(const string &s) {
    unsigned int i;
    double temp;
    const char *pointer;
    pointer = s.c_str();
    istringstream ist(pointer);
    ist >> temp;
    if (ist == NULL) {
        cerr << "\"" << s << "\" is no unsigned integer!!" << endl;
        ERROR("Integer expected", exception);
    }
    if ((temp - static_cast<unsigned int> (temp)) != 0) {
        cerr << "\"" << s << "\" is no unsigned integer!!" << endl;
        ERROR("Wrong format", exception);
    }
    if (temp < 0) {
        ERROR("Negative unsigned integer", exception);
    }
    i = static_cast<unsigned int> (temp);
    return i;
}

/**
 * @Description string to long
 * @param string to convert
 */
long stol(const string &s) {
    long l;
    long double temp;
    const char *pointer;
    pointer = s.c_str();
    istringstream ist(pointer);
    ist >> temp;
    if (ist == NULL) {
        cerr << "\"" << s << "\" is no long integer!!" << endl;
        ERROR("Long expected", exception);
    }
    if ((temp - static_cast<long> (temp)) != 0) {
        cerr << "\"" << s << "\" is no long integer!!" << endl;
        ERROR("Wrong format", exception);
    }
    l = static_cast<long> (temp);
    return l;
}

/**
 * @Description string to float
 * @param string to convert
 */
float stof(const string &s) {
    float f;
    float temp;
    const char *pointer;
    pointer = s.c_str();
    istringstream ist(pointer);
    ist >> temp;
    if (ist == NULL) {
        cerr << "\"" << s << "\" is no float!!" << endl;
        ERROR("Float expected", exception);
    }
    if ((temp - static_cast<float> (temp)) != 0) {
        cerr << "\"" << s << "\" is no float!!" << endl;
        ERROR("Wrong format", exception);
    }
    f = temp;
    return f;
}

/**
 * @Description string to double
 * @param string to convert
 */
double stod(const string &s) {
    double d;
    double temp; // temporary integer
    const char *pointer;
    string temp_s;
    if (s == "e-99") {
        temp_s = "1e-99";
        pointer = temp_s.c_str();
    }
    else pointer = s.c_str();

    istringstream ist(pointer);

    ist >> temp;

    if (ist == NULL) {
        cerr << "\"" << s << "\" is no double!!" << endl;
        ERROR("Wrong format", exception);
    }
    d = temp;
    return d;
}

/** 
 * @name of function sToVectorOfInt
 * @Description 
 * changes string to vector of integer. Format: n1,n2,n3, ...
 * precondition: spaces between numbers. String can also be enclosed in "..."
 * postcondition: the return vector is filled with the apropriate numbers
 * error behaviour: quits, if precondition is not fullfilled.
 */
vector<int> sToVectorOfInt(const string& s_orig) {
    string s = s_orig;
    vector<int> v;
    if (s[1] == '\'') {
        DEBUG_MSG("Removing colons!");
        ASSERT(s[s.size() - 1] == '\'', exception);
#ifdef G_COMPILER
        s.remove(0U, 1U); // remove the first letter
        s.remove(s.size() - 1U, 1U); // remove last letter
#else
        s.erase(0, 1); // remove the first letter
        s.erase(s.size() - 1, 1);
#endif
    }
    istringstream ist(s.c_str()); // see Stroustrup chapter 21.5.3
    string w; // word
    while (ist >> w) {
        v.push_back(stoi(w));
        DUMP(w);
    }
    return v;
}

/** 
 * @name of function sToVectorOfUInt
 * @Description  changes string to vector of unsigned integer. Format: n1,n2,n3, ...
 * precondition: spaces between numbers. String can also be enclosed in "..."
 * postcondition: the return vector is filled with the apropriate numbers
 * error behaviour: quits, if precondition is not fullfilled.
 */
vector<unsigned int> sToVectorOfUInt(const string& s_orig) {
    string s = s_orig;
    vector<unsigned int> v;
    if (s[1] == '\'') {
        DEBUG_MSG("Removing colons!");
        ASSERT(s[s.size() - 1] == '\'', exception);
#ifdef G_COMPILER
        s.remove(0U, 1U); // remove the first letter
        s.remove(s.size() - 1U, 1U); // remove last letter
#else
        s.erase(0, 1); // remove the first letter
        s.erase(s.size() - 1, 1);
#endif
    }
    istringstream ist(s.c_str()); // see Stroustrup chapter 21.5.3
    string w; // word
    while (ist >> w) {
        v.push_back(stoi(w));
        DUMP(w);
    }
    return v;
}

/**
 * @Description changes integer into sting:
 */

string itos(const int &i) {
    string s;

    ostringstream os;
    os << i << '\0';
    return s = os.str();
}

/**
 * @Description changes unsigned integer into sting:
 * */
string uitos(const unsigned int &i) {
    string s;

    ostringstream os;
    os << i << '\0';
    return s = os.str();
}

/**
 * @Description changes long into sting:
 * */
string ltos(const long &l) {
    string s;

    ostringstream os;
    os << l << '\0';
    return s = os.str();
}

/**
 * @Description  changes float into sting:
 * */
string ftos(const float &f) {
    string s;

    ostringstream os;
    os << f << '\0';
    return s = os.str();
}

/**
 * @Description  changes double into sting:
 * */
string dtos(const double &d) {
    string s;

    ostringstream os;
    os << d; // << '\0';
    return s = os.str();
}

/**
 * @Description  tokenize text 
 * */
vector<string> getTokens(const string& text) {
    istringstream ist(text.c_str());
    char* charLine = new char[text.size() + 1]; // size of string
    vector<string> v;
    string s;
    while (!ist.eof()) {
        ist >> charLine;
        s = charLine; // assignment of c-strings to string!
        //    DUMP(s);
        if (s != "") {
            v.push_back(s);
        }
    }
    delete[] charLine;
    return v;
}

/**
 * @Description replace each occurence of c1 by c2 
 */
string
translate(const string& s, char c1, char c2) {
    string t = s;
    for (unsigned int i = 0; i < t.size(); ++i) {
        if (t[i] == c1) {
            t[i] = c2;
        }
    }
    return t;
}

/**
 * @Description return vector with positions, at which character c is found in string 
 */
vector<unsigned int> findPositions(const string& s, char c) {
    vector<unsigned int> result;
    for (unsigned int i = 0; i < s.size(); ++i) {
        if (s[i] == c) {
            result.push_back(i);
        }
    }
    return result;
}
