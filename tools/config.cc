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
@Description */
#include <config.h>

string trenner = " \t,;"; // these chars disconnect 
// the name and the value of the parameter

vector<parameter> param;

/**
 * @Method:        Constructor
 * @Class:         config
 * @Project Name:  Lobo
 * @Description: Create a config object which includes all parameters from the
 *               config file "fileName" and it's includes.
 *               The default constructor is using the default name "config.cfg"
 *
 * */
config::config(string fileName) {
    param = scan_configFile(fileName);
}

config::config() {
    string fileName = "config.cfg";
    param = scan_configFile(fileName);
}

/**
 * @Method:        Constructor
 * @Class:         config

 * @Reviewed By:   -
 * @Description:  returns 1 if the parameter exists otherwise 0
 */
int
config::existParameter(string paramName) {
    unsigned int i = 0;
    while (i < param.size()) {
        if (param[i++].name.compare(paramName) == 0) {
            return 1;
        }
    }
    return 0;
}

/**
 * @Method:        Constructor
 * @Class:         config

 * @Reviewed By:   -
 * @Description:   Returns the value of the parameter "paramName".
 *                 The return value is from the same type as the parameter tmp.
 */
template<class Tget> Tget config::getParameter(string paramName, Tget tmp) {
    unsigned int i = 0;
    while ((i < param.size()) && (param[i].name.compare(paramName) != 0)) {
        i++;
    }
    if (i == param.size()) {
        ERROR("parameter not found", exception);
    } // Entry not found
    // Now extract the typeid of the string
    return convertString2Type(param[i].value, tmp);
}

template int config::getParameter<int>(string paramName, int value);
template unsigned int config::getParameter<unsigned int>(string paramName, unsigned int value);
template long config::getParameter<long>(string paramName, long value);
template float config::getParameter<float>(string paramName, float value);
template double config::getParameter<double>(string paramName, double value);
template string config::getParameter<string>(string paramName, string value);

int config::convertString2Type(string value, int tmp) {
    return stoiDEF(value);
}

unsigned int config::convertString2Type(string value, unsigned int tmp) {
    return stouiDEF(value);
}

double config::convertString2Type(string value, double tmp) {
    return stodDEF(value);
}

float config::convertString2Type(string value, float tmp) {
    return stofDEF(value);
}

long config::convertString2Type(string value, long tmp) {
    return stolDEF(value);
}

string config::convertString2Type(string value, string tmp) {
    return value;
}

/**
 * @Method:        Constructor
 * @Class:         config
 * @Reviewed By:   -
 * @Description:    Set a new value of the parameter. The new value is NOT permanent set
 *                 in the file. Only while the object exists.
 */
template <class Tset> int config::setParameter(string paramName, Tset value) {
    unsigned int i = 0;
    int found = 0;
    while ((i < param.size()) && !found) {
        if (param[i++].name.compare(paramName) == 0) {
            found = 1;
        }
    }
    if (!found) {
        return 0;
    } // Entry not found
    param[--i].value = convertTypes(value);
    return 1;
}

template int config::setParameter<int>(string paramName, int value);
template int config::setParameter<unsigned int>(string paramName, unsigned int);
template int config::setParameter<double>(string paramName, double);
template int config::setParameter<string>(string paramName, string);
template int config::setParameter<long>(string paramName, long);
template int config::setParameter<float>(string paramName, float);

string config::convertTypes(int value) {
    return itosDEF(value);
}

string config::convertTypes(unsigned int value) {
    return uitosDEF(value);
}

string config::convertTypes(float value) {
    return ftosDEF(value);
}

string config::convertTypes(string value) {
    return value;
}

string config::convertTypes(long value) {
    return ltosDEF(value);
}

string config::convertTypes(double value) {
    return dtosDEF(value);
}

template int config::changeParameter<int>(string paramName, int);
template int config::changeParameter<unsigned int>(string paramName, unsigned int);
template int config::changeParameter<double>(string paramName, double);
template int config::changeParameter<string>(string paramName, string);
template int config::changeParameter<long>(string paramName, long);
template int config::changeParameter<float>(string paramName, float);

/**
 * @Method:        Constructor
 * @Class:         config
 * @Description:  Set a new value. It changes the config file too. 
 *                 The new value is permanent changed !
 */

template<class Tchange> int config::changeParameter(string paramName, Tchange value) {
    unsigned int i = 0;
    string tmpline;
    string aktline;
    int finish = 0;
    if (!setParameter(paramName, value)) {
        return 0;
    } // first update the parameter at the heap
    while ((i < param.size()) && (param[i].name.compare(paramName) != 0)) {
        i++;
    }
    if (i == param.size()) {
        return 0;
    } // Entry not found 
    fstream datei(param[i].file.c_str(), ios::in | ios::out);
    datei.seekp(0, ios::beg); // go to the first line in the file
    while (getline(datei, tmpline) && !finish) {
        string::size_type anfpos = tmpline.find_first_not_of(" ");
        string::size_type pos = tmpline.find('#');
        string::size_type anf, end;
        if ((pos != anfpos) && !tmpline.empty()) // If no comment and not empty
        {
            // cout << "a value was found in a config file " << fileName << endl;
            anf = tmpline.find_first_not_of(trenner);
            end = tmpline.find_first_of(trenner, anf);
            if (tmpline.substr(anf, end - anf).compare(paramName) == 0) // substract the name
            {
                finish = 1;
                aktline = tmpline; //because tmpline will be lost in the while statement
            }
        }
    }
    if (finish) {
        datei.seekp(-(aktline.size() + tmpline.size() + 2), ios::cur);
        string delstr(aktline.size(), ' ');
        datei.write(delstr.c_str(), delstr.size()); // delete line
        datei.seekp(-delstr.size(), ios::cur);
        datei.write((param[i].name + " " + param[i].value).c_str(), param[i].name.size() + 1 + param[i].value.size());
        datei.close();
        return 1; // all successfully done
    } else {
        datei.close();
        return 0;
    }
}

template int config::newParameter<int>(string fileName, string paramName, int);
template int config::newParameter<unsigned int>(string fileName, string paramName, unsigned int);
template int config::newParameter<double>(string fileName, string paramName, double);
template int config::newParameter<string>(string fileName, string paramName, string);
template int config::newParameter<long>(string fileName, string paramName, long);
template int config::newParameter<float>(string fileName, string paramName, float);

/**
 * @Method:        Constructor
 * @Class:         config
 * @Description:   To create a new parameter. The parameter is also written in the file.
 */


template<class Tnew> int config::newParameter(string fileName, string paramName, Tnew value) {
    string tmpline;
    parameter tmpparam;
    if (existParameter(paramName)) {
        return 0;
    }
    tmpparam.name = paramName;
    tmpparam.file = fileName;
    param.push_back(tmpparam);
    if (!existParameter(paramName) || (!setParameter(paramName, value))) {
        return 0;
    } // well there is something wrong !
    fstream datei(fileName.c_str(), ios::app);
    if (!datei) {
        ERROR("Error while file handling ! The parameter was not created.",
                exception);
    }
    datei.put('\n'); // create a new line
    if (typeid (Tnew) == typeid (string)) {
        tmpline = param[param.size() - 1].name + " \"" + param[param.size() - 1].value + "\"";
    } else {
        tmpline = param[param.size() - 1].name + " " + param[param.size() - 1].value;
    }
    for (unsigned int wr = 0; wr < tmpline.size(); wr++) {
        datei.put(tmpline[wr]);
    }
    datei.put('\n');
    datei.close();
    return 1; // all successfully done
}

/**
 * @Method:        Constructor
 * @Class:         config
 -
 * @Description:  Delete a parameter. This change is also in the file done.
 */

int config::delParameter(string paramName) {
    long dpos;
    unsigned int count;
    unsigned int i = 0;
    string tmpline;
    string aktline;
    int finish = 0;
    int found = 0;

    while ((i < param.size()) && (!found)) {
        if (param[i++].name.compare(paramName) == 0) {
            found = 1;
        }
    }
    if (!found) {
        return 0;
    } // Entry not found 
    fstream datei(param[--i].file.c_str(), ios::in | ios::out);
    if (!datei) {
        return 0;
    }
    datei.seekp(0, ios::beg); // go to the first line in the file
    dpos = datei.tellg();

    count = 0;
    while ((!finish) && (getline(datei, tmpline))) {
        string::size_type anfpos = tmpline.find_first_not_of(' ');
        string::size_type pos = tmpline.find('#');
        string::size_type anf, end;
        if ((pos != anfpos) && !tmpline.empty()) // If no comment and not empty
        {
            anf = tmpline.find_first_not_of(trenner);
            end = tmpline.find_first_of(trenner, anf);
            if (tmpline.substr(anf, end - anf).compare(paramName) == 0) // substract the name
            {

                count = tmpline.size();
                finish = 1;
            } else {
                dpos = datei.tellg();
            }
        }
    }
    if (finish) {
        datei.seekp(dpos);
        datei.seekg(dpos);
        char tc;
        for (unsigned int wr = 0; wr < count; wr++) {
            datei.get(tc);
            datei.seekp(dpos + wr);
            datei.seekg(dpos + wr);

            datei.put(' '); // delete line
        }
        datei.close();
        // delete parameter in the existing object
        param[i] = param[param.size() - 1];
        param.pop_back();
        return 1; // all successfully done
    } else {
        datei.close();
        return 0;
    } // well there`s something wrong :(
}

/**
 * @Method:        Constructor
 * @Class:         config
 * @Description: Scans the file for parameters and prepares it for the other functions.
 *                 Normally you�ll not use these method.
 */

vector<parameter>
config::scan_configFile(string fileName) {
    vector<parameter> tmpparam;
    parameter tmp;
    vector<parameter> paramIncludeReturn;
    string tmpstr;

    ifstream is(fileName.c_str(), ios::in);
    if (!is) {
        ERROR("Could not open the config file !", exception);
    }
    while (getline(is, tmpstr)) {
        string::size_type pos = tmpstr.find("#include<");

        if (pos != string::npos) {
            string::size_type posend = tmpstr.find(">");
            if (posend == string::npos) {
                ERROR("Error in #include filename.", exception);
            }

            string tmpinclude(tmpstr, pos + 9, posend - (9 + pos));
            paramIncludeReturn = scan_configFile(tmpinclude);

            for (unsigned int i = 0; i < paramIncludeReturn.size(); i++) {
                tmpparam.push_back(paramIncludeReturn[i]);
            }
        } else {
            string::size_type anfpos = tmpstr.find_first_not_of(" ");
            string::size_type pos = tmpstr.find('#');
            string::size_type anf, end;
            if ((pos != anfpos) && !tmpstr.empty()) // If no comment and not empty
            {
                anf = tmpstr.find_first_not_of(trenner);
                end = tmpstr.find_first_of(trenner, anf);
                tmp.name = tmpstr.substr(anf, end - anf); // substract the name
                anf = tmpstr.find_first_not_of(trenner, end);
                if (tmpstr[anf] == '"') {
                    anf++;
                    end = tmpstr.find_first_of("\"", anf); // the value is a sting

                } else {
                    end = tmpstr.find_first_of(trenner, anf); // substract the normal value
                }
                if (anf == string::npos) {
                    ERROR("Parameter without value.", exception);
                }
                tmp.value = tmpstr.substr(anf, end - anf); // substract the value
                tmp.file = fileName;
                tmpparam.push_back(tmp);
            }
        }
    }
    is.close();
    return tmpparam;
}
