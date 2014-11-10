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
#include "GetArg.h"
#include <string.h>
#include <limits.h>

/**
 * @Name of the method: searchValueOfOption
 * @Description  evaluate an argv and define the input and the option according to argv
 *
 */
void searchValueOfOption(char option[80], string &input, char *argv) {
    if (argv[0] == '-') {
        if (argv[1] == '-') {
            int k = 0;
            while (argv[k] != '\0') {
                option[k] = argv[k + 1];
                k++;
            }
            input = "";
        }
        else {
            option[0] = argv[1];
            option[1] = '\0';
            if (argv[2] != '\0') {
                char tempName[80];

                int k = 1;
                while (argv[k] != '\0') {
                    tempName[k - 1] = argv[k + 1];
                    k++;
                }
                input = tempName;
            } else {
                input = "";
            }
        }
    }
    else {
        input = argv;
    }
}

/**
 * @Description:  interprets the inputs of the commando line, assigne a  variable a commando line input
 *
 * @Precondition: option of the variable
 *
 * @Postcondition: the value of the variable
 * the value can be
 * - a string
 * - a integer
 * - a long integer
 * - a float
 * - a double
 * - a vector of integers
 * - a vector of doubles
 * 
 * Errors: if an option is chosen several times, an argument is
 * defined severals times or if an option is chosen but no argument is
 * given
 */
bool getArg(const char *givenOption, string &value,
        const int &argc, char *argv[], string defaultValue,
        bool justOption) {
    int countOption = 0, countValue = 0;
    string input;
    char option[80];

    option[0] = '\0';
    for (int i = 1; i < argc; i++) {
        searchValueOfOption(option, input, argv[i]);
        if (argv[i][0] == '-') {
            if (strcmp(givenOption, option) == 0) {
                countOption++;
            } else {
                continue;
            }
        }
        if (input != "") {
            if (strcmp(givenOption, option) == 0) {
                countValue++;
                value = input;
            } else {
                continue;
            }
        }
    }


    if (countOption > 1) {
#ifndef CONTINUE
        cerr << "You have given too often the option \""
                << givenOption << "\"" << endl;
        ERROR("Too many options", exception);
#else
        return false;
#endif
    }
    if (countValue > 1) {
#ifndef CONTINUE
        cerr << "You have given too many input of the option \""
                << givenOption << "\"" << endl;
        ERROR("Too many options", exception);
#else
        return false;
#endif
    }
    if (justOption == false) {
        if ((countOption == 1) && (countValue == 0)) {
#ifndef CONTINUE
            cerr << "You have chosen the option \"" << givenOption
                    << "\" but havn't specified the value!" << endl;
            ERROR("Missing value", exception);
#else
            return false;
#endif
        }
        if (value == "") {
            value = defaultValue;
        }
#ifdef PATH
        else {

            if (value[value.size() - 1] == '/') {
                value = value + defaultValue;
            }
        }
#endif
    }
    else {
        if ((countOption == 1) && (countValue == 0) && (value == "")) {
            value = defaultValue;
        }
        if ((countOption == 0) && (countValue == 0)) {
            value = "";
        }
#ifdef PATH

        if (value[value.size() - 1] == '/') {
            value = value + defaultValue;
        }
#endif	  
    }
    return true;
}

bool getArg(const char *givenOption, int &value,
        const int &argc, char *argv[], int defaultValue) {
    bool result;
    string tempInt;

    value = defaultValue;
    result = getArg(givenOption, tempInt, argc, argv);
    if (tempInt != "") {
        value = stoiDEF(tempInt);
    }
    return result;
}

bool getArg(const char *givenOption, unsigned int &value,
        const int &argc, char *argv[], unsigned int defaultValue) {
    bool result;
    string tempInt;

    value = defaultValue;
    result = getArg(givenOption, tempInt, argc, argv);
    if (tempInt != "") {
        value = stouiDEF(tempInt);
    }
    return result;
}

bool getArg(const char *givenOption, long &value,
        const int &argc, char *argv[], long defaultValue) {
    bool result;
    string tempInt;

    value = defaultValue;
    result = getArg(givenOption, tempInt, argc, argv);
    if (tempInt != "") {
        value = stolDEF(tempInt);
    }
    return result;
}

bool getArg(const char *givenOption, float &value,
        const int &argc, char *argv[], float defaultValue) {
    bool result;
    string tempInt;

    value = defaultValue;
    result = getArg(givenOption, tempInt, argc, argv);
    if (tempInt != "") {
        value = stofDEF(tempInt);
    }
    return result;
}

bool getArg(const char *givenOption, double &value,
        const int &argc, char *argv[], double defaultValue) {
    bool result;
    string tempInt;

    value = defaultValue;
    result = getArg(givenOption, tempInt, argc, argv);
    if (tempInt != "") {
        value = stodDEF(tempInt);
    }
    return result;
}

bool getArg(const char *givenOption, vector<int> &value,
        const int &argc, char *argv[], bool twiceSameNumber) {
    int countOption = 0;
    string input;
    char option[80];

    option[0] = '\0';
    for (int i = 1; i < argc; i++) {
        searchValueOfOption(option, input, argv[i]);
        // if there's an option, get this option:
        if (argv[i][0] == '-') {
            if (strcmp(givenOption, option) == 0) {
                countOption++;
            } else {
                continue;
            }
        }
        if (input != "") {
            if (twiceSameNumber == true) {

                if (strcmp(givenOption, option) == 0) {
                    value.push_back(stoiDEF(input));
                } else {
                    continue;
                }
            } else {

                if (strcmp(givenOption, option) == 0) {
                    bool alreadyAdded = false;
                    for (unsigned int i = 0; i < value.size(); i++) {
                        if (stoiDEF(input) == value[i]) {
                            alreadyAdded = true;
                            break;
                        }
                    }
                    if (alreadyAdded == false) {
                        value.push_back(stoiDEF(input));
                    }
                }
                else {
                    continue;
                }
            }
        }
    }

    if ((countOption > 0) && (value.size() == 0)) {
#ifndef CONTINUE
        cerr << "You have chosen the option \""
                << givenOption << "\" but haven't specified the value!" << endl;
        ERROR("Missing value(s)", exception);
#else
        return false;
#endif
    }
    return true;
}

bool getArg(const char *givenOption, vector<unsigned int> &value,
        const int &argc, char *argv[], bool twiceSameNumber) {
    int countOption = 0;
    string input;
    char option[80];

    option[0] = '\0';
    for (int i = 1; i < argc; i++) {
        searchValueOfOption(option, input, argv[i]);

        if (argv[i][0] == '-') {
            if (strcmp(givenOption, option) == 0) {
                countOption++;
            } else {
                continue;
            }
        }

        if (input != "") {
            if (twiceSameNumber == true) {

                if (strcmp(givenOption, option) == 0) {
                    value.push_back(stouiDEF(input));
                } else {
                    continue;
                }
            }
            else {

                if (strcmp(givenOption, option) == 0) {
                    bool alreadyAdded = false;
                    for (unsigned int i = 0; i < value.size(); i++) {
                        if (stouiDEF(input) == value[i]) {
                            alreadyAdded = true;
                            break;
                        }
                    }
                    if (alreadyAdded == false) {
                        value.push_back(stouiDEF(input));
                    }
                }
                else {
                    continue;
                }
            }
        }
    }

    if ((countOption > 0) && (value.size() == 0)) {
#ifndef CONTINUE
        cerr << "You have chosen the option \""
                << givenOption << "\" but haven't specified the value!" << endl;
        ERROR("Missing value(s)", exception);
#else
        return false;
#endif
    }
    return true;
}

bool getArg(const char *givenOption, vector<double> &value,
        const int &argc, char *argv[], bool twiceSameNumber) {
    int countOption = 0;
    string input;
    char option[80];

    option[0] = '\0';
    for (int i = 1; i < argc; i++) {
        searchValueOfOption(option, input, argv[i]);

        if (argv[i][0] == '-') {
            if (strcmp(givenOption, option) == 0) {
                countOption++;
            } else {
                continue;
            }
        }

        if (input != "") {
            if (twiceSameNumber == true) {

                if (strcmp(givenOption, option) == 0) {
                    value.push_back(stodDEF(input));
                } else {
                    continue;
                }
            }
            else {

                if (strcmp(givenOption, option) == 0) {
                    bool alreadyAdded = false;
                    for (unsigned int i = 0; i < value.size(); i++) {
                        if (stodDEF(input) == value[i]) {
                            alreadyAdded = true;
                            break;
                        }
                    }
                    if (alreadyAdded == false) {
                        value.push_back(stodDEF(input));
                    }
                }
                else {
                    continue;
                }
            }
        }
    }

    if ((countOption > 0) && (value.size() == 0)) {
#ifndef CONTINUE
        cerr << "You have chosen the option \""
                << givenOption << "\" but haven't specified the value!" << endl;
        ERROR("Missing value(s)", exception);
#else
        return false;
#endif
    }
    return true;
}

bool getArg(const char *givenOption, vector<string> &value,
        const int &argc, char *argv[], bool twiceSameNumber) {
    int countOption = 0;
    string input;
    char option[80];

    option[0] = '\0';
    for (int i = 1; i < argc; i++) {
        searchValueOfOption(option, input, argv[i]);

        if (argv[i][0] == '-') {
            if (strcmp(givenOption, option) == 0) {
                countOption++;
            } else {
                continue;
            }
        }

        if (input != "") {
            if (twiceSameNumber == true) {

                if (strcmp(givenOption, option) == 0) {
                    value.push_back(input);
                } else {
                    continue;
                }
            }
            else {

                if (strcmp(givenOption, option) == 0) {
                    bool alreadyAdded = false;
                    for (unsigned int i = 0; i < value.size(); i++) {
                        if (input.compare(value[i]) == 0) {
                            alreadyAdded = true;
                            break;
                        }
                    }
                    if (alreadyAdded == false) {
                        value.push_back(input);
                    }
                }
                else {
                    continue;
                }
            }
        }
    }

    if ((countOption > 0) && (value.size() == 0)) {
#ifndef CONTINUE
        cerr << "You have chosen the option \""
                << givenOption << "\" but haven't specified the value!" << endl;
        ERROR("Missing value(s)", exception);
#else
        return false;
#endif
    }
    return true;
}

/**
 * @Name of the method: getArgs
 * @Description:interprets the inputs of the commando line, assigne a
 * variable a commando line input, several words will be put in one
 * string
 *
 * @Precondition: option of the variable
 *
 * @Postcondition: the value of the variable
 * 
 **/
bool getArgs(const char *givenOption, string &value,
        const int &argc, char *argv[], string defaultValue) {
    int countOption = 0;
    char option[80];
    string input;

    for (int i = 1; i < argc; i++) {
        searchValueOfOption(option, input, argv[i]);
        if (input != "") {
            if (strcmp(givenOption, option) == 0) {
                value = value + " " + input;
            } else {
                continue;
            }
        }
    }

    if ((countOption > 0) && (value.size() == 0)) {
#ifndef CONTINUE
        cerr << "You have chosen the option \""
                << givenOption << "\" but haven't specified the value!" << endl;
        ERROR("Missing value(s)", exception);
#else
        return false;
#endif
    }
    if (value.size() == 0) {
        value = defaultValue;
    }
    return true;
}

/**
 * @Name of the method: getArg
 * @Description: check, if an option is chosen
 */

bool getArg(const char *givenOption, const int &argc, char *argv[],
        string text, bool stop) {
    char option[80];
    string str;

    option[0] = '\0';
    for (int i = 1; i < argc; i++) {
        searchValueOfOption(option, str, argv[i]);
        if (strcmp(givenOption, option) == 0) {
            if (text != "") {
                cout << text << endl;
            }
            if (stop == true) {
                exit(0);
            }
            return true;
        }
    }
    return false;
}

/**
 * @Name of the method: checkArg
 * @Description:check, if there are chosen wrong options, if an
 * unknown option is chosen, then error, if you give a text, then show
 * this text
 */

bool checkArg(const int &argc, char *argv[], int optionNumber ...) {
    bool rightOption = false;
    char option[80];
    string str;
    va_list listOfOptions;

    option[0] = '\0';
    if (argc <= 1) {
        return true;
    }
    // check, if the first input is an option:
    if (argv[1][0] != '-') {
#ifndef CONTINUE 
        ERROR("First input isn't an option", exception);
#else
        return false;
#endif 
    }
    for (int i = 1; i < argc; i++) {
        rightOption = false;
        searchValueOfOption(option, str, argv[i]);
        // define the listOfOptions:
        va_start(listOfOptions, optionNumber);
        // every given Options is a *char, get all givenOptions:
        for (int i = 0; i < optionNumber; i++) {
            char *givenOption = va_arg(listOfOptions, char*);
            if (givenOption == 0) {
                break;
            }
            if (strcmp(givenOption, option) == 0) {
                rightOption = true;
                break; // go to the next option
            }
        }

        va_end(listOfOptions);
        if (rightOption == false) {
#ifndef CONTINUE 
            cerr << "\"" << option << "\" is a wrong option" << endl;
            ERROR("Wrong option", exception);
#else
            return false;
#endif 
        }
    }
    return true;
}

// show text:

bool checkArg(string text, const int &argc, char *argv[],
        int optionNumber ...) {
    bool rightOption = false;
    char option[80];
    string str;
    va_list listOfOptions;

    option[0] = '\0';
    if (argc <= 1) {
        return true;
    }
    // check, if the first input is an option:
    if (argv[1][0] != '-') {
#ifndef CONTINUE 
        cerr << text << endl;
        ERROR("First input isn't an option", exception);
#else
        return false;
#endif 
    }
    for (int i = 1; i < argc; i++) {
        rightOption = false;
        searchValueOfOption(option, str, argv[i]);
        // define the listOfOptions:
        va_start(listOfOptions, optionNumber);
        // every given Options is a *char, get all givenOptions:
        for (int i = 0; i < optionNumber; i++) {
            char *givenOption = va_arg(listOfOptions, char*);
            if (givenOption == 0) {
                break;
            }
            if (strcmp(givenOption, option) == 0) {
                rightOption = true;
                break; // go to the next option
            }
        }

        va_end(listOfOptions);
        if (rightOption == false) {
#ifndef CONTINUE 
            cerr << "\"" << option << "\" is a wrong option" << endl
                    << text << endl;
            ERROR("Wrong option", exception);
#else
            return false;
#endif 
        }
    }
    return true;
}

void getOption(double& param, char* optName, int nArgs, char* argv[],
        bool verbose) {
    double tmp = -100000.0;
    getArg(optName, tmp, nArgs, argv, -100000.0);
    if (tmp > -100000.0) {
        param = tmp;
        if (verbose)
            cout << optName << " = " << tmp << "\n";
    }
}

void getOption(unsigned int& param, char* optName, int nArgs, char* argv[],
        bool verbose) {
    unsigned int tmp = INT_MAX;
    getArg(optName, tmp, nArgs, argv, INT_MAX);
    if (tmp != INT_MAX) {
        param = tmp;
        if (verbose)
            cout << optName << " = " << tmp << "\n";
    }
}

void getOption(float& param, char* optName, int nArgs, char* argv[],
        bool verbose) {
    float tmp = -100000.0;
    getArg(optName, tmp, nArgs, argv, -100000.0);
    if (tmp > -100000.0) {
        param = tmp;
        if (verbose)
            cout << optName << " = " << tmp << "\n";
    }
}

bool getOption(char* optName, int nArgs, char* argv[],
        bool verbose) {
    bool tmp = getArg(optName, nArgs, argv);
    if ((tmp) && (verbose))
        cout << optName << "\n";
    return tmp;
}
