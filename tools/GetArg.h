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
* @Name of the method: getArg(s)
* @Description: 
* - interprets the inputs of the command line, assigne a
*   variable a command line input, possibly set default value,
* - check, if an option is chosen, if stop=true, the program will
*   end, if the option is chosen, if the text is defined, it will be
*   shown
* - if bool justOption=true, then you can just chose the option and
*   as input will be taken the default value, if the option is not
*   chosen, then the value will be ""
* - if you define "PATH", then the default name is added at the input,
*   if the input is just a path
* - if bool twiceSameNumber=false, the program will remove every
*   number, that is already member of the vector
* Call method by: getArg( "a", a_name, argc, argv [, default name] )
*                 getArg( "-all", all_name, argc, argv )    
* if an option has several characters, it must begin with a "-"
* (called at the command line with --all)
* 
* ERROR: 
* - if an opton is called several times or a value is specified
*   several times (except for the method "getArgs") 
* - if an option is called but no value is given (except if the bool
*   justOption is set on true)
* This can be hold back if "CONTINUE" is defined (return false)
*
*
* Name of the method: showText
* @Decription: show the text, if the option is chosen and ends, 
* if the option "CONTINUE" is defined, then the program doesn't end
* Call method by: showText( "t", text, argc, argv ) 
*                 showText( "-text", text, argc, argv )
*
*
* Name of the method: checkArg
* 
* Description: check, if there are chosen wrong options
* Call method by: checkArg( [text], argc, argv, 5, 
*                           "a", "b", "c", "-test", "-help" ) 
*
* ERROR: if a wrong option is chosen or the first input isn't an
* option
* This can be hold back if "CONTINUE" is defined (return false)
*/
#ifndef __GET_ARG_H__
#define __GET_ARG_H__
#include <iostream>
#include <cstdarg>
#include <string>
#include <vector>
#include <Debug.h>
#include <String2Number.h> // convert number to string and and vice versa
#include <FileName.h>
#include <limits>
// "value" has only one "word" as input, set default value (if
// wanted), set gotten path AND default value (if wanted):
bool getArg( const char *option, string &value, 
	     const int &argc, char *argv[], string defaultValue = "",
	     bool justOption = false );
bool getArg( const char *option, int &value, 
	     const int &argc, char *argv[], int defaultValue = 0 );
bool getArg( const char *option, unsigned int &value, 
	     const int &argc, char *argv[], unsigned int defaultValue = 0 );
bool getArg( const char *option, long &value, 
	     const int &argc, char *argv[], long defaultValue  = 0);
bool getArg( const char *option, float &value, 
	     const int &argc, char *argv[], float defaultValue = 0 );
bool getArg( const char *option, double &value, 
	     const int &argc, char *argv[], double defaultValue = 0 );
bool getArg( const char *option, vector<int> &value, 
	     const int &argc, char *argv[], 
	     bool twiceSameNumber = true );
bool getArg( const char *option, vector<unsigned int> &value, 
	     const int &argc, char *argv[], 
	     bool twiceSameNumber = true );
bool getArg( const char *option, vector<double> &value, 
	     const int &argc, char *argv[],
	     bool twiceSameNumber = true );
bool getArg( const char *option, vector<string> &value, 
	     const int &argc, char *argv[],
	     bool twiceSameNumber = true );

// several words with the same option are put together in the string
// "value":
bool getArgs( const char *option, string &value, 
	      const int &argc, char *argv[], string defaultValue = "" );

// check, if an option is chosen:
bool getArg( const char *givenOption, const int &argc, char *argv[],
	     string text = "", bool stop = false );

// check command line, if all options are right:
bool checkArg( const int &argc, char *argv[], int optionNumber ... );
bool checkArg( string text, const int &argc, char *argv[], 
	       int optionNumber ... );

void getOption(double& param, char* optName, int nArgs, char* argv[], 
	       bool verbose = false);

void getOption(unsigned int& param, char* optName, int nArgs, char* argv[], 
	       bool verbose = false);

void getOption(float& param, char* optName, int nArgs, char* argv[], 
	       bool verbose = false);
bool getOption(char* optName, int nArgs, char* argv[], bool verbose = true);

#endif // __GET_ARG_H__
