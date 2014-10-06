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
#include "FileName.h"

/**
* @Name of the method: fileNamePath
* @Description:get the path (if existing!!) of the file-name
*
* Precondition: file-name of which you get the path
*
* Postcondition: path of the file-name
*/
string fileNamePath( const string &fileName )
{
  string path, tempFileName;

  tempFileName = fileName;
  // first check, if there is a path, that means, if there is the
  // letter '/':
  if( tempFileName.find_last_of('/') < tempFileName.size() )
    {
      // remove everything but the path:
#ifdef G_COMPILER
      tempFileName.remove( (tempFileName.find_last_of('/')+1) ); 
#else
      tempFileName.erase( (tempFileName.find_last_of('/')+1) ); 
#endif
      path = tempFileName;
    }
  return path;
}

/**
* @Name of the method: fileNameEnding
* @Description: get the ending (if existing!!!) of the file-name
*
* Precondition: file-name of which you get the ending
*
* Postcondition: ending of the file-name
*/

string fileNameEnding( const string &fileName )
{
  string ending, tempFileName;

  tempFileName = fileName;
  // first check, if there is an ending, that means, if there is the
  // letter '.':
  if( tempFileName.find_last_of('.') < tempFileName.size() )
    {
      // now check, if the point is not part of the path, therefore
      // look first, IF there exists an path and then, If this path is
      // written with a point:
      if( (tempFileName.size() < tempFileName.find_last_of('/')) ||
	  (tempFileName.find_last_of('/') < tempFileName.find_last_of('.')) )
	{
#ifdef G_COMPILER
	  tempFileName.remove( 0, tempFileName.find_last_of('.') ); 
#else
	  tempFileName.erase( 0, tempFileName.find_last_of('.') ); 
#endif
	  ending = tempFileName;
	}
    }
  return ending;
}

/**
* @Name of the method: fileNameBase
* @Description:get the name-base (that means, the name without path
* and ending) of the file-name
*
* Precondition: file-name
*
* Postcondition: file-name without path and ending
*/

string fileNameBase( const string &fileName )
{
  string fileNameBase, tempFileName;

  tempFileName = fileName;
  // first check, if there is a path, that means, if there is the
  // letter '/':
  if( tempFileName.find_last_of('/') < tempFileName.size() )
    {
      // remove the path:
#ifdef G_COMPILER
      tempFileName.remove( 0, (tempFileName.find_last_of('/')+1) ); 
#else
      tempFileName.erase( 0, (tempFileName.find_last_of('/')+1) ); 
#endif
      fileNameBase = tempFileName;
    }
  // then check, if there is an ending, that means, if there is the
  // letter '.'
  if( tempFileName.find_last_of('.') < tempFileName.size() )
    {
      // now check, if the point is not part of the path, therefore
      // look first, IF there exists an path and then, If this path is
      // written with a point:
      if( (tempFileName.size() < tempFileName.find_last_of('/')) ||
	  (tempFileName.find_last_of('/') < tempFileName.find_last_of('.')) )
	{
#ifdef G_COMPILER
	  tempFileName.remove( tempFileName.find_last_of('.') );  
#else
	  tempFileName.erase( tempFileName.find_last_of('.') ); 
#endif
	}
    }
  fileNameBase = tempFileName;
  return fileNameBase;
}

/**
* @Name of the method: fileName
* @Description: test, if a file with the given file-name already
* exists, if it exists make new file-name by adding an integer at the
* end (before the ending!!!!)
*
* Precondition: file-name which should be tested
*
* Postcondition: file-name which doesn't exists so far
*/
string fileName( const string &fileName )
{
  int countTestFile = 1;
  bool makeNewOutputName = true;
  string tempFileName, newFileName, ending, testName, tempTestName;

  tempFileName = fileName;
  // first get ending of the file-name:
  ending = fileNameEnding( tempFileName );
  // remove ending (IF existing) of tempFileName:
  if( ending != "" )
    {
      // remove ending:
#ifdef G_COMPILER
      tempFileName.remove( tempFileName.find_last_of('.') ); 
#else
      tempFileName.erase( tempFileName.find_last_of('.') ); 
#endif
    } 
  tempTestName = tempFileName;
  while( makeNewOutputName == true )
    {
      testName = tempTestName + ending;
      ifstream testFile( testName.c_str() );
      if( testFile ) 
	{ tempTestName = tempFileName + "__" + itos( countTestFile++ ); }
      else
	{ makeNewOutputName = false; }
      testFile.close();
    } 
  newFileName = tempTestName + ending;
  return newFileName;
}










