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
 * @Description This program clusters the data contained in a Ramachandran distribution file
 */
#include <string>
#include <GetArg.h>
#include <RamachandranData.h>

using namespace Biopool;

int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))
    cout << "This program clusters the data contained in a Ramachandran"
	 << " distribution file.\n"
	 << "Usage: \n" 
	 << " LoopTableTest -A<src1> -B<src2> -O<file> -R<rama> -S<size> " 
	 << "[-v] [-h] \n" 
	 << "\t -i<file> \t Name of input file.\n"
	 << "\t -o<file> \t Name of the output file.\n"
	 << "\t -c<size> \t Cutoff value.\n"
	 << "\t -h \t\t Display this help.\n"
	 << endl;
  
  string inputFile, outputFile;
  double cutoff;

  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "!");
  getArg( "c", cutoff, nArgs, argv, 0.0);
  
  if ((inputFile == "!") || (outputFile == "!") || (cutoff == 0.0))     {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  
  RamachandranData* ramaData = new RamachandranData;
  ifstream ramaIn(inputFile.c_str());
  if (!ramaIn)
    ERROR("Rama filename invalid.", exception);
  ramaData->load(ramaIn);
  
  ramaData->cluster(cutoff);

  ofstream ramaOut(outputFile.c_str());
  if (!ramaOut)
    ERROR("Rama filename invalid.", exception);
  ramaData->save(ramaOut);

  return 0;
}
