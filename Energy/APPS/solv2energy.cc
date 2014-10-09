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
 @Description calculates the Energy using the solvation distribution
 */
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
using namespace Victor::Biopool; 
using namespace Victor;

const unsigned int MAX_BINS = 15;
const unsigned int BIN_POLAR = 1;

void sShowHelp(){
  cout << "Solvation to Energy\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input solvation distribution file\n"
       << "\t-o <filename> \t\t Output file for solvation histogram\n"
       << "\t[-p] \t\t\t Polarity check mode\n"
       << "\n";
}


int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv)){
      sShowHelp();
      return 1;
    };

  string inputFile, outputFile;
  bool polar = getArg( "p", nArgs, argv);
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "test.out");

  if (inputFile == "!") {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  vector<vector<int> > sum;
  vector<vector<vector<int> > > sumPolar;

  vector<int> tmp;
  for (unsigned int i = 0; i < MAX_BINS; i++)
    tmp.push_back(0);

  for (unsigned int i = 0; i < AminoAcid_CODE_SIZE; i++)
    sum.push_back(tmp);

  for (unsigned int i = 0; i < (MAX_BINS/BIN_POLAR); i++)
    sumPolar.push_back(sum);

  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);

  while (inFile)    {
      unsigned int num;
      string type;
      if (polar){
	  string tmp;
	  inFile >> tmp;
	  if (!inFile)
	    break;

	  if (tmp.substr(0,1) != "P")
	    ERROR("Invalid data found in input file.", exception);
	}

      inFile >> num >> type;

      AminoAcidCode code = aminoAcidThreeLetterTranslator(type);
      if ((num > 1000) || (code == XXX))
	continue;

      for (unsigned int i = 0; i < num; i++)	{
	  unsigned int tmp, tmpPolar;
	  inFile >> tmp;
	  if (polar)  {
	      inFile >> tmpPolar;
	      
	      unsigned int binPolar = tmpPolar / BIN_POLAR;
	      if (binPolar >= (MAX_BINS/BIN_POLAR))
		binPolar = (MAX_BINS/BIN_POLAR)-1;
	      if (tmp >= MAX_BINS)
		tmp = MAX_BINS-1;

	      sumPolar[binPolar][code][tmp]++;
	      sumPolar[binPolar][sum.size()-1][tmp]++;
	    }
	  else  {
	      if (tmp >= MAX_BINS)
		tmp = MAX_BINS-1;
	      
	      sum[code][tmp]++;
	      sum[sum.size()-1][tmp]++;
	    }
	}
      skipToNewLine(inFile);
    }

  ofstream outFile(outputFile.c_str());
  if (!outFile)
    ERROR("File not found.", exception);

  for (unsigned int i = 0; i < sum.size(); i++) {
      unsigned int count = 0;

      if (polar){
	  for (unsigned int k = 0; k < sumPolar.size(); k++)   {
	      count = 0;

	      for (unsigned int j = 0; j < sumPolar[k][i].size(); j++)
		count += sumPolar[k][i][j];
	      
	      outFile << setw(3) << count << "  ";
	    }
	}
      else{
	  for (unsigned int j = 0; j < sum[i].size(); j++)
	    count += sum[i][j];

	  outFile << setw(3) << count << "  ";
	}

      outFile << aminoAcidThreeLetterTranslator(static_cast<AminoAcidCode>(i)) 
	      << "   "; 

      if (polar){
	  for (unsigned int k = 0; k < sumPolar.size(); k++)
	    for (unsigned int j = 0; j < sumPolar[k][i].size(); j++)
	      outFile << "  " << sumPolar[k][i][j];
	}
      else{
	  for (unsigned int j = 0; j < sum[i].size(); j++)
	    outFile << "  " << sum[i][j];
	}

      outFile << "\n";
    }

}
