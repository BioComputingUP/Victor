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

#include <string>
#include <GetArg.h>
#include <IoTools.h>
#include <StatTools.h>

/**
 @Description Tools used for calculating the correlation coefficient of  two sets od data contained in a unique file.
 */

void sShowHelp(){
  cout << "Correlation\n"
       << "Tool used for calculating the correlation coefficient of"
       << " two sets of data contained in a unique file.\n"
       << "The file must contain only the data organized in column and the"
       << "first column must be a name.\n\n"
       << "The correlation coefficient considered are:\n\n"
       << "1.\t\tThe Pearson Correlation coefficient ( or parametric ).\n"
       << "2.\t\tThe Spearman Correlation coefficient ( or not parametric ).\n\n"
       << "\t-i <filename> \t\t Input file\n"
       << "Options:\n"
       << "\t[-p]          \t\t to have only Pearson Correlation.\n"
       << "\t[-s]          \t\t to have only Spearman Correlation.\n"
       << "\t without options the program return both.\n\n"
       << endl;
}

int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv)) {
      sShowHelp();
      return 1;
    };

  string inputFile;

  getArg( "i", inputFile, nArgs, argv, "!");
  bool pearson = getArg( "p", nArgs, argv);
  bool spearman = getArg( "s", nArgs, argv);

  if (inputFile == "!") {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);
   
  vector<long double> data1, data2;
  int tot = 0;
  string tmp;
  while (inFile) {
      long double num1 = 0.00;
      long double num2 = 0.00;
      
      inFile >>tmp >>num1>>num2;
      cout<<num1<<" "<<num2<<"\n";
      if (!inFile)
	break;
      tot++;

      data1.push_back( num1 );
      data2.push_back( num2 );
    }
  
  long double Pearson = pearsonCorrelation ( data1, data2);
  long double Spearman = spearmanCorrelation ( data1, data2);

  cout << "---------------------------------------------------\n";
  
  cout << "Number of data processed" << setw(4) <<tot<<"\n\n";
  if ( !pearson )
    cout<< "Spearman Correlation : " <<Spearman<<"\n";
  if ( !spearman ) 
    cout<< "Pearson Correlation : " <<Pearson<<"\n";
 
  cout << "---------------------------------------------------\n";

  return 0;
}

