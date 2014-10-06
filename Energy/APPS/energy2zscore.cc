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
 *@description Program to convert the  Energy to Z-score
 */
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <IoTools.h>
#include <StatTools.h>
    
using namespace Biopool;

void sShowHelp(){
  cout << "Energy to Z-score\n"
       << "For FRST results.\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file\n"
       << "\t-c <column>   \t\t Select column for ranking (default = 1)\n"
       << "\t[-r <filename>]\t\t Input RMSD file for F.E. & C.C. (optional)\n"
       << "\t[-n <filename>]\t\t Name of native structure for R1 and Z-score\n"
       << "\t[-s]          \t\t Short output format (for parsing)\n"
       << "\t[-v]          \t\t Verbose mode\n"
       << "\t[--rmsd]          \t\t RMSD in 1st column mode\n"
       << "\t[--max]          \t\t Maximize score (e.g. GDT_TS or MaxSub)\n"
       << endl;
}

double sCorrelation(vector<double>& zs1,vector<double>& zs2 ){
    double tmp = 0.0;
    for (unsigned int j = 0; j < zs1.size(); j++)
	tmp += zs1[j] * zs2[j];

    tmp /= zs1.size();
    return tmp;
}


int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv)) {
      sShowHelp();
      return 1;
    };

  string inputFile, rmsdFile, nameNative;
  unsigned int selected;

  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "r", rmsdFile, nArgs, argv, "!");
  getArg( "n", nameNative, nArgs, argv, "!");
  getArg( "c", selected, nArgs, argv, 1);
  bool shortFormat = getArg( "s", nArgs, argv);
  bool verbose = getArg( "v", nArgs, argv);
  if (verbose)
    shortFormat = false;
  bool rmsdPresent = getArg( "-rmsd", nArgs, argv);
  bool maxMode = getArg( "-max", nArgs, argv);

  if (inputFile == "!") {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);

  unsigned int num = 9999;
  unsigned n_prot = 0;
  vector<string> nameData;
  vector<double> energyData, rmsdData;
  double avg = 0.0;
  double sd = 0.0;
  vector<double> zs, energyZS, rmsdZS;
  vector<double> natData;

  while (inFile) {
      string tmp = "xxx";
      long double tmpRmsd = 9999, n1;

      if (!rmsdPresent)
	inFile >> tmp;
      else
	inFile >> tmpRmsd;

      long double tmpN;
      for (unsigned int i = 1; i <= selected; i++)
	inFile >> tmpN;
      n1 = tmpN;
      
      if (!inFile)
	break;
      
      if (tmp == nameNative){
	  num = n_prot;
	  natData.push_back(n1);
	}
      else{
	  nameData.push_back(tmp);
	  energyData.push_back(n1);
	  rmsdData.push_back(tmpRmsd);
	}

      n_prot++;
      skipToNewLine(inFile);
    }

  if ((num == 9999) && (nameNative != "!"))
    ERROR("Native structure not found in decoy set. Check name.", exception);
  
  avg = average(energyData);
  sd = standardDeviation(energyData, avg);

  if (nameNative != "!")
    zs = Zscore(natData, avg, sd);
  else
    zs.push_back(0.0);

  double pearsonCC = 0.0;
  double fracEnr = 0.0;
  unsigned int numFrac = rmsdData.size() / 10;
  unsigned int fractionE = 0;
  double pBest1 = 0.0;
  double pBest10 = 0.0;

  if ((rmsdFile != "!") || (rmsdPresent))  {
      if (rmsdFile != "!") {// RMSD data has to be read from file
	  // read rmsd file:
	  ifstream inRmsdFile(rmsdFile.c_str());
	  if (!inRmsdFile)
	    ERROR("File not found.", exception);
      
	  while (inRmsdFile) {
	      string tmp;
	      double tmpRmsd;
	      inRmsdFile >> tmp >> tmpRmsd;
	      
	      if (!inRmsdFile)
		break;
	      
	      if (verbose){
		  cout << tmp << ";;\t\t;;" << tmpRmsd << ";;\n";
		}

	      for (unsigned int i = 0; i < nameData.size(); i++)
		if (nameData[i] == tmp) {
		    rmsdData[i] = tmpRmsd;
		    break;
		  }
	    }
	  
	  vector<string> nameData2;
	  vector<double> rmsdData2, energyData2;

	  // check for spurious data:
	  for (unsigned int i = 0; i < nameData.size(); i++)
	    if (rmsdData[i] != 9999) {
		nameData2.push_back(nameData[i]);
		energyData2.push_back(energyData[i]);
		rmsdData2.push_back(rmsdData[i]);
	    }

	  rmsdData = rmsdData2;
	  energyData = energyData2;
	  nameData = nameData2;
	}

      if (verbose){
	  cout << "Table data summary:\n";
	  for (unsigned int i = 0; i < rmsdData.size(); i++)
	    cout << nameData[i] << "  \t" << rmsdData[i] << "  \t" 
		 << energyData[i] << "\n";
	}

     // then calculate correlation coefficient:
      energyZS = Zscore(energyData);
      rmsdZS = Zscore(rmsdData);
      pearsonCC = sCorrelation(energyZS, rmsdZS);
      // calculate 10% by 10% fraction enrichment:
      vector<double> sortRmsd;
      for (unsigned int i = 0; i < rmsdData.size(); i++){
	  if (maxMode)
	    rmsdData[i] = -rmsdData[i];
	  sortRmsd.push_back(rmsdData[i]);
	}

      sort(&sortRmsd[0], &sortRmsd[sortRmsd.size()]);

      if (verbose){
	  cout << "Sorted data summary:\n";
	  for (unsigned int i = 0; i < sortRmsd.size(); i++)
	    cout << sortRmsd[i]  << "  \t\t";
	  cout << "\n\n";
	}

      double threshold = sortRmsd[numFrac];
      for (unsigned int i = 0; i <= numFrac; i++)
	if (rmsdData[i] <= threshold)
	  fractionE++;

      fracEnr = (static_cast<double>(fractionE) / numFrac) * 100.0;

      // calculate logP_B1:
      double val1 = rmsdData[0];
      unsigned int r1 = 9999;
      for (unsigned int i = 0; i < sortRmsd.size(); i++)
	if (sortRmsd[i] == val1) {
	    r1 = i;
	    break;
	  }
      
      pBest1 = log10(static_cast<double>(r1+1) / rmsdData.size());

       // calculate logP_B10:
      double val10 = rmsdData[0];
      for (unsigned int i = 1; i < 10; i++)
 	if (rmsdData[i] < val10)
 	  val10 = rmsdData[i];
      
      unsigned int r10 = 9999;
      
      for (unsigned int i = 0; i < sortRmsd.size(); i++)
	if (sortRmsd[i] == val10) {
	    r10 = i;
	    break;
	  }
      pBest10 = log10(static_cast<double>(r10+1) / rmsdData.size());
  }


  if (!shortFormat)    {
      cout << "---------------------------------------------------\n";
      
      if (nameNative != "!")
	cout << nameNative << "\t\t N proteins = " << setw(4) << n_prot << "\n"
	     << "Energy:\t average = " << setw(7) << setprecision(5) 
	     << avg << " (+/- " << setw(7) << setprecision(5) << sd
	     << ") \n\t Z-score = " << setw(3) << setprecision(2)
	     << zs[0] << "\t (rank = " << setw(3) << num+1 << ")" << "\n";
      else
	cout << "\t\t\t N proteins = " << setw(4) << n_prot << "\n"
	     << "Energy:\t average = " << setw(7) << setprecision(5) 
	     << avg << " (+/- " << setw(7) << setprecision(5) << sd
	     << ") \n\t Z-score = " << setw(3) << setprecision(2)
	     << zs[0] << "\n";

      if ((rmsdFile != "!") || (rmsdPresent))
	cout << "Corr.coef. (Pearson) = \t\t" << pearsonCC << "\n"
	     << "Fraction Enrichment (10 by 10) = " 
	     << setw(5) << setprecision(3) << fracEnr << "% ("
	     << fractionE << "/" << numFrac << ")\n"
	     << "log P_best1 = \t\t\t" << setw(4) << pBest1 << "\n"
	     << "log P_best10 = \t\t\t" << setw(4) << pBest10 << "\n";

      cout << "---------------------------------------------------\n";
    }
  else  {
      cout << nameNative << "\t\t" << setw(4) << n_prot << "\t"
	   << setw(3) << setprecision(2) << zs[0] << "\t"
	   << setw(3) << num+1;

      // have to write entries for compatibility:
      cout << "\t" << setw(4) << pearsonCC << "\t" << setw(4) << fracEnr
	   << "\t" << setw(4) << pBest1 << "\t" << setw(4) << pBest10
	   << "\n";      
    }
  return 0;
}
