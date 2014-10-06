/** @section LICENCE
*  This file is part of Victor.
*
*    Victor is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Victor is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.

*    You should have received a copy of the GNU General Public License
*    along with Victor.  If not, see <http://www.gnu.org/licenses/>.
*/
/** 
 *@Class:              RankAnalyzer 
  * 
 *@Description:This class implements functions to analyze the ranking output of lobo
  *      
*/
// Includes:
#include <RankAnalyzer.h>
#include <IoTools.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool{
    
static void sLine(){
  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

/**
 *@Description Calculates the correlation between two sets of data
*@param  reference to each of the two sets of data(vector<double>&)
 *@return  corresponding value(static double)
*/
static double sCorrelation(vector<double>& zs1,vector<double>& zs2 ){
    double tmp = 0.0;
    for (unsigned int j = 0; j < zs1.size(); j++)
	tmp += zs1[j] * zs2[j];

    tmp /= zs1.size();
    return tmp;
}


// CONSTRUCTORS/DESTRUCTOR:
/**
 *@Description basic constructor  
*/
RankAnalyzer::RankAnalyzer() {
  PRINT_NAME;
  for (unsigned int i = 0; i < MAX_COL; i++)  {
      vector<double> tmp;
      col.push_back(tmp);
      vector<double> tmp2;
      zScore.push_back(tmp2);
      avg[i] = 0.0;
      sd[i] = 0.0;
  }

  for (unsigned int i = 0; i < MAX_LOOP; i++)  {
      vector<double> tmp;
      result.push_back(tmp);
      vector<double> tmp2;
      resScore.push_back(tmp);
  }

}
/**
* Constructor based on an original object copied
*@param  original object(RankAnalyzer& )
*/
RankAnalyzer::RankAnalyzer(const RankAnalyzer& orig){
  PRINT_NAME;
  this->copy(orig);
}
/**
*  basic destructor
*/
RankAnalyzer::~RankAnalyzer(){
  PRINT_NAME;
  for (unsigned int i = 0; i < MAX_COL; i++)  {
      col[i].clear();
      zScore[i].clear();
  }
  for (unsigned int i = 0; i < MAX_LOOP; i++)  {
      result[i].clear();
      resScore[i].clear();
  }

}

// PREDICATES:
/**
 *@description Method to print the headers
*@param  none
*/
void RankAnalyzer::printHeader(){
    cout << "\t 1     \t 2     \t  3    \t  4    \t  5    \t 6     \t 7   \t"
	 << "  8    \t  9    \t  10 \t  11  \n";
    cout << "\tRMSD   \tScore  \tendRms \tenergy \tsol_en \tsecpre \tpack \t"
	 << "hydrog \tcompct \tflank \tprop \n\n";
}
/**
 *@Description Prints the correlation with a corresponding introduccion
*@param  text to print as introduction (char *), the index where the value to calculate the correlation is(unsigned int)
*/
void RankAnalyzer::printCorrelation(char* intro, unsigned int index){
  cout << intro;
  for (unsigned int i = 0; i < MAX_COL; i++)
      cout << " \t" << setw(5) << setprecision(3) 
	   << sCorrelation(zScore[index],zScore[i]); 
  cout << "\n";
}
/**
 *@Description Prints the correlation with a corresponding introduccion
*@param  text to print as introduction (char *), the vector where the value to calculate the correlation is(vector<double> )
*/
void RankAnalyzer::printCorrelation(char* intro, vector<double> data){
  cout << intro;
  for (unsigned int i = 0; i < MAX_COL; i++)
      cout << " \t" << setw(5) << setprecision(3) 
	   << sCorrelation(data,zScore[i]); 
  cout << "\n";
}

/**
 *@Description prints some results
*@param index of the topvalue(unsigned int), name of file where the top results will be put(string) ,Maximum possible score(double)
 *@return  
*/
void RankAnalyzer::printTopXResults(unsigned int top, string topFile, double maxScore){
    // do file I/O:
    ofstream outFile(topFile.c_str(), ios::app);
    if (!outFile)
	ERROR("Could not write file.", exception);

    unsigned int maxIndex = 0;
    for (unsigned int i = 0; i < MAX_LOOP; i++)
	if (result[i].size() == 0)
	    break;
	else
	    maxIndex++;

    vector<double> bestX;
    vector<double> bestXSc;

    unsigned int notCovered = 0;

    // find Top X answer for every single loop:
    for (unsigned int i = 0; i < maxIndex; i++)    {
	double tmp = result[i][0];
	double tmpSc = resScore[i][0];

	if (tmpSc > maxScore)	{
	    notCovered++;
	    continue;
	}

	unsigned int maxTmp = ( top < result[i].size() ? top 
				: result[i].size() );

	for (unsigned int j = 1; j < maxTmp; j++)	{
	    if (result[i][j] < tmp)
		tmp = result[i][j];

	    if (j >= resScore[i].size())
		ERROR("Index out of scope.", exception);
	    if (resScore[i][j] < tmpSc)
		tmpSc = resScore[i][j];
	}

	bestX.push_back(tmp);
	bestXSc.push_back(tmpSc);
    }

    // now calculate the statistics:

    // average:
    double avgT = 0.0;
    for (unsigned int j = 0; j < bestX.size(); j++)
	avgT += bestX[j];
    avgT = avgT / bestX.size();
    
    double avgSc = 0.0;
    for (unsigned int j = 0; j < bestXSc.size(); j++)
	avgSc += bestXSc[j];
    avgSc = avgSc / bestXSc.size();


    //standard deviation:
    double sdT = 0.0;
    for (unsigned int j = 0; j < bestX.size(); j++)
	sdT += (bestX[j]-avgT) * (bestX[j]-avgT);
    sdT = sqrt(sdT / bestX.size());
    
    double sdSc = 0.0;
    for (unsigned int j = 0; j < bestXSc.size(); j++)
	sdSc += (bestXSc[j]-avgSc) * (bestXSc[j]-avgSc);
    sdSc = sqrt(sdSc / bestXSc.size());
    

    // if topFile is ''valid'', do file output, otherwise screen output:
    if (topFile != "!")
	outFile << setw(3) << maxScore << "\t" << setw(3) 
		<< top << "\t" << setw(5) << setprecision(3)
		<< avgT << "\t" << setw(5) << setprecision(3) 
		<< sdT << "\t" << setw(5) << setprecision(3)
		<< avgSc << "\t" << setw(5) << setprecision(3) 
		<< sdSc << "\t" << setw(5) << setprecision(3) 
		<< 100 * static_cast<double>(maxIndex - notCovered)/maxIndex 
		<< "\n";
    else
	cout << "Top " << setw(3) << top << " results:   avg = " 
	     << setw(5) << setprecision(3)
	     << avgT << " \t sd = " << setw(5) << setprecision(3) 
	     << sdT << "\t score:   avg = "
	     << setw(5) << setprecision(3)
	     << avgSc << "\t sd = " << setw(5) << setprecision(3) 
	     << sdSc << "\t cov% = " << setw(5) << setprecision(3) 
	     << 100 * static_cast<double>(maxIndex - notCovered)/maxIndex 
	     << "\n";
}


// MODIFIERS:
/**
 *@Description loads the file data considering a number of "records"
*@param  reference to the input file(istream&), number of the "records" to consider(unsigned int)
 *@return  
*/
void RankAnalyzer::load(istream& inFile, unsigned int select){
  // result has to be cleared here

  result.clear();
  vector<double> tmpVec;
  result.push_back(tmpVec);
  unsigned int line = 0;

  // read data from file:
  while (inFile)  {
      eatComment(inFile);
      if ((!checkForKeyword(inFile, "LOOP")) && (inFile))
	  ERROR("Invalid data found in file.", exception);

      bool res_set = true;
      unsigned int len, numSol;
      inFile >> len >> numSol;
      skipToNewLine(inFile);

      for (unsigned int c = 0; c < numSol; c++)      {
	  for (unsigned int i = 0; i < MAX_COL; i++)	  {
	      string tmpStr;
	      double tmp = 0.0;
	      inFile >> tmpStr;

	      if (!inFile)  
		  break;   

	      // converts string to extract value:
	      if (tmpStr.length() < 3)
		ERROR("String read from file is too short ?!", exception);
	      cout << "----> " << tmpStr << endl;

	      if (tmpStr[1] == ':')
		tmp = stod(tmpStr.substr(2, tmpStr.length()-2));
	      else if (tmpStr[2] == ':')
		tmp = stod(tmpStr.substr(3, tmpStr.length()-3));
	      else 
		tmp = stod(tmpStr);

	      cout << "-.-.# " << tmp << endl; 
	       
	      if ((select > 0) && (len != select))	      {
		  res_set = false;
		  continue;
	      }
	      else	      {
		  if (i == 0)
		      result[line].push_back(tmp);
		  else if (i == 1)
		      resScore[line].push_back(tmp);

		  col[i].push_back(tmp);
	      }
	  }
	  skipToNewLine(inFile);
      }

      if (res_set)
	  line++;
  }
}

/**
 *@Description calculates the statistics like RMS, energy , etc
*@param  none
 *@return  changes are made internally(void)
*/
void RankAnalyzer::calcStatistics(){
  sLine();
  cout << "\t RMSD \t\t combScore  \t endRms \t energy \t solv_en " 
       << "\t secpref \t packing \t hydrogen \t compact \t flank   \t prop \n";

 // average:
  for (unsigned int i = 0; i < MAX_COL; i++)
      for (unsigned int j = 0; j < col[i].size(); j++)
	  avg[i] += col[i][j];
  
  for (unsigned int i = 0; i < MAX_COL; i++)
      avg[i] = avg[i] / col[i].size();

  cout << " # entries = " << col[0].size() << "\n";
  cout << "Avg: ";
  for (unsigned int i = 0; i < MAX_COL; i++)
      cout << " \t " << setw(6) << avg[i];
  cout << "\n";

  //standard deviation:
  for (unsigned int i = 0; i < MAX_COL; i++)
      for (unsigned int j = 0; j < col[i].size(); j++)
	  sd[i] += (col[i][j]-avg[i]) * (col[i][j]-avg[i]);

  for (unsigned int i = 0; i < MAX_COL; i++)
      sd[i] = sqrt(sd[i] / col[i].size());

  cout << "SD: ";
  for (unsigned int i = 0; i < MAX_COL; i++)
      cout << " \t " << setw(6) << sd[i];
  cout << "\n";
  sLine();

  //Z-scores:
  for (unsigned int i = 0; i < MAX_COL; i++)
      for (unsigned int j = 0; j < col[i].size(); j++)
	  zScore[i].push_back( (col[i][j]- avg[i]) / sd[i] );

}

/**
 *@Description copies the information in an object into another
*@param  reference to the original object
 *@return  changes are made internally(void)
*/
void RankAnalyzer::copy(const RankAnalyzer& orig){
  for (unsigned int i = 0; i < MAX_COL; i++)  {
      col[i].clear();
      for (unsigned int j = 0; j < orig.col[i].size(); j++)
	  col[i].push_back(orig.col[i][j]);

      zScore[i].clear();
      for (unsigned int j = 0; j < orig.col[i].size(); j++)
	  zScore[i].push_back(orig.zScore[i][j]);

      avg[i] = orig.avg[i];
      sd[i] = orig.sd[i];
  }

  for (unsigned int i = 0; i < MAX_LOOP; i++)  {
      result[i].clear();
      resScore[i].clear();
      for (unsigned int j = 0; j < orig.result[i].size(); j++)      {
	  result[i].push_back(orig.result[i][j]);
	  resScore[i].push_back(orig.resScore[i][j]);
      }
  }
}


// HELPERS:
 
}
