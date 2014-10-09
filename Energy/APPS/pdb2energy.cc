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
 @Description Calculates the energy of a pdb template
 */
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <SolvationPotential.h>
#include <PolarSolvationPotential.h>
#include <RapdfPotential.h>
#include <TorsionPotential.h>
#include <PhiPsi.h>
#include <PhiPsiOmegaChi1Chi2PreAngle.h>
#include <StatTools.h>
using namespace Victor::Biopool;
using namespace Victor::Energy;
using namespace Victor;

const long double W_RAPDF = 2.5;
const long double W_SOLV = 500.0;
const long double W_HYDB = -50.0;
const long double W_TORS = 350.0;

void sShowHelp(){
  cout << "PDB to Energy\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file for PDB template\n"
       << "\t[-c <id>]  \t ID of chain to load from PDB file\n"
       << "\t[-o <filename>] \t\t Output PDB file with FRST value B-factors\n"
       << "\t[-T] \t\t\t New (Mk2) torsion angle potential \n"
       << "\t[-S] \t\t\t New (Mk2) solvation potential \n"
       << "\t[-r <filename>] \t\t Output file for RAPDF histogram\n"
       << "\t[-s <filename>] \t\t Output file for solvation histogram\n"
       << "\t[-t <filename>] \t\t Output file for torsion histogram\n"
       << "\t[--comp <filename>] \t\t Output file for composite histogram\n"
       << "\t[-C] \t\t\t Cafasp evaluation mode\n"
       << "\t[-w <value>] \t\t\t Sliding window size for histogram (def = 0)\n"
       << "\t[-v] \t\t\t Verbose mode\n"
       << "\t[--start <residue>] \t PDB residue for starting energy computation (def = first)\n"
       << "\t[--end <residue>]  \t PDB residue for ending energy computation (def = last)\n"
       << "\t[--nmr] \t\t Calculate average over NMR ensemble\n"
       << "\n";
}

inline void fillLine(){
  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

void sWriteData(ofstream& out, vector<double> data, unsigned int window, 
Spacer& sp, unsigned int startRes){
  if (data.size() - 1 < window) {
      cerr << "Warning: Invalid window size selected. Reducing "
	   << "to window of 0.\n";
      for (unsigned int i = 0; i < data.size(); i++){
	  out << "  " << setw(4) << sp.getPdbNumberFromIndex(i+startRes) 
	      << "\t" << setw(7) << setprecision(4) << data[i] << "\n";
	}
    }

  for (unsigned int i = 0; i < window; i++) {
      out << "  " << setw(4) << sp.getPdbNumberFromIndex(i+startRes) << "\t" 
	  << setw(7) << setprecision(4); 
      double tmp = 0;
      double count = 0;
      for (unsigned int j = 0; j <= i + window; j++){
	  tmp += data[j];
	  count++;
	}
      out << (tmp / count) << "\n";
    }

  for (int i = window; i < (int)data.size() - 1 - (int)window; i++) {
      out << "  " << setw(4) << sp.getPdbNumberFromIndex(i+startRes) << "\t" 
	  << setw(7) << setprecision(4); 
      double tmp = 0;
      double count = 0;
      for (unsigned int j = i - window; j <= i + window; j++){
	  tmp += data[j];
	  count++;
	}
      out << (tmp / count) << "\n";
    }

  for (int i = data.size() - 1 - window; i < (int)data.size(); i++) {
      out << "  " << setw(4) << sp.getPdbNumberFromIndex(i+startRes) << "\t" 
	  << setw(7) << setprecision(4); 
      double tmp = 0;
      double count = 0;
      for (unsigned int j = i - window; j < data.size(); j++){
	  tmp += data[j];
	  count++;
	}
      out << (tmp / count) << "\n";
    }
}


double sHydrogen(Spacer& sp){
  int count = 0;

  for (int i = 0; i < static_cast<int>(sp.sizeAmino()); i++)
    for (int j = 0; j < static_cast<int>(sp.sizeAmino()); j++)
      if (abs(i-j) > 1)	{
	  double dist = sp.getAmino(i)[N].distance(sp.getAmino(j)[O]);
	  double dist2 = sp.getAmino(i)[N].distance(sp.getAmino(j)[C]);
	  double dist3 = sp.getAmino(i)[CA].distance(sp.getAmino(j)[O]);
	  
	  if ((dist <= 4.0) && (dist >= 2.0) && (dist < dist2)
	      && (dist < dist3)) {
	      count++;
	      break;
	    }
	}
  return count;
}


int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv)) {
      sShowHelp();
      return 1;
    };
  vector<char> allCh; 
  string inputFile, outputFile, solvOut, rapdfOut, torOut, compOut, chainID;
  bool verbose = getArg( "v", nArgs, argv);
  bool cafasp = getArg( "C", nArgs, argv);
  unsigned int window;
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "!");
  getArg( "c", chainID, nArgs, argv, "!");
  getArg( "s", solvOut, nArgs, argv, "!");
  getArg( "r", rapdfOut, nArgs, argv, "!");
  getArg( "t", torOut, nArgs, argv, "!");
  getArg( "-comp", compOut, nArgs, argv, "!");
  getArg( "w", window, nArgs, argv, 0);


  int startRes, endRes;
  getArg( "-start", startRes, nArgs, argv, -99999);
  getArg( "-end", endRes, nArgs, argv, 99999);

  bool nmr = getArg( "-nmr", nArgs, argv);
  bool newTorsion = getArg( "T", nArgs, argv);
  bool newSolvation = getArg( "S", nArgs, argv);
   
  if (inputFile == "!") {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  Potential* solv = NULL;

  if (!newSolvation)
    solv = new SolvationPotential;
  else
    solv = new PolarSolvationPotential;

  RapdfPotential rapdf;
  TorsionPotential* tors = NULL;
  if (!newTorsion)
    tors =  new PhiPsi(10);//default ARCSTEP = 10
  else
    tors =  new PhiPsiOmegaChi1Chi2PreAngle(20);//ARCSTEP 20, ARCSTEP2 40 

  Spacer *sp;
  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);

  
  PdbLoader pl(inFile);
  pl.setNoHAtoms();
  pl.setNoVerbose();
  allCh = pl.getAllChains(); 
  for (unsigned int i = 0; i < allCh.size(); i++)
    cout << "\t," << allCh[i] << ",";
  cout << "\n";
 
  /*check on validity of chain: 
  if user select a chain then check validity
   else select first valid one by default*/
  if (chainID != "!")  {
      bool validChain=false;
      for (unsigned int i = 0; i < allCh.size(); i++ ) {
        if (allCh[i]==chainID[0])  {
          pl.setChain(chainID[0]);
          cout << "Loading chain " << chainID << "\n";
          validChain=true;
          break;
        }
      }
      if (!validChain) {
         cout << "Chain " << chainID << " is not available\n";
        return -1;
      }
      
  }
  else{ chainID[0]=allCh[0];
        cout << "Using chain " << chainID <<"\n";
  }
  Protein prot;
  prot.load(pl);
  
  
  /////////
  

  unsigned int max = pl.getMaxModels();

  if ((!nmr) || (max < 2))   {
      sp=prot.getSpacer(chainID[0]);
      
      // check limits for start & end residue 
      // and transform from PDB to index numbering:
      
      if ((startRes == -99999) || 
	  (startRes < static_cast<int>(sp->getStartOffset())))
	startRes = 1;
      else 	{
	  if (startRes > static_cast<int>(sp->maxPdbNumber()))
	    startRes = sp->sizeAmino() - 1;
	  else
	    startRes = sp->getIndexFromPdbNumber(startRes);
	}
      
      if ((endRes == 99999) || (endRes > (int)sp->maxPdbNumber()))
	endRes = sp->sizeAmino() - 2;
      else	{
	  if (endRes < (int)sp->getStartOffset())
	    endRes = 0;
	  else
	    endRes = sp->getIndexFromPdbNumber(endRes);
	}
      
      if (cafasp)	{
	  long double res = rapdf.calculateEnergy(*sp);
	  long double sres = solv->calculateEnergy(*sp);
	  long double hres = - 1.0 *  sHydrogen(*sp);
	  long double tres = tors->calculateEnergy(*sp);
	  
	  long double sum = W_RAPDF * res + W_SOLV * sres + W_HYDB * hres 
	    + W_TORS * tres;
	  
	  cout.setf(ios::fixed, ios::floatfield);
	  
	  cout << inputFile << "\t" << setw(9) << setprecision(4) 
	       << sum << "\t" << setw(9) << setprecision(4) << res << "\t" 
	       << setw(9) << setprecision(4) << sres << "\t" << setw(9) 
	       << setprecision(4) << hres << "\t" << setw(9) << setprecision(4) 
	       << tres << endl;
	  
	  exit(0);
	}
      
      if (verbose)
	cout << "Calculating energy over residues from " 
	     << startRes << " to " 
	     << endRes << ".\n";
      
      int num = endRes - startRes;
      
      vector<double> rapdfData;
      long double res = 0;
      for (int i = startRes; i <= endRes; i++)	{
	  double tmp = 0;
	  AminoAcid& aa = sp->getAmino(i);
	  string aaType = aa.getType();
	  for (int ii = i+1; ii <= endRes; ii++) {
	      AminoAcid& aa2 = sp->getAmino(ii);
	      string aaType2 = aa2.getType();
	      for (unsigned int j = 0; j < aa.size(); j++)
		for (unsigned int k = 0; k < aa2.size(); k++)
		  tmp += rapdf.calculateEnergy(aa[j], aa2[k], aaType, aaType2);
	    }
	  rapdfData.push_back(tmp);
	  res += tmp;
	}
      
      if (verbose)
	cout << "Rapdf energy =   \t\t";
      if (verbose)
	cout << res << "  \t (avg = " << res/num << " per AA)\n";
      else
	cout << res << "  \t " << res/num << "\n";
      
      vector<double> solvData;
      long double sres = 0;
      for (int i = startRes; i <= endRes; i++)	{
 
	  double tmp = solv->calculateEnergy(sp->getAmino(i), *sp);
	  solvData.push_back(tmp);
	  sres += tmp;
	}
      
      if (verbose)
	cout << "Solvation energy = \t\t";
      if (verbose)
	cout << sres << "  \t (avg = " << sres/num << " per AA)\n";
      else
	cout << sres << "  \t " << sres/num << "\n";
      
      if (verbose)
	cout << "Mainchain hydrogen bonds = \t";
      cout << sHydrogen(*sp) << "\n";
      
      vector<double> torData;
      long double tres = 0;
      for (int i = startRes; i <= endRes; i++){
	  double tmp = tors->calculateEnergy(sp->getAmino(i), *sp);
	  torData.push_back(tmp);
	  tres += tmp;
	}
      
      if (verbose)
	cout << "Torsion energy = \t\t";
      if (verbose)
	cout << tres << "  \t (avg = " << tres/num << " per AA)\n";
      else
	cout << tres << "  \t " << tres/num << "\n";
      
      
      vector<double> compData;
      long double compres = 0;
      for (unsigned int i = 0; i <= rapdfData.size() - 1; i++){
	  double tmp = W_RAPDF * rapdfData[i] + W_SOLV * solvData[i] 
	    + W_TORS * torData[i];
	  compData.push_back(tmp);
	  compres += tmp;
	}
      
      if (verbose)
	cout << "Composite FRST energy = \t";
      if (verbose)
	cout << setw(10) << setprecision(8) << compres << "  \t (avg = " 
	     << setw(10) << setprecision(8) << compres/num << " per AA)\n";
      else
	cout << setw(10) << setprecision(8) << compres << "  \t " 
	     << setw(10) << setprecision(8) << compres/num << "\n";
      
      if (rapdfOut != "!"){
	  cout << "Writing RAPDF histogram to file: " << rapdfOut << "\n";
	  ofstream rOut(rapdfOut.c_str());
	  if (!rOut)
	    ERROR("Could not open file for writing.", exception);
	  
	  sWriteData(rOut, rapdfData, window, *sp, startRes);
	}
      
      if (solvOut != "!"){
	  cout << "Writing solvation histogram to file: " << solvOut << "\n";
	  ofstream sOut(solvOut.c_str());
	  if (!sOut)
	    ERROR("Could not open file for writing.", exception);
	  
	  sWriteData(sOut, solvData, window, *sp, startRes);
	}
      
      if (torOut != "!"){
	  cout << "Writing torsion histogram to file: " << torOut << "\n";
	  ofstream tOut(torOut.c_str());
	  if (!tOut)
	    ERROR("Could not open file for writing.", exception);
	  
	  sWriteData(tOut, torData, window, *sp, startRes);
	}
      
      if (compOut != "!"){
	  cout << "Writing composite histogram to file: " << compOut << "\n";
	  ofstream cOut(compOut.c_str());
	  if (!cOut)
	    ERROR("Could not open file for writing.", exception);
	  
	  sWriteData(cOut, compData, window,*sp, startRes);
	}

      // write output PDB file with FRST scores as B factors
      if (outputFile != "!"){
	  double avgB = average(compData);
	  double sdB = standardDeviation(compData);
	  
	  for (unsigned int i = 0; i < sp->sizeAmino(); i++)
	    for (unsigned int j = 0; j < sp->getAmino(i).size(); j++)
		sp->getAmino(i)[j].setBFac(0.0);

	  for (unsigned int i = 0; i < compData.size(); i++) {
	      double tmp = 50 + 20 * ((compData[i] - avgB) / sdB);
	      for (unsigned int j = 0; j < sp->getAmino(i+startRes).size(); j++)
		sp->getAmino(i+startRes)[j].setBFac(tmp);
	    }

	  ofstream outFile(outputFile.c_str());
	  if (!outFile)
	    ERROR("Could not open file for writing.", exception);
 	  PdbSaver ps(outFile);
	  sp->save(ps);
	}
    }
  else    {
      cout << "NMR ensemble mode.\n";
      cout << "Number of models in structure: " << max << "\n";
      fillLine();
      
      vector<long double> res;
      vector<long double> sres;
      vector<long double> hres;
      vector<long double> tres;
      vector<long double> sum;

      cout << "Energy of model:\n";
      for (unsigned int i = 1; i <= max; i++)	{
	  pl.setModel(i);
	  Spacer *sp;
          prot.load(pl);
	  sp=prot.getSpacer(i);
	  res.push_back(rapdf.calculateEnergy(*sp));
	  sres.push_back(solv->calculateEnergy(*sp));
	  hres.push_back(- 1.0 *  sHydrogen(*sp));
	  tres.push_back(tors->calculateEnergy(*sp));
	  
	  sum.push_back(W_RAPDF * res[i-1] + W_SOLV * sres[i-1] 
			+ W_HYDB * hres[i-1] + W_TORS * tres[i-1]);

	  cout.setf(ios::fixed, ios::floatfield);
	  
	  cout << setw(2) << i
	       << "\t" << setw(9) << setprecision(4) << sum[i-1] 
	       << "\t" << setw(9) << setprecision(4) << res[i-1] 
	       << "\t" << setw(9) << setprecision(4) << sres[i-1] 
	       << "\t" << setw(9) << setprecision(4) << hres[i-1] 
	       << "\t" << setw(9) << setprecision(4) << tres[i-1] << "\n";
	}
      
      fillLine();

      long double minR = minimum(res);
      long double maxR = maximum(res);
      long double avgR = average(res);
      long double sdR = standardDeviation(res);

      long double minS = minimum(sres);
      long double maxS = maximum(sres);
      long double avgS = average(sres);
      long double sdS = standardDeviation(sres);

      long double minH = minimum(hres);
      long double maxH = maximum(hres);
      long double avgH = average(hres);
      long double sdH = standardDeviation(hres);

      long double minT = minimum(tres);
      long double maxT = maximum(tres);
      long double avgT = average(tres);
      long double sdT = standardDeviation(tres);

      long double minF = minimum(sum);
      long double maxF = maximum(sum);
      long double avgF = average(sum);
      long double sdF = standardDeviation(sum);

      cout << "Ensemble data:\n"
	   << "Min"
	   << "\t" << setw(9) << setprecision(4) << minF 
	   << "\t" << setw(9) << setprecision(4) << minR 
	   << "\t" << setw(9) << setprecision(4) << minS 
	   << "\t" << setw(9) << setprecision(4) << minH
	   << "\t" << setw(9) << setprecision(4) << minT << "\n"
	   << "Max"
	   << "\t" << setw(9) << setprecision(4) << maxF 
	   << "\t" << setw(9) << setprecision(4) << maxR 
	   << "\t" << setw(9) << setprecision(4) << maxS 
	   << "\t" << setw(9) << setprecision(4) << maxH
	   << "\t" << setw(9) << setprecision(4) << maxT << "\n"
	   << "Avg"
	   << "\t" << setw(9) << setprecision(4) << avgF 
	   << "\t" << setw(9) << setprecision(4) << avgR 
	   << "\t" << setw(9) << setprecision(4) << avgS 
	   << "\t" << setw(9) << setprecision(4) << avgH
	   << "\t" << setw(9) << setprecision(4) << avgT << "\n"
	   << "SD "
	   << "\t" << setw(9) << setprecision(4) << sdF 
	   << "\t" << setw(9) << setprecision(4) << sdR 
	   << "\t" << setw(9) << setprecision(4) << sdS 
	   << "\t" << setw(9) << setprecision(4) << sdH
	   << "\t" << setw(9) << setprecision(4) << sdT << "\n";

      fillLine();
    }
  delete solv;
  delete tors;
}
