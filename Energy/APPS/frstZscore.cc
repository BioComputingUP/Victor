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
 @Description FRST - Function of Rapdf, Solvation and Torsion potentials 
      *Z-Score -  The program is also able to calculate the Z-Score value  
       *as:\n (frst - mean)/standard_deviation.  
       *This program calculates a pseudo-energy to evaluate the quality  
      *of a given protein structural model, as expressed in a single 
      *(real) number. 
       *In Z-score evaluation mean and standard deviation values are calculated on a set of 100 sequences  
       *derived from the original one by inverting and randomly permutating  
      *the original amino acid sequence.
 */
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <SolvationPotential.h>
#include <EffectiveSolvationPotential.h>
#include <RapdfPotential.h>
#include <TorsionPotential.h>
#include <PhiPsi.h>
#include <PhiPsiOmegaChi1Chi2PreAngle.h>
#include <stdlib.h>
#include <time.h>
#include <numeric>  
#include <math.h>


using namespace Biopool;

const long double W_RAPDF = 2.5;
const long double W_SOLV = 500.0;
const long double W_HYDB = -50.0;
const long double W_TORS = 350.0;

//const long double RND_MODELS = 99;

void sShowHelp(){
  cout << "FRST - Function of Rapdf, Solvation and Torsion potentials\n"
       << "Z-Score -  The program is also able to calculate the Z-Score value "
       << "as:\n (frst - mean)/standard_deviation.\n\n" 
       << "This program calculates a pseudo-energy to evaluate the quality "
       << "of a given protein structural model, as expressed in a single "
       << "(real) number.\n"
       << "In Z-score evaluation mean and standard deviation values are calculated on a set of 100 sequences "
       << "derived from the original one by inverting and randomly permutating "
       << "the original amino acid sequence.\n"
       << "Contact the author for further details.\n"
       << "Developed by Silvio Tosatto <silvio.tosatto@unipd.it> \n"
       << "   Options: \n"
       << "\t-I <filelist> \t\t Input *filelist* file for PDB templates\n"
       << "\t-i <filename> \t\t Input file for PDB template\n"
       << "\t-c <id>  \t\t ID of chain to load from PDB file\n"
       << "\t[-z] \t\t\t Get Z-Score Value\n"
       << "\t[--perc <double>] \t Percentace of residues to permute\n"
       << "\t[--rnd <int>] \t\t Number of random models to generate in order to calculate Zscore value\n"
       << "\t[-T] \t\t\t New (Mk2) torsion angle potential \n"
       << "\t[-S] \t\t\t New (Mk2) solvation potential \n"
       << "\t[-v] \t\t\t Verbose mode\n"
       << "\t[-p] \t\t\t Per residue energy calculation\n"
      
       << "\n";
}


double sHydrogen(Spacer& sp){
  int count = 0;
  
  for (unsigned int i = 0; i < sp.sizeAmino(); i++)
    for (unsigned int j = 0; j < sp.sizeAmino(); j++)
      if (abs(i-j) > 1)	{
	  double dist = sp.getAmino(i)[N].distance(sp.getAmino(j)[O]);
	  double dist2 = sp.getAmino(i)[N].distance(sp.getAmino(j)[C]);
	  double dist3 = sp.getAmino(i)[CA].distance(sp.getAmino(j)[O]);
	  
	  if ((dist <= 4.0) && (dist >= 2.0) && (dist < dist2) && (dist < dist3)) {
	      count++;
	      break;
	    }
	}
  return count;
}

//
// Evaluation of the frst value of a given spacer
//

long double frstEval(Spacer& sp, bool perResidue, bool verbose, 
        RapdfPotential &rapdf, Potential* solv, Potential* tors){
  long double res = rapdf.calculateEnergy(sp);
  long double sres = solv->calculateEnergy(sp);
  long double hres = - 1.0 *  sHydrogen(sp);
  long double tres = tors->calculateEnergy(sp);
  
  long double sum = W_RAPDF * res + W_SOLV * sres + W_HYDB * hres 
    + W_TORS * tres;
  
  double divisor = (perResidue == true) ? sp.sizeAmino() : 1;

  if (verbose){
      cout << "rapdf value:" << "\t\t" << setw(9) << setprecision(4) << (res / divisor) << "\n" 
	   << "solvation value:" << "\t" << setw(9) << setprecision(4) << (sres / divisor) << "\n" 
	   << "hydrogen:" << "\t\t" << setw(9) << setprecision(4) << (hres / divisor) << "\n" 
	   << "torsion value:" << "\t\t" << setw(9) << setprecision(4) << (tres / divisor) << "\n\n"
	   << "frst value:" << "\t\t" << setw(9) << setprecision(4) << (sum / divisor) << "\n\n";
    }

  return (sum / divisor);
}


//
// swap aminoacids in positions x and y
//
void swapAmino(Spacer& sp, int x, int y){
  char xa;
  xa = sp.getAmino(x).getType1L(); 
  sp.getAmino(x).setType1L(sp.getAmino(y).getType1L());
  sp.getAmino(y).setType1L(xa);
}

//
// Create a new spacer with amino sequence inverted
//
Spacer reverseAminoSequence(Spacer& sp){
  unsigned int aminoCount = sp.sizeAmino();
  Spacer rev;
  rev.copy(sp);
  
  for (unsigned int i = 0; i < (unsigned int)floor(aminoCount/2); i++)
    swapAmino(rev, i, (aminoCount-1)-i); 

  return rev;
}


//
// generate a new spacer whose aminoacid sequence
// is a random permutation of the original one
// N.B: about 80% of aminoacids are permutated
//
void permuteAminoSequence(Spacer& sp, double permute){
  unsigned int length = sp.sizeAmino();
  unsigned int threshold = (int)(length*permute);
  unsigned int aminoA, aminoB;
  //swap about 80% of aminoacids
  for (unsigned int aa = 0; aa < threshold; aa++)  {
      aminoA = rand()%length;
      aminoB = rand()%length;
      swapAmino(sp, aminoA, aminoB);
    }
}
  
//
//main method
//
int main(int nArgs, char* argv[]){ 
  // Initialize random number generator.
  // generation of a different number each time rand()
  // is called is assured
  srand(time(0));  

  if (getArg( "h", nArgs, argv)){
      sShowHelp();
      return 1;
  };

  string inputFile, inputFilelist, solvOut, rapdfOut, torOut, chainID;
  double pPercentage;
  unsigned int rndModels;
  bool newTorsion = getArg( "T", nArgs, argv);
  bool newSolvation = getArg( "S", nArgs, argv);
  bool verbose = getArg( "v", nArgs, argv);
  bool perResidue = getArg( "p", nArgs, argv);
  bool zScore = getArg( "z", nArgs, argv);
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "I", inputFilelist, nArgs, argv, "!");
  getArg( "c", chainID, nArgs, argv, "!");
  // percentage of permutated aminoacids in order to generate a randon structure
  // default is 80%
  getArg( "-perc", pPercentage, nArgs, argv, 0.8);
  // number of structures to be generated to calculate Zscore, default is 99
  getArg( "-rnd", rndModels, nArgs, argv, 20);
  
  if ((inputFile == "!") && (inputFilelist == "!")){
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  if ((inputFile != "!") && (inputFilelist != "!")){
       cout << "Please choose between filelist and file mode. Aborting. "
	    << "(-h for help)" << endl;
      return -2;     
    }

  Potential* solv = NULL;

  if (!newSolvation)
    solv = new SolvationPotential;
  else
    solv = new EffectiveSolvationPotential;

  RapdfPotential rapdf;

  TorsionPotential* tors = NULL;
  if (!newTorsion)
    tors =  new PhiPsi(10);//default ARCSTEP = 10
  else
    tors =  new PhiPsiOmegaChi1Chi2PreAngle(20);//ARCSTEP 20, ARCSTEP2 40 
   
  ifstream inFile(inputFilelist.c_str());
  if ((!inFile) && (inputFilelist != "!"))
    ERROR("File not found.", exception);

  while ((inFile) || (inputFile != "!")){
      if (inputFilelist != "!")	{
	  inFile >> inputFile;
	  if (!inFile)
	    break;
	}

  Spacer *sp;
  ifstream inFile2(inputFile.c_str());
  if (!inFile2)
	ERROR("File not found.", exception);
  PdbLoader pl(inFile2);
  pl.setNoHAtoms();
  pl.setNoVerbose();
  pl.setPermissive();
  pl.setNoHAtoms();
  vector<char> allCh;
  allCh = pl.getAllChains(); 
  for (unsigned int i = 0; i < allCh.size(); i++)
      cout << "\t," << allCh[i] << ",";
  cout << "\n";
 
  /*check on validity of chain: 
  if user select a chain then check validity
   else select firs valid one by default*/
  if (chainID != "!")  {
      bool validChain=false;
      for (unsigned int i = 0; i < allCh.size(); i++ ) {
        if (allCh[i]==chainID[0]) {
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
    else{
       chainID[0]=allCh[0];
       cout << "Using chain " << chainID <<"\n";
    }
    Protein prot;
    prot.load(pl);
    sp=prot.getSpacer(chainID[0]);
      if (!pl.isValid()){
	  cout << "Warning: Invalid PDB file found.\n";
	  inputFile = "!";
	  continue;
	}
      if (verbose){
	  cout << "\n" << inputFile << "\n\n";
        }
      cout.setf(ios::fixed, ios::floatfield);      
      if (!zScore){
	  if (verbose)
	    frstEval(*sp, perResidue, verbose, rapdf, solv, tors);
	  else{
	      double frstValue = frstEval(*sp, perResidue, verbose, rapdf, solv, tors);
	      cout << setw(9) << setprecision(4) << frstValue << "\n";
	  }
	}
      else{
	vector<long double> energies;
	Spacer reverseModel;
	Spacer randomModel;
	long double revEn;
	long double rndEn;
	//evaluate energy of original spacer
	long double frstValue = frstEval(*sp, perResidue, false, rapdf, solv, tors);
	//remove sidechains of all aminoacids in spacer
	//add CB atom to allow evaluation of torsion potential
	for (unsigned int a = 0; a < sp->sizeAmino(); a++) {	
	    AminoAcid amino = sp->getAmino(a);
	    amino.removeSideChain();
	    amino.patchBetaPosition();
	  }
	//reverse sequence
	reverseModel = reverseAminoSequence(*sp);
	revEn = frstEval(reverseModel, perResidue, false, rapdf, solv, tors);
	energies.push_back(revEn);
	//generate RND_MODELS permutations of aminoacid sequence and calculate their energies
	for (unsigned int m = 1; m <= rndModels; m++) {
	    permuteAminoSequence(reverseModel, pPercentage);
	    rndEn = frstEval(reverseModel, perResidue, false, rapdf, solv, tors);
	    energies.push_back(rndEn);
	  }
	unsigned int length = energies.size();
	//calculate mean value
	long double sum = 0;
	for( unsigned int i = 0; i < length; i++)
	  sum += energies[i];
	long double mean = sum/length;
	//calculate standard deviation
	long double nvariance = 0;
	for( unsigned int i = 0; i < length; i++)
	  nvariance += pow((energies[i]-mean), 2);
	long double stdDev = sqrt(nvariance/length);
	//calculate Z-Score
	long double zScore = (frstValue - mean)/stdDev;
	cout << setw(9) << setprecision(4) << zScore << "\n\n" ;
	if (verbose)  {
	    cout << "mean value:" << "\t\t\t" << setw(9) << setprecision(4) << mean << "\n";
	    cout << "standard deviation:" << "\t\t" << setw(9) << setprecision(4) << stdDev << "\n";
	    cout << "energy frst value of sequence:" << "\t" << setw(9) << setprecision(4) << frstValue << "\n\n";
	  }
      }
      cout << endl;
      // reset variable to trigger break condition:
      inputFile = "!";
   }
  delete tors;
}
