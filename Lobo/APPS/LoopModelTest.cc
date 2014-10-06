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
 * @Description This program generates a model for a single loop.
 */

#include <string>
#include <GetArg.h>
#include <LoopModel.h>
#include <limits.h>
#include <LoopTableEntry.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <XyzSaver.h>
#include <IntCoordConverter.h>
#include <VectorTransformation.h>
    
using namespace Biopool;

void sShowHelp(){
  LoopModel lm;
  cout << "Loop Model Test\n"
       << "This program generates a model for a single loop.\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file\n"
       << "\t-c <chain> \t\t pdb chain\n"   
       << "\t-s <n-start> \t\t Aminoacid where loop starts\n" 
       << "\t-e <n-end> \t\t Aminoacid where loop ends\n"
       << "\t[-o <filename>] \t Output file (default = test.pdb)\n"
       << "\t[--sol1 <n-solutions>] \t Number of primary solutions"
       << " (default= " << LoopModel::MAX_ITER_SOL <<")\n"
       << "\t[--sol2 <n-solutions>] \t Number of secondary solutions"
       << " (default= 1)\n"
       << "\t[--cluster] \t\t Cluster results.\n" 
       << "\t[--norankRms] \t\t Do not rank or filter results according "
       << "to end AA RMS.\n"
       << "\t[--noFullModel] \t Do not write full model, write only"
       << " loop instead.\n"
       << "\t[--noWrite]\t\t Do not write output file.\n"
       << "\t[--endrmsWeigth] \t Choose weigthing of endrms"
       << " (default= " << lm.getENDRMS_WEIGHT() << ")\n"
       << "\t[--energyWeigth] \t Choose weigthing of energy"
       << " (default= " << LoopModel::ENERGY_WEIGTH << ")\n"
       << "\t[--secprefWeigth] \t Choose weigthing of secondary preference"
       << " (default= " << LoopModel::SECPREF_WEIGTH << ")\n"
       << "\t[--secprefTol] \t\t Choose tolerance for secondary preference"
       << " (default= " << LoopModel::SECPREF_TOL << ")\n"
       << "\t[--packingWeigth] \t Choose weigthing of packing density"
       << " (default= " << LoopModel::PACKING_WEIGTH << ")\n"
       << "\t[--weigthEP] \t\t Choose weigthing of EP compared to ED & EN"
       << " (default= " << LoopTableEntry::LAMBDA_EP << ")\n"
       << "\t[--maxSearch] \t\t Choose max fraction of tables to search "
       << "for best result (default= " << LoopTable::MAX_FACTOR << ")\n"
       << "\t[--vdwLimit] \t\t Choose threshold for VDW filter"
       << " (default= " << LoopModel::VDW_LIMIT << ")\n"
       << "\t[--energyLimit] \t Choose threshold for energy filter"
       << " (default= " << LoopModel::ENERGY_LIMIT << ")\n"
       << "\t[--simLimit] \t\t Choose threshold for clustering similarity"
       << " (default= " << LoopModel::SIM_LIMIT << ")\n"
       << "\t[--seqFile] \t\t Sequence file (in one letter code)\n"
       << "\t[--useZScore]\t\t Use Z score ranking scheme instead of "
       << "raw scores.\n"
       << "\t[--refine] \t\t Refines the solutions (ie. reduced endRms)\n"
       << "\t[--optimize] \t\t Local optimization of the top solution\n"
       << "\t[--optall] \t\t Local optimization of all solutions pre-ranking\n"
       << "\t[--opt1] \t\t Number of outer iterations for optimization\n"
       << "\t[--opt2] \t\t Number of inner iterations for optimization\n"
       << "\t[--optnum] \t\t Number of solutions to optimize\n"
       << "\t[--withOxygen] \t\t Include Oxygen atoms in RMSD calculation\n"
       << "\t[--interpol] \t\t Use interpolated RAPDF energy\n"
       << "\t[--scatter <filename>] \t Write a scatter plot file.\n";
}


void sGetOption(double& param, char* optName, int nArgs, char* argv[], 
bool verbose = true){
  double tmp = -100000.0;
  getArg(optName, tmp, nArgs, argv, -100000.0);
  if (tmp > -100000.0)    {
      param = tmp;
      if (verbose)
	cout << optName << " = " << tmp << "\n";
    }
}

void sGetOption(unsigned int& param, char* optName, int nArgs, char* argv[], 
bool verbose = true){
  unsigned int tmp = INT_MAX;
  getArg(optName, tmp, nArgs, argv, INT_MAX);
  if (tmp != INT_MAX)    {
      param = tmp;
      if (verbose)
	cout << optName << " = " << tmp << "\n";
    }
}

void sGetOption(float& param, char* optName, int nArgs, char* argv[], 
bool verbose = true){
  float tmp = -100000.0;
  getArg(optName, tmp, nArgs, argv, -100000.0);
  if (tmp > -100000.0)    {
      param = tmp;
      if (verbose)
	cout << optName << " = " << tmp << "\n";
    }
}


string sGetSequence(string sequenceFile){
  cout << "seqFile = " << sequenceFile << "\n";

  ifstream seqFile(sequenceFile.c_str());
  if (!seqFile)
    ERROR("Sequence file not found.", exception);

  string tmp;
  seqFile >> tmp;
  return tmp;
}



int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))    {
      sShowHelp();
      return 1;
    };

  string inputFile, outputFile, sequenceFile, scatterFile, chainID;
  vector<char> allCh; 
  unsigned int pdbIndex1, pdbIndex2, num, num2;

  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "c", chainID, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "test.pdb");
  getArg( "-seqFile", sequenceFile, nArgs, argv, "null");
  getArg( "s", pdbIndex1, nArgs, argv, 0);
  getArg( "e", pdbIndex2, nArgs, argv, 0);
 
  getArg( "-sol1", num, nArgs, argv, LoopModel::MAX_ITER_SOL);
  getArg( "-sol2", num2, nArgs, argv, 1);

  getArg( "-scatter", scatterFile, nArgs, argv, "!");

  ofstream scatFile(scatterFile.c_str());
  if (!scatFile)
      ERROR("File not found.", exception);
  
  LoopModel lm;  
  lm.setVerbose();

   
  double tmpEW = lm.getENDRMS_WEIGHT();
  sGetOption(tmpEW, "-endrmsWeigth", nArgs, argv);
  lm.setENDRMS_WEIGHT(tmpEW);
  sGetOption(LoopModel::ENERGY_WEIGTH, "-energyWeigth", nArgs, argv);
  sGetOption(LoopModel::SECPREF_WEIGTH, "-secprefWeigth", nArgs, argv);
  sGetOption(LoopModel::SECPREF_TOL, "-secprefTol", nArgs, argv);
  sGetOption(LoopModel::PACKING_WEIGTH, "-packingWeigth", nArgs, argv);
  sGetOption(LoopTableEntry::LAMBDA_EP, "-weigthEP", nArgs, argv);
  sGetOption(LoopTableEntry::LAMBDA_ED, "-weigthED", nArgs, argv);
  sGetOption(LoopTableEntry::LAMBDA_EN, "-weigthEN", nArgs, argv);
  sGetOption(LoopTable::MAX_FACTOR, "-maxSearch", nArgs, argv);
  sGetOption(LoopModel::VDW_LIMIT, "-vdwLimit", nArgs, argv);
  sGetOption(LoopModel::ENERGY_LIMIT, "-energyLimit", nArgs, argv);
  sGetOption(LoopModel::SIM_LIMIT, "-simLimit", nArgs, argv);
  sGetOption(LoopModel::OPT_MAX1, "-opt1", nArgs, argv);
  sGetOption(LoopModel::OPT_MAX2, "-opt2", nArgs, argv);
  sGetOption(LoopModel::OPT_NUM, "-optnum", nArgs, argv);

  bool cluster = getArg( "-cluster", nArgs, argv);
  if (cluster)
    cout << "cluster\n";

  bool norankRms = getArg( "-norankRms", nArgs, argv);
  if (norankRms)
    cout << "norankRms\n";

  bool noFullModel = getArg( "-noFullModel", nArgs, argv);
  if (noFullModel)
    cout << "noFullModel\n";

  bool noWrite = getArg( "-noWrite", nArgs, argv);
  if (noWrite)
    cout << "noWrite\n";

  bool useZScore = getArg( "-useZScore", nArgs, argv);
  if (useZScore)
    cout << "useZScore\n";

  bool refine = getArg( "-refine", nArgs, argv);
  if (refine)
    cout << "refine\n";

  bool optimize = getArg( "-optimize", nArgs, argv);
  if (optimize)
    cout << "optimize\n";

  bool optall = getArg( "-optall", nArgs, argv);
  if (optall)
    cout << "optall\n";

  bool withOxygen = getArg( "-withOxygen", nArgs, argv);
  if (withOxygen)
    cout << "withOxygen\n";

  bool inter = getArg( "-interpol", nArgs, argv);
  if (inter)    {
      cout << "interpol\n";
      lm.setInterpolated();
    }

  if ((inputFile == "!") || (pdbIndex1 == 0) || (pdbIndex2 == 0)) 
    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  cout << "Start, input: " << inputFile << " \n";
  Spacer *sp;

  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);
  PdbLoader pl(inFile);
  pl.setNoHAtoms();
  allCh = pl.getAllChains(); 
  for (unsigned int i = 0; i < allCh.size(); i++)
           cout << "\t," << allCh[i] << ",";
  cout << "\n";

    /*check on validity of chain: 
    if user select a chain then check validity
     else select first valid one by default*/
  if (chainID != " ")  {
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
    
  pl.setPermissive();
  Protein prot;
  prot.load(pl);
  sp=prot.getSpacer(chainID[0]);

  cout << "loaded.\n" << endl;

  //check if anchor residue indexes are valid:
  if (sp->isGap(pdbIndex1))
    ERROR("Invalid N-terminal anchor residue specified.", exception);

  if (sp->isGap(pdbIndex2))
    ERROR("Invalid C-terminal anchor residue specified.", exception);

  // if sequence file was passed, check if any AAs in the loop are missing:
  if (sequenceFile != "null")     {  
      // load seq file
      string seq = sGetSequence(sequenceFile);

      // add missing AAs
      for (unsigned int i = pdbIndex1+1; i < pdbIndex2; i++)
	if (sp->isGap(i))	  { 
	    string tmpType = oneLetter2ThreeLetter(seq[i-1]);

	    if (sp->isGap(i-1))	      {
		ERROR("Internal error adding missing amino acids.", exception);
	      }
	    else	      {
		//sp.insertAminoAfter() (sp.getIndexFromPdbNumber(i-1), tmpType);
		sp->removeGap(i);
	      }
	  }
    }
  else {// check that all loop AAs are valid: 
   
      for (unsigned int i = pdbIndex1+1; i < pdbIndex2; i++)
	if (sp->isGap(i))
	  ERROR("Invalid loop residue specified without sequence template.", 
		exception);
    }

  // adapt indexes to internal AA count:
  unsigned int index1 = sp->getIndexFromPdbNumber(pdbIndex1); 
  unsigned int index2 = sp->getIndexFromPdbNumber(pdbIndex2);
  if ((index1 >= sp->sizeAmino()) || (index2 >= sp->sizeAmino()))
    ERROR("Invalid index encountered.", exception); 

  if (sp->getAmino(index1)[C].distance(sp->getAmino(index1+1)[N]) == 0)    {
      cout << "--> NB: Patching N-terminal anchor region.\n";
      index1--;
    };

  IntCoordConverter icc;

  cout << "---------------------------------------------------\n";
  cout << "Sequence & B-factors:\n";
 
  for (unsigned int i = index1+1; i < index2+1; i++)    {
    cout << "   " << sp->getAmino(i).getType();
    for (unsigned int j = 0; j < sp->getAmino(i).sizeBackbone(); j++)
      cout << "   " << setw(5) << setprecision(3) 
	   << sp->getAmino(i)[j].getBFac();
    cout << "\n";
    }

  cout << "---------------------------------------------------\n";
 
  sp->setStateFromTorsionAngles();
 
  for (unsigned int i = 0; i < sp->sizeAmino(); i++)    {
      if (sp->getAmino(i).getState() == HELIX)
	cout << "H";
      else if (sp->getAmino(i).getState() == STRAND)
	cout << "E";
      else
	cout << "-";
      if ((i+1) % 60 == 0)
	cout << "\n";
    }
  cout << "\n";
  for (unsigned int i = 0; i < sp->sizeAmino(); i++)    {
    if ((i >= index1) && (i <= index2))
      cout << "#";
    else
      cout << ".";
      if ((i+1) % 60 == 0)
	cout << "\n";
    }
  cout << "\n";

  cout << "---------------------------------------------------\n";

  vector<string> typeVec;

  for (unsigned int i = index1+1; i < index2+1; i++)
    typeVec.push_back(sp->getAmino(i).getType());

  vector<Spacer> vsp;
  cout<<sp->getAmino(index1).getType() <<" "<<sp->getAmino(index1+1).getType()<<" ";
  cout<< sp->getAmino(index2).getType()<<" "<<sp->getAmino(index2+1).getType()<<" "<< index1<<" "<< index2<<" "<< num<<" "<<num2<<" "<<typeVec[0];
  vsp =
          lm.createLoopModel(sp->getAmino(index1), 
	   sp->getAmino(index1+1)[N].getCoords(), sp->getAmino(index2), 
	   sp->getAmino(index2+1)[N].getCoords(), index1, index2, num, 
	   num2, typeVec);

  cout << "ranking...\n";
  if (refine)
    {
      // refine all solutions
      lm.refineModel(*sp, index1, index2, vsp);
    }

  if (optall)     { // attempt local minimization of all solutions
      LoopModel::OPT_NUM = 1000;
      lm.optimizeModel(*sp, index1, index2, vsp);
    }

  if (!norankRms)    {
      if (useZScore)

	lm.rankRawScore(*sp, index1, index2, vsp);

      lm.doScatterPlot(*sp, index1, index2, vsp);
    }

  cout << "-----------------------------------------\n";
  cout << " Orig Loop energy = " << lm.calcOrigEnergy(*sp, index1, index2) 
       << "  ( ?! ) \n";
  cout << "-----------------------------------------\n";


  cout << " Results:  \t\t\t\t\t\t    1.35     121     180\n";

  for (unsigned int i = 0; i < vsp.size(); i++)    {
      cout << setw(3) << i << "   ";
      lm.calculateRms(*sp, index1, index2, vsp[i], true, withOxygen);
    }
  cout << endl;

  if (cluster)    {
      lm.clusterLoops(vsp);
      
      cout << "-----------------------------------------\n";
      
      cout << " Results: (clustered)\t\t\t\t\t    1.35     121     180\n";
      
      for (unsigned int i = 0; i < vsp.size(); i++)	{
	  cout << setw(3) << i << "   ";
	  lm.calculateRms(*sp, index1, index2, vsp[i], true, withOxygen);
	}
      cout << endl;
    }

  if (optimize)      // attempt local minimization of top solution
      lm.optimizeModel(*sp, index1, index2, vsp);

  if (!noWrite)
    for (unsigned int i = 0; i < vsp.size(); i++)      {
	string tmpStr = outputFile + "." + itos(i);
	
	ofstream outFile(tmpStr.c_str());
	if (!outFile)
	  ERROR("File not found.", exception);
	PdbSaver ps(outFile);
	
	if (!noFullModel)	  {
	    lm.setStructure(*sp, vsp[i], index1, index2);
	    sp->save(ps);
	  }
	else	  {
	    vsp[i].save(ps);
	  }
      };

  cout << "-----------------------------------------\n";

 return 0;
}


