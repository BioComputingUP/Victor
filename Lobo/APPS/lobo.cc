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
 * @Description This program generates a model for a single loop or indel.
 */
#include <string>
#include <GetArg.h>
#include <LoopModel.h>
#include <LoopTableEntry.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <IntCoordConverter.h>
#include <VectorTransformation.h>

#include <LoboTools.h>
    

int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))    {
      showLoboHelp();
      return 1;
    };
  if (getArg( "-help", nArgs, argv))    {
      showMoreLoboHelp();
      return 1;
    };

// -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
// treat command-line arguments first:
// -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-

  string inputFile, outputFile, sequenceFile, scwrlFile, scatterFile, 
    tableFile,chainID;
  vector<char> allCh; 
  unsigned int pdbIndex1, pdbIndex2, num, num2, maxWrite;

  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "target.pdb");
  getArg( "-seq", sequenceFile, nArgs, argv, "null");
  getArg( "-scwrl", scwrlFile, nArgs, argv, "null");
  getArg( "s", pdbIndex1, nArgs, argv, 0);
  getArg( "e", pdbIndex2, nArgs, argv, 0);
  getArg( "c", chainID, nArgs, argv, " ");
  getArg( "-sol1", num, nArgs, argv, LoopModel::MAX_ITER_SOL);
  getArg( "-sol2", num2, nArgs, argv, 1);
  getArg( "-maxWrite", maxWrite, nArgs, argv, 9999);
  getArg( "-scatter", scatterFile, nArgs, argv, "!");
  getArg( "-table", tableFile, nArgs, argv, "!");

  bool verbose = getArg( "-verbose", nArgs, argv);
  bool test = getArg( "-test", nArgs, argv);

  LoopModel lm;

  if (scatterFile != "!")    {
      ofstream scatFile(scatterFile.c_str());
      if (!scatFile)
	ERROR("Could not open file for writing. " + scatterFile, exception);
      lm.setScatterPlot(&scatFile);
    }
  
  if (tableFile != "!")
    lm.setTableFileName(tableFile);

  lm.setVerbose(2);

  double tmpEW = lm.getENDRMS_WEIGHT();
  getOption(tmpEW, "-endrmsWeigth", nArgs, argv);
  lm.setENDRMS_WEIGHT(tmpEW);
  treatSpecialOptions(nArgs, argv);

  getOption(LoopModel::OPT_MAX1, "-opt1", nArgs, argv);
  getOption(LoopModel::OPT_MAX2, "-opt2", nArgs, argv);
  getOption(LoopModel::OPT_NUM, "-optnum", nArgs, argv);

  bool cluster = getOption( "-cluster", nArgs, argv);
  bool norankRms = getOption( "-norankRms", nArgs, argv);
  bool noFullModel = getOption( "-noFullModel", nArgs, argv);
  bool noWrite = getOption( "-noWrite", nArgs, argv);
  bool refine = getOption( "-refine", nArgs, argv);
  bool optimize = getOption( "-optimize", nArgs, argv);
  bool optall = getOption( "-optall", nArgs, argv);
  bool withOxygen = getOption( "-withOxygen", nArgs, argv);
  bool inter = getOption( "-interpol", nArgs, argv);
  if (inter)
    lm.setInterpolated();
  
  if ((inputFile == "!") || (pdbIndex1 == 0) || (pdbIndex2 == 0)) 
    {
      cout << "Missing file specification. Aborting. (-h for help, --help "
	   << "for long help)" << endl;
      return -1;
    }

// -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
// >>>>> END >>>>> treat command-line arguments
// -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-

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
      else 
	chainID=allCh[0];    
    
      pl.setPermissive();
      Protein prot;
      prot.load(pl);
      sp=prot.getSpacer(chainID[0]);
  cout << "Loaded protein.\n" << endl;

  //check if anchor residue indexes are valid:
  if (sp->isGap(pdbIndex1))
    ERROR("Invalid N-terminal anchor residue specified.", exception);

  if (sp->isGap(pdbIndex2))
    ERROR("Invalid C-terminal anchor residue specified.", exception);

  // if sequence file was passed, check if any AAs in the loop are missing:
  if (sequenceFile != "null") 
    {  
      // load seq file
      string seq = sGetSequence(sequenceFile, silent);

      // add missing AAs
      for (unsigned int i = pdbIndex1+1; i < pdbIndex2; i++)
	if (sp->isGap(i))	  { 
	    string tmpType = oneLetter2ThreeLetter(seq[i-1]);

	    if (sp->isGap(i-1))	      {
		ERROR("Internal error adding missing amino acids.", exception);
	      }
	    else
	      {
		sp->insertAminoAfter(tmpType, sp->getIndexFromPdbNumber(i-1));
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

  treatProlineBug(index1, index2, *sp);

  IntCoordConverter icc;

  printSeqAndBfactors(index1, index2,*sp);
  fillLine();

  sp->setStateFromTorsionAngles();
 
  if (verbose)
    printSecStructureAndLoop(index1, index2, *sp);

  if (!test)    {
      vector<string> typeVec;
      for (unsigned int i = index1+1; i < index2+1; i++)
	typeVec.push_back(sp->getAmino(i).getType());
      
      vector<Spacer> vsp;
      vsp = lm.createLoopModel(sp->getAmino(index1), 
	   sp->getAmino(index1+1)[N].getCoords(), sp->getAmino(index2), 
	   sp->getAmino(index2+1)[N].getCoords(), index1, index2, num, 
	   num2, typeVec);

      cout << "ranking...\n";
      
      if (refine)  // refine all solutions
	  lm.refineModel(*sp, index1, index2, vsp);
      
      if (optall)  {    // attempt local minimization of all solutions
	
	  LoopModel::OPT_NUM = 1000;
	  lm.optimizeModel(*sp, index1, index2, vsp);
	}
      
      if (!norankRms)	{
	  lm.rankRawScore(*sp, index1, index2, vsp, maxWrite);
	  lm.doScatterPlot(*sp, index1, index2, vsp);
	}
      
      cout << "-----------------------------------------\n";
      cout << " Orig Loop energy = " << lm.calcOrigEnergy(*sp, index1, index2) 
	   << "  ( ?! ) \n";
      printResults(index1, index2, lm, *sp, vsp, withOxygen, maxWrite);
      
      if (cluster)	{
	  lm.clusterLoops(vsp);
	  printResults(index1, index2, lm, *sp, vsp, withOxygen, maxWrite);
	}
      
      if (optimize)      // attempt local minimization of top solution
	lm.optimizeModel(*sp, index1, index2, vsp);
      
      if (!noWrite)
	saveLoopEnsemble(index1, index2, lm, *sp, vsp, outputFile, maxWrite, 
			 noFullModel);
    }
  else    {
      PdbSaver ps(cout);
      sp->save(ps);
    }

  // if SCWRL output file was requested, write out:
  if (scwrlFile != "null")     {  
      ofstream scwrlOut(scwrlFile.c_str());
      string tmp = lm.getSCWRLConservedSequence(*sp, index1, index2);
      scwrlOut << tmp << "\n";
    }

  fillLine();
  cout << endl;

  return 0;
}
