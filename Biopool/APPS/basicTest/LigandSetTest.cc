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
@Description   A simple program to test class LigandSets features.
                   */
#include <GetArg.h>
#include <Spacer.h>
#include <vector3.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <LigandSet.h>

using namespace Biopool;

void sShowHelp()
{
  cout << "Ligand Test\n"
       << " Options: \n"
       << "\t-i <filename>   \t\t Input PDB file\n"
       << "\t-o <filename>   \t\t Output to file\n"
       << "\t[-d] <distance> \t\t Max distance for ligand 'contact'\n"
       << "\t[-c] <id>       \t\t Chain identifier to read\n"
       << "\t[-m] <number>   \t\t Model number to read (NMR only)\n"
       << "\t[-w]            \t\t To select only het metal ions  \n"   
       << "\n";
}

void sAddLine()
{
  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

int main(int nArgs, char* argv[])
{ 
  if (getArg( "h", nArgs, argv))
    {
      sShowHelp();
      return 1;
    };

  string inputFile,chainID,outputFile;
  double maxDist;
  unsigned int modelNum = 1;
  vector<char> allCh;
  
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "!");
  getArg( "d", maxDist, nArgs, argv, 4.0);
  getArg( "c", chainID, nArgs, argv, "");
  getArg( "m", modelNum, nArgs, argv, 999);
  char* op = "w";
  bool onlyMetalHetAtoms =  getOption( op, nArgs, argv);
  
  if (inputFile == "!")
    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);
   
  PdbLoader pl(inFile);
  pl.setModel(modelNum);
    
  cout << "Number of models in structure: " << pl.getMaxModels() << "\n"
       << "Available chains: ";
  allCh = pl.getAllChains(); 
  for (unsigned int i = 0; i < allCh.size(); i++)
    cout << "\t," << allCh[i] << ",";
  cout << "\n";
  if (chainID != "")
  {
      bool validChain=false;
      for (unsigned int i = 0; i < allCh.size(); i++ )
      {
        if (allCh[i]==chainID[0])
        {
          pl.setChain(chainID[0]);
          cout << "Loading chain " << chainID << "\n";
          validChain=true;
          break;
        }
      }
      if (!validChain)
      {
        cout << "Chain " << chainID << " is not available\n";
        return -1;
      }   
  }
  else 
  {   
      pl.setChain(allCh[0]); 
      cout << "Loading chain " << allCh[0] << "\n";
  }
  if (onlyMetalHetAtoms)
      pl.setOnlyMetalHetAtoms();
  Spacer sp;
  sp.load(pl);
  cout<<"Spacer Loaded"<<endl;
  LigandSet ls;
  ls.load(pl);
  cout<<"LigandSet Loaded"<<endl;

  sAddLine();

  cout << "# of Ligands: " << ls.sizeLigand() <<  "\n";
  cout << "Adjacent residues:\n";

  for (unsigned int i = 0; i < ls.sizeLigand(); i++)
   if (ls.getLigand(i).isCommonMetal())
   {
      cout << " " << i << "  " << ls.getLigand(i).getType() << "\n";

      for (unsigned int j = 0; j < ls.getLigand(i).size(); j++)
      {
	  cout << "\t" << j << "  "<< ls.getLigand(i)[j].getType() << "   ";
	  
	  for (unsigned int k = 0; k < sp.sizeAmino(); k++)
	      for (unsigned int l = 0; l < sp.getAmino(k).size(); l++)
		  if (sp.getAmino(k)[l].distance(ls.getLigand(i)[j]) 
		      <= maxDist)
		      cout << "\t" << k << " (" 
			   << sp.getAmino(k)[l].getType() << ") ";
	  cout << "\n";
      }

   }
//    PdbSaver ps(cout);
//    sp.save(ps);
//    ls.save(ps);
//    ps.endFile();

  sAddLine();
  cout<<"Writing output in Pdb format (list of ligands for all chain)..."<<endl;
  
          
  if (outputFile == "!")
    {
      cout << "Missing output file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  
  ofstream os(outputFile.c_str());
  if (!os)
    ERROR("Could not open file for writing.", exception);
  
  for(unsigned int i=0; i<allCh.size();i++)      //for each chain, save ligands in pdb format
  {
      //PdbLoader pl(inFile2);
      //pl.setModel(modelNum);
      //pl.setNoHAtoms(); //ha effetto solo per loadSpacer
      os<<"Catena "<<allCh[i]<<":"<<endl;
      pl.setChain(allCh[i]);
      
      LigandSet ls;
      ls.load(pl);
      
      PdbSaver ps(os);
      ps.setWriteAtomOnly();
      ls.save(ps);
      os << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
  }
  cout << "EXIT LigandTest" << endl;
  return 0;
}
