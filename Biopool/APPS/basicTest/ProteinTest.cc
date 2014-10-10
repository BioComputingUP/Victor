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
 
* @Description:     A simple program to test class Protein's features.
*/
#include <GetArg.h>
#include <Protein.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <IoTools.h>
#include <iostream>
using namespace Victor;using namespace Victor::Biopool;

void sShowHelp()
{
  cout << "Protein test\n"
       << "Options: \n"
       << "\t-i <filename> \t Input PDB file\n"
       << "\t-o <filename> \t Output to file\n"
       << "\t-c <id>       \t Chain identifier to read(default is first chain)\n"
       << "\t--all         \t or type '-all' to select all chains analisys)\n"
       << "\t-m <number>   \t Model number to read (NMR only, default is first model)\n"
       << "\t--hetatm      \t type '-hetatm' if you want hetAtm\n"
       << "\t--metalOnly   \t setOnlyMetalHetAtoms; require --hetatm option.\n"
       << "\t--water       \t to enable water selection; require --hetatm option(and metalOnly not).\n"   
       << "\n"
       << "\t Default options: verbose=true; noSecondary=true; noWater= true;\n";
}

void sAddLine()
{
  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

int main(int argc, char* argv[])
{ 
  if (getArg( "h", argc, argv)){
      sShowHelp();
      return 1;
    }
  cout<<"Remember to use -c or --all parameters for chains, and --hetatm to enable Ligands"<<endl;
  string inputFile,outputFile,chainID;
  unsigned int modelNum;
  vector<char> allCh;
  bool hetatm, all, metalOnly, water;
  
  getArg( "i", inputFile, argc, argv, "!");
  getArg( "o", outputFile, argc, argv, "!");
  getArg( "c", chainID, argc, argv, " ");
  all = getArg("-all", argc, argv);
  getArg( "m", modelNum, argc, argv, 999);
  hetatm = getArg("-hetatm", argc, argv);
  metalOnly = getArg("-metalOnly", argc, argv);
  water = getArg("-water", argc, argv);
  if (inputFile == "!"){
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);     
  
  PdbLoader pl(inFile);
  
  pl.setVerbose();
  //pl.setNoSecondary();
  pl.setModel(modelNum);
  //pl.setNoHAtoms();                     //noHAtoms=true

  if (!hetatm)                  
      pl.setNoHetAtoms();               //noHetAtoms=true
  if (all)
      pl.setAllChains();              //allChains=true
  if (metalOnly)
      pl.setOnlyMetalHetAtoms();      //not all cofactors (just metal ions)
  if (water)
      pl.setWater();
  
  
  if (chainID!=" ")  {
      if (chainID.size()>1)
          ERROR("You can choose only 1 chain",error);
      pl.setChain(chainID[0]);
  } 
  cout << "Number of models in structure: " << pl.getMaxModels() << "\n"
       << "Available chains: ";
  allCh = pl.getAllChains(); 
  for (unsigned int i = 0; i < allCh.size(); i++)
    cout << "\t," << allCh[i] << ",";
  cout << "\n";     
  
  Protein prot;
  prot.load(pl);
  cout<<"Protein Loaded\n"
      <<"...\n";
  
 
  if (outputFile == "!"){
      cout << "Missing output file specification. Aborting. (-h for help)" << endl;
      return -1;
  }
  ofstream os(outputFile.c_str());
  if (!os){
    ERROR("Could not open file for writing.", exception);
  }
  
  cout << "Size protein objects: " << prot.sizeProtein() << "\n";
  
    Spacer* sp=NULL;

    for (unsigned int i=0; i<prot.sizeProtein();i++){                 
        sp = prot.getSpacer(i);
        
        cout << i << "\n";
        sp->printGaps();
    }
  
  PdbSaver ps(os);
  ps.setWriteAtomOnly();
  prot.save(ps);
  cout<<"...\n"
      <<"Protein saved"<<endl;
  
  sAddLine();
  cout<<"EXIT LigandTest" << endl;
  return 0;

}
