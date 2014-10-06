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
 * @Description This program runs through filelist and generates a commented 
 * list of PDB IDs and loop start / end residues.
 */
#include <string>
#include <GetArg.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <LoopExtractor.h>
#include <String2Number.h>
#include <globalStatistic.h>

using namespace Biopool;

void sShowHelp(){
  cout << "createLoopTestset\n"
       << "This program runs through filelist and generates a commented " 
       << "list of PDB IDs and loop start / end residues.\n"
       << " Options: \n"
       << "\t-i <Input file> \t\t List of PDB-files (default = " 
       << "filelist)\n"
       << "\t-o <Outfput file> \t\t Output file with loop definitions\n"
       << "\t[--verbose]       \t\t Verbose mode\n";
  exit(1);
}


int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))
    sShowHelp();
  
  unsigned int index1 = 0; // contains the start index of the loop
  unsigned int index2 = 0; // contains the end index of the loop

  int start = 0;           // contains the beginning of another loop in 
                           // the pdb-file
  int end = 0;             // contains the end of another loop in the pdb-file
  string filename;         // contains the name of the current pdb file
  globalStatistic gs;      // contains the overall statistic for the whole 
                           // thing (e.g. propensities) 

   string outputFile, inputFile, endrmsFile,chainID;
  vector<char> allCh; 
   getArg("i", inputFile, nArgs, argv, "filelist");
   getArg("o", outputFile, nArgs, argv, "!");
   
   cout << "input file: " << inputFile << "\n";

   if ( outputFile == "!")     {
      cout << "Missing file specification. Aborting. (-h for help)" << "\n";
      return -1;
    }

  bool verbose = getArg( "-verbose", nArgs, argv);

  cout << "Start\n";
  cout << "---------------------------------------------\n";
   
  // open the file in which all the names of the relevant pdb-files are stored
  ifstream pdb_Filelist(inputFile.c_str());
  if (!pdb_Filelist)
    ERROR("Input file not found.", exception);

  ofstream out(outputFile.c_str());
  if (!out)
    ERROR("Could not create file.", exception);

  // process all files in filelist
  while(pdb_Filelist)    {
      // open the current file
      string pdbname = "";
      pdb_Filelist >> pdbname;
      filename = pdbname + ".pdb";
      ifstream inFile(filename.c_str());
      cout << "\nOpening file named: " << filename.c_str() << "\n";
      if (!inFile)	{
	  cout << "file from filelist not found! Name: " << filename << "\n";
	  continue;
	}
      
      // open the file for the statistic output for each pdb file
      string::size_type pos = filename.rfind("pdb"); 
      if (pos == string::npos) 
	{
	  cout << "Not a *.pdb file, skipping to next file" << "\n";
	  continue;
	}
      
      Spacer *sp;
      PdbLoader pl(inFile);
      pl.setNoHAtoms();
       allCh = pl.getAllChains(); 
      
      cout << "Using chain\n"<<chainID;
      
    /*check on validity of chain: 
    if user select a chain then check validity
     else select first valid one by default*/
      chainID=allCh[0];
    
    
      pl.setPermissive();
      Protein prot;
      prot.load(pl);
      sp=prot.getSpacer(chainID[0]);
      cout<<sp->sizeAmino();
      // write secondary structure signature:
      for (unsigned int i = 0; i < sp->sizeAmino(); i++)	{
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
      
      
      LoopExtractor le;
      le.setSpacer(sp);
      
      while (start != -1) 	{
	  le.nextLoop (start, end);
	  if (verbose)
	    cout << "Start: " << start << " end: " << end << "\n";

	  if (end - start > gs.LOOP_MAX - 1 || end - start < gs.LOOP_MIN - 1 ) 	    {
	      if (verbose)
		cout << "Skipping loop due to length, l= " << end-start << "\n";
	      continue;
	    }
	  if (start <= 1 || (unsigned int)end == sp->sizeAmino() - 1)	    {
	      if (verbose)
		cout << "Skipping loop due to position.\n";
	      continue;
	    }
	  
	  // check for reasonable B-factors in loop:
	  double maxB = 0.0;
	  
	  for ( int i = start-2; i < end+1; i++)
	    for (unsigned int j = 0; j < sp->getAmino(i).sizeBackbone(); j++)
	      if (sp->getAmino(i)[j].getBFac() > maxB)
		maxB = sp->getAmino(i)[j].getBFac();
	  
	  if (maxB > 25.0)	    {
	      if (verbose)
		cout << "Bad B-factor (" << maxB<< ") found.\n";
	      continue;
	    }
	  
	  // now we have to set index2 and index1
	  index2 = end;                     // here we had end+1
	  index1 = start - 1;               // here we had start
	  cout << "index1 (-s): " << index1 << " index2 (-e) " << index2 << "\n";
	  double looplaenge = index2 - index1 - 2;              
	  // this stuff is soleley used for the propensities
	  if (looplaenge <= 0)                                  
	    looplaenge = 1;                                     
	  // just that we don't divide by something stupid later on
	  
	  
	  // write loop signature
	  if (verbose)	    {
	      for (unsigned int i = 0; i < sp->sizeAmino(); i++)	{
		  if ((i >= index1) && (i <= index2))
		    cout << "#";
		  else
		    cout << ".";
		  if ((i+1) % 60 == 0)
		    cout << "\n";
		}
	      cout << "\n";

	      cout << "Sequence & B-factors:\n";
	      for (unsigned int i = index1+1; i < index2+1; i++){
		  cout << "   " << sp->getAmino(i).getType();
		  for (unsigned int j = 0; j < sp->getAmino(i).sizeBackbone(); 
		       j++)
		    cout << "   " << setw(5) << setprecision(3) 
			 << sp->getAmino(i)[j].getBFac();
		  cout << "\n";
		}
	    }

	  // write output data to file:

	  out << pdbname << "\t\t" << sp->getPdbNumberFromIndex(index1) 
	      << "\t" << sp->getPdbNumberFromIndex(index2) << "\n";
	  
	}
      start = 0;

    }
  out.close();
    
  exit(0);
}
