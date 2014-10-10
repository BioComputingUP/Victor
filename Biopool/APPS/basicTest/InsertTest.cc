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
@Description */

#include <GetArg.h>
#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <vector3.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>

using namespace Victor;using namespace Victor::Biopool;

void sShowHelp()
{
  cout << "Spacer Test\n"
       << "This program tests the spacer generation methods.\n"
       << " Options: \n"
       << "\t-i <filename>      \t\t Input file\n"
       << "\t-o <filename>      \t\t Output to file\n"
       << "\t[-m <number>]      \t\t Model number to read (NMR only)\n"
       << "\t[-s]               \t\t Test sidechain integrity\n"
       << "\t[-w]               \t\t Write output PDB file to COUT\n"
       << "\t[-a]               \t\t Test insertAminoAfter functionality\n"
       << "\t[--sp]             \t\t Test insertSubSpacerAfter functionality\n"
       << "\t[--sub] <filename> \t\t New spacer to insert\n" 
       << endl;
}



int main(int nArgs, char* argv[])
{ 
  if (getArg( "h", nArgs, argv))
    {
      sShowHelp();
      return 1;
    };

  string inputFile, outputFile, chainID, spacerFile;
  unsigned int modelNum = 1;

  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "!");
  getArg( "-sub", spacerFile, nArgs, argv, "!");
  getArg( "m", modelNum, nArgs, argv, 999);
  bool insertSpacer =  getOption( "-sp", nArgs, argv );
  bool insertAmino = getOption( "a", nArgs, argv );
  bool sidechainTest =  getOption( "s", nArgs, argv);
  bool writeOutput =  getArg( "w", nArgs, argv);
  
  
  if (inputFile == "!") 
    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  Spacer sup;
  Spacer sp;
  //Spacer spn;
  ifstream inFile(inputFile.c_str());
 
  if (!inFile)
    ERROR("File not found.", exception);
  PdbLoader pl(inFile);
  pl.setNoHAtoms();
  pl.setModel(modelNum);
  sp.load(pl);
  
  cout << "loaded.\n" << endl;

  fillLine(cout);

  //cout << "Insert now\n";
  
  //unsigned int size = sp.sizeAmino()-1;
  
  if(insertAmino)
    {
      //unsigned int size;
      //sp.makeFlat();
      //sp.insertAminoBefore(aminoAcidThreeLetterTranslator(ALA));
      //for (unsigned int i = 0; i < 5; i++)
      for (AminoAcidCode code = ALA; code <= TYR; code++)
        {
	  //size = sp.sizeAmino()-1;
	  sp.insertAminoAfter(aminoAcidThreeLetterTranslator(code));
	  sp.insertAminoBefore(aminoAcidThreeLetterTranslator(code));
	}

      //spn.insertFirstAmino(sp.getAmino(1).getType(), sp.getAmino(1).getPhi(), sp.getAmino(1).getPsi(),
      //		   sp.getAmino(1).getOmega());
      //spn.insertAminoAfter(sp.getAmino(2).getType());

      cout << "Insertion complete\n";
    }

  if (insertSpacer)
    {
      if (spacerFile == "!") 
	{
	  cout << "Missing file specification. Aborting. (-h for help)" << endl;
	  return -1;
	}
      Spacer nw;
 
      ifstream spFile(spacerFile.c_str());
      
      if (!spFile)
	ERROR("File not found.", exception);
      
      PdbLoader pl(spFile);
      pl.setNoHAtoms();
      pl.setModel(modelNum);
      nw.load(pl);
      
      //sup.insertFirstSpacer(&sp);
      
      //sup.insertSubSpacerAfter(&nw, sup.sizeSpacer()-1);
      
      sup.makeFlat();	
    }
  
  if (sidechainTest)
    {
      for (unsigned int i = 0; i < sp.sizeAmino(); i++)
	{
	  cout << sp.getPdbNumberFromIndex(i) << "\t"
	       << sp.getAmino(i).getType() << "\t";
	  for (unsigned int j = 0; j < sp.getAmino(i).getMaxChi(); j++)
	    cout << sp.getAmino(i).getChi(j) << "\t";
	  cout << "\n";
	}

      fillLine(cout);
    }

  cout << "-------------------------------------------------------\n";
  
  if (insertSpacer)
    {
      if (outputFile != "!")
	{
	  ofstream os(outputFile.c_str());
	  if (!os)
	    ERROR("Could not open file for writing.", exception);
	  PdbSaver ps(os);
	  ps.setWriteAtomOnly();
	  sup.save(ps);
	}
      else if (writeOutput)
	{
	  PdbSaver ps(cout);
	  ps.setWriteAtomOnly();
	  sup.save(ps);
	  
	  fillLine(cout);
	}
    }
  else
    if (outputFile != "!")
      {
	ofstream os(outputFile.c_str());
	if (!os)
	  ERROR("Could not open file for writing.", exception);
	PdbSaver ps(os);
	ps.setWriteAtomOnly();
	sp.save(ps);
      }
    else if (writeOutput)
      {
	PdbSaver ps(cout);
	ps.setWriteAtomOnly();
	sp.save(ps);
	
	fillLine(cout);
      }
  
  cout << "EXIT InsertTest" << endl;
  return 0;
      
}
