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
       << "\t-i <filename> \t\t Input file\n"
       << "\t-o <filename> \t\t Output to file\n"
       << "\t-c <id>     \t\t Chain identifier to read\n"
       << "\t-m <number> \t\t Model number to read (NMR only)\n"
      // << "\t[-s]         \t\t Test sidechain integrity\n"
      // << "\t[-t]         \t\t Output main chain torsion angles\n"
       << "\t--deep      \t\t Deep Controll over all chains and angles\n"
       << "\t[-w]         \t\t Write output PDB file to COUT\n"
       << endl;
}

void sDeepControll(PdbLoader pl, ifstream& inFile, string inputFile)
{
  string A ="";
  vector<int> len;
  vector<char> chains = pl.getAllChains();   
  if ( chains.size() == 0 )
    chains.push_back(A[0]);
  
  for ( unsigned int i = 0; i < chains.size(); i++ )    //for each chain I load a spacer
    {
      PdbLoader pl(inFile);
      pl.setNoHAtoms();
      pl.setPermissive();
      Spacer sp;
      if ( chains[i] != A[0] )
	pl.setChain(chains[i]);
      sp.load(pl);
      len.push_back(sp.sizeAmino());   //sizeAmino is: size of components + size of subspacerlist, which are not used at the moment
      for ( unsigned int j = 1; j < sp.sizeAmino(); j++ )
	{
	  sp.getAmino(i).getPhi();
	  sp.getAmino(i).getPsi();
	  sp.getAmino(i).getOmega();
	  int a = sp.getAmino(i).getMaxChi();
	  if ( a == 1)
	    sp.getAmino(i).getChi(0);
	  else if ( a > 1 )
	    {
	      sp.getAmino(i).getChi(0);
	      sp.getAmino(i).getChi(1);
	    }
	}
     
    }
  for ( unsigned int i = 0; i < chains.size(); i++ )
    {
      cout <<inputFile<<" Chain: "<<chains[i]<<" length: "<<len[i];
      cout <<"\tAll angles are OK\n";
    }
}

int main(int nArgs, char* argv[])
{ 
  if (getArg( "h", nArgs, argv))
    {
      sShowHelp();
      return 1;
    };

  string inputFile, outputFile, chainID;
  unsigned int modelNum = 1; 

  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "!");
  getArg( "c", chainID, nArgs, argv, "");
  getArg( "m", modelNum, nArgs, argv, 999);
 // bool sidechainTest =  getOption( "s", nArgs, argv);
 // bool torsionOut =  getOption( "t", nArgs, argv);
  bool writeOutput =  getArg( "w", nArgs, argv);
  bool deep = getArg ( "-deep", nArgs,argv);
  
  if ((inputFile == "!")) 
    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

	//while ( true )
	//{
  Spacer sp;

  ifstream inFile(inputFile.c_str()); 
 
  if (!inFile)
    ERROR("File not found.", exception);
  PdbLoader pl(inFile);
  pl.setNoHAtoms();
  pl.setPermissive();
  pl.setModel(modelNum);
  if ( deep )
    {
      sDeepControll( pl , inFile , inputFile);
      return 0;
    }
  
  pl.setAltAtom('B');
  cout << "Number of models in structure: " << pl.getMaxModels() << "\n"
       << "Available chains: ";
  vector<char> allCh = pl.getAllChains(); 
  for (unsigned int i = 0; i < allCh.size(); i++)
    cout << "\t," << allCh[i] << ",";
  cout << "\n";
 
  /*check on validity of chain: 
  if user select a chain then check validity
   else select firs valid one by default*/
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
  
  sp.load(pl);   
  
  cout << "loaded.\n" << endl;
  Spacer arrMemoryEater[100];
     for ( int i=0; i<100; ++i )
	arrMemoryEater[i]=sp;
     cout << "copied" <<endl;
//
//  fillLine(cout);
//
//  if (sidechainTest)
//    {
//      for (unsigned int i = 0; i < sp.sizeAmino(); i++)
//	{
//	  cout << sp.getPdbNumberFromIndex(i) << "\t"
//	       << sp.getAmino(i).getType() << "\t";
//	  for (unsigned int j = 0; j < sp.getAmino(i).getMaxChi(); j++)
//	    cout << sp.getAmino(i).getChi(j) << "\t";
//	  cout << "\n";
//	}
//
//      fillLine(cout);
//    }
//
//  cout << "-------------------------------------------------------\n";
// 
//  if (torsionOut)
//    {
//      for (unsigned int i = 0; i < sp.sizeAmino(); i++)
//       {
//	 cout << setw(3) << sp.getPdbNumberFromIndex(i) 
//	      << "  " << sp.getAmino(i).getType() 
//	      << "   " << setw(7) << setprecision(4) << sp.getAmino(i).getPhi()
//	      << "   " << setw(7) << setprecision(4) << sp.getAmino(i).getPsi()
//	      << "\n";
//       }
//
//      cout << "-------------------------------------------------------\n";
//    }
//
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
  else 
      cout<<"Output not specified\n";
  cout << "EXIT SpacerTest" << endl;
  //}
  return 0;
      
}
