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
 @Description Find the contacts from the pdb template
 */
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <Spacer.h>
using namespace Victor::Biopool; 
using namespace Victor;

void sShowHelp(){
  cout << "PDB to Contacts\n"
       << "Tool used for finding the contacts in a pdb file.\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file for PDB template\n"
       << "\t-c <filename> \t chainID input file\n"
       << "\t[-f <filename>] \t CONpred input file\n"
       << "\n";
}


const double cont[] = { 9.5, 11.5,  9.5,  8.4, 10.5,
		        9.0,  9.0, 10.5, 8.75, 10.5,
		       10.0, 8.75,  9.5,  9.0,  9.0,
		        9.0,  9.5, 10.5,10.75, 10.5,  9.5};

int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))    {
      sShowHelp();
      return 1;
    };
  string chainID;
 
  vector<char> allCh;
  string inputFile, conFile; 
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "c", chainID, nArgs, argv, "!");
  getArg( "f", conFile, nArgs, argv, "!");
  if ((inputFile == "!"))     {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

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
  sp=prot.getSpacer(chainID[0]);
  string conpred = "";

  if (conFile != "!")   {
      ifstream input(conFile.c_str());
      if (!input)
	ERROR("File not found.", exception);

      string tmp;
      while (input){
	  tmp = readLine(input);
	  for (unsigned int i = 0; i < tmp.length(); i++)
	    if ((tmp[i] == '-') || (tmp[i] == '+'))
	      conpred += tmp[i];	    
	    else
	      break;
	}
    }

  const double DISTANCE = 8.0;
  double corr = 0.0;
  
   for (unsigned int i = 0; i < sp->sizeAmino(); i++) {
      double count = 0.0;
      for (unsigned int j = 0; j < i; j++)
	if (sp->getAmino(i)[CA].distance(sp->getAmino(j)[CA]) <= DISTANCE)
	  count++;
      for (unsigned int j = i+1; j < sp->sizeAmino(); j++)
	if (sp->getAmino(i)[CA].distance(sp->getAmino(j)[CA]) <= DISTANCE)
	  count++;

      if (count > cont[sp->getAmino(i).getCode()]){
	  cout << "+";
	  if ((i < conpred.length()) && (conpred[i] == '+'))
	    corr++;
	}
      else {
	  cout << "-";
	  if ((i < conpred.length()) && (conpred[i] == '-'))
	    corr++;
	}

      if ((i+1) % 60 == 0)
	cout << "\n";
    }
  cout << "\n";

  if (conpred.length() > 0) {
      cout << "Percentage match with CONpred = " 
	   << setw(6) << setprecision(3) << sp->sizeAmino() * 100 
	   << "%\t (length= " << sp->sizeAmino() << ")\n";
      
      for (unsigned int i = 0; i < sp->sizeAmino(); i++){
	  cout << conpred[i];
	  
	  if ((i+1) % 60 == 0)
	    cout << "\n";
        };
      cout << "\n";
    }
   
}
