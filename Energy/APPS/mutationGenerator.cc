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
 @DescriptionThis program give a file containing a protein sequence with a single mutations. 
 *    The mutation must not change the structure of the protein.
 *      This file must be used with scwr for have the structure of the new protein
 */
 
#include <string>
#include <GetArg.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <PdbSaver.h>

using namespace Biopool;

void sShowHelp(){
  cout << "MutationGenerator.\n"
       << "This program give a file containing a protein sequence with a single mutations.\n "
       << "The mutation must not change the structure of the protein.\n"
       << "This file must be used with scwr for have the structure of the new protein "
       << "Contact the author for further details.\n"
       << "Developed by Silvio Tosatto <silvio.tosatto@unipd.it> and Alessandro Albiero <alex@cribi.unipd.it> \n"
       << "\n\n"
       << "\t -i <filename>\t\tThe pdb file of the native structure.\n"
       << "\t --id \t\tThe pdb Id position of the mutation.\n"
       << "\t -c  \t\t ID of chain to load from PDB file\n"
       << "\t -a \t\tThe new aa to insert (three letter code in capital letter).\n"
       << "\t -o \t\t The output file."
       << endl; 
}

string sSCWRLSequence(const Spacer& sp, unsigned int index1, unsigned int index2, char mut){
  string tmp = "";
  
  for (unsigned int i = 0; i < index1; i++)
    tmp += static_cast<char>(tolower(threeLetter2OneLetter(sp.getAmino(i).getType())));
  
  for (unsigned int i = index1; i < index2; i++)
    tmp += toupper(mut);
  
  for (unsigned int i = index2; i < sp.sizeAmino(); i++)
    tmp += static_cast<char>(tolower(threeLetter2OneLetter(sp.getAmino(i).getType())));
  
  return tmp;
}



int main(int nArgs, char* argv[]){ 
    string chainID;
  vector<char> allCh; 
  if (getArg( "h", nArgs, argv)) {
      sShowHelp();
      return 1;
    };
  if (getArg( "-help", nArgs, argv)) {
      sShowHelp();
      return 1;
    };

  string inputFile, scwrlFile, aa;
  int ID;
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "a", aa , nArgs, argv, "!");
  getArg( "o", scwrlFile, nArgs, argv, "!");
  getArg( "-id", ID, nArgs, argv, 0);
  getArg( "c", chainID, nArgs, argv, "!");
  if (inputFile == "!") {
      cout << "Missing input file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  if ( ID == 0) {
      cout << "Missing the position of the mutation. Aborting. (-h for help)" << endl;
      return -1;
    }
  if ( aa == "!" )  {
      cout << "Missing mutation. Aborting. (-h for help)" << endl;
      return -1;
    }
  if (scwrlFile == "!") {
      cout << "Missing output file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);
  
  AminoAcidCode code = aminoAcidThreeLetterTranslator(aa);
  char mut = aminoAcidOneLetterTranslator(code);
  Spacer *sp;

  
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
  sp=prot.getSpacer(chainID[0]);
  inFile.close();
  
  unsigned int index1 = sp->getIndexFromPdbNumber(ID);
  unsigned int index2 = index1 + 1 ;
  
  
  ofstream scwrlOut(scwrlFile.c_str());
  string tmp = sSCWRLSequence(*sp, index1, index2, mut);
  scwrlOut << tmp << "\n";
}
