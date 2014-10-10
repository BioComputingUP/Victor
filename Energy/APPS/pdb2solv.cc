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
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <SolvationPotential.h>
#include <PolarSolvationPotential.h>

using namespace Victor;

using namespace Victor::Energy;
using namespace Victor::Biopool;
void sShowHelp(){
  cout << "Pdb to Solvation Energy\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file for PDB template\n"
       << "\t-o <filename> \t\t Output file for solvation distribution\n"
       << "\t-c <id>  \t ID of chain to load from PDB file\n"
       << "\t[-p] \t\t\t Polarity check mode\n"
       << "\t[-v] \t\t\t Verbose mode\n"
       << "\t[-m] \t\t\t Check mode\n"
       << "\n";
}


int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))  {
      sShowHelp();
      return 1;
    };
  vector<char> allCh;
  string inputFile, outputFile,chainID;
  bool verbose = getArg( "v", nArgs, argv);
  bool check = getArg( "m", nArgs, argv);
  bool polar = getArg( "p", nArgs, argv);
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "test.out");
  getArg( "c", chainID, nArgs, argv, "!");
  if (inputFile == "!")  {
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
  
  if ((verbose) || (check))
    cout << "Input: " << inputFile << " \t Output: " << outputFile << "\n";

  if (verbose)  {
      cout << "loaded." << endl;
      cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
      cout << "Solvation statistics:\n";
    }

  vector<vector<int> > sum;
  vector<vector<int> > sumPolar;
  vector<int> tmp;
  for (unsigned int i = 0; i < AminoAcid_CODE_SIZE; i++)
    sum.push_back(tmp);
  for (unsigned int i = 0; i < AminoAcid_CODE_SIZE; i++)
    sumPolar.push_back(tmp);

  for (unsigned int i = 0; i < sp->sizeAmino(); i++) {
      AminoAcid aa = sp->getAmino(i);

      if (!aa.getSideChain().isMember(CB)){ 
	  if (verbose)
	    cout << "\t X";
	}
      else{
	  int count = 0;
	  int count2 = 0;
	  
	  for (unsigned int j = 0; j < i; j++)
	    if (sp->getAmino(j).getSideChain().isMember(CB))
	      if (aa.getSideChain()[CB].distance( sp->getAmino(j).getSideChain()[CB]) 
		       <= (polar ? SOLVATION_CUTOFF_DISTANCE_POLAR: SOLVATION_CUTOFF_DISTANCE)){
		  if (isPolar((AminoAcidCode)(sp->getAmino(j).getCode())) )
		    count++;
		  else
		    count2++;
		}

	  for (unsigned int j = i+1; j < sp->sizeAmino(); j++)
	    if (sp->getAmino(j).getSideChain().isMember(CB))
	      if (aa.getSideChain()[CB].distance(
		     sp->getAmino(j).getSideChain()[CB]) 
		       <= (polar ? SOLVATION_CUTOFF_DISTANCE_POLAR 
			: SOLVATION_CUTOFF_DISTANCE)){
		  if (isPolar((AminoAcidCode)sp->getAmino(j).getCode()))
		    count++;
		  else
		    count2++;
		}

	  if (verbose)  {
	      cout << "\t " << (count+count2);
	      if (polar)
		cout << " " << count;
	    }

	  sum[aa.getCode()].push_back(count+count2);
	  sumPolar[aa.getCode()].push_back(count);
	}

      if (verbose)
	if ((i + 1 ) % 8 == 0)
	  cout << "\n";
    }
  if (verbose)   {
      cout << "\n";
      cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
    }

  ofstream outFile(outputFile.c_str());
  if (!outFile)
    ERROR("File not found.", exception);

  for (unsigned int i = 0; i < sum.size()-1; i++)   {
      if (polar)
	outFile << "P  ";
      outFile << setw(3) << sumPolar[i].size() << "  " 
	      << aminoAcidThreeLetterTranslator(static_cast<AminoAcidCode>(i)) 
	      << "   "; 
      if (polar){
	  for (unsigned int j = 0; j < sum[i].size(); j++)
	    outFile << "\t" << sum[i][j] << " " << sumPolar[i][j];
	}
      else{
	  for (unsigned int j = 0; j < sum[i].size(); j++)
	    outFile << "  " << sum[i][j];
	}
      outFile << "\n";
    }

}
