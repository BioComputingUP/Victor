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
 * @Description his program analyzes the backbone geometry of a PDB file in
 * terms of bond lengths and bond angles.
 */
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
#include <StatTools.h>

using namespace Biopool;

void sShowHelp(){
  cout << "Backbone Analyzer\n"
       << "This program analyzes the backbone geometry of a PDB file in"
       << "terms of bond lengths and bond angles.\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file\n"
       << "\t-c <chain> \t\t chain ID\n"
       << "\t[-s <n-start>] \t\t Aminoacid where analysis starts\n" 
       << "\t[-e <n-end>] \t\t Aminoacid where analysis ends\n"
       << "\t[--verbose] \t\t Verbose mode\n"
       << endl;
}

int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))    {
      sShowHelp();
      return 1;
    };

  string inputFile, chainID;
  vector<char> allCh; 
  unsigned int pdbIndex1, pdbIndex2;

  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "c", chainID, nArgs, argv, " "); 
  getArg( "s", pdbIndex1, nArgs, argv, 1);
  getArg( "e", pdbIndex2, nArgs, argv, 9999);

  bool verbose = getArg( "-verbose", nArgs, argv);

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
          chainID =allCh[0];
    
    
      pl.setPermissive();
      Protein prot;
      prot.load(pl);
      sp=prot.getSpacer(chainID[0]);
  pdbIndex2 = (pdbIndex2 >= sp->sizeAmino()-1) ? sp->sizeAmino()-2 : pdbIndex2;

  cout << "loaded.\n" << endl;

  cout << "-------------------------------------------------------\n";
  cout << "\t      Bond Lengths\t\t\tBond Angles\n";
  cout << "Num\t N->CA\t CA->C'\t C'->N\t\t N->CA\t CA->C'\t C'->N\n";
  cout << "-------------------------------------------------------\n";
  IntCoordConverter ic;

  vector<double> bl_n_ca;
  vector<double> bl_ca_c;
  vector<double> bl_c_n;
  vector<double> ba_n_ca;
  vector<double> ba_ca_c;
  vector<double> ba_c_n;

  for (unsigned int i = pdbIndex1; i <= pdbIndex2; i++)    {
      double tmp1 = ic.getBondLength(sp->getAmino(i)[N], sp->getAmino(i)[CA]);
      double tmp2 = ic.getBondLength(sp->getAmino(i)[CA], sp->getAmino(i)[C]);
      double tmp3 = ic.getBondLength(sp->getAmino(i)[C], sp->getAmino(i+1)[N]);

      double tmp4 = RAD2DEG * ic.getBondAngle(sp->getAmino(i-1)[C], 
				    sp->getAmino(i)[N], sp->getAmino(i)[CA]);
      double tmp5 = RAD2DEG * ic.getBondAngle(sp->getAmino(i)[N], 
				    sp->getAmino(i)[CA], sp->getAmino(i)[C]);
      double tmp6 = RAD2DEG * ic.getBondAngle(sp->getAmino(i)[CA],
                                    sp->getAmino(i)[C], sp->getAmino(i+1)[N]);

      bl_n_ca.push_back(tmp1);
      bl_ca_c.push_back(tmp2);
      bl_c_n.push_back(tmp3);

      ba_n_ca.push_back(tmp4);
      ba_ca_c.push_back(tmp5);
      ba_c_n.push_back(tmp6);

      if (verbose)
	cout << i << "\t" << setw(6) << setprecision(4) << tmp1
	          << "\t" << setw(6) << setprecision(4) << tmp2
	          << "\t" << setw(6) << setprecision(4) << tmp3 
	        << "\t\t" << setw(6) << setprecision(2) << tmp4
	          << "\t" << setw(6) << setprecision(2) << tmp5
	          << "\t" << setw(6) << setprecision(2) << tmp6	<< endl;
    }
  
  if (verbose)
    cout << "-------------------------------------------------------\n";

  unsigned int oldPrec = cout.precision();
  ios::fmtflags oldFlags = cout.flags();
  cout.setf(ios::fixed, ios::floatfield);

  cout << "Min:" 
       << "\t" << setw(6) << setprecision(4) << minimum(bl_n_ca) 
       << "\t" << setw(6) << setprecision(4) << minimum(bl_ca_c)
       << "\t" << setw(6) << setprecision(4) << minimum(bl_c_n)
       << "\t" 
       << "\t" << setw(6) << setprecision(2) << minimum(ba_n_ca) 
       << "\t" << setw(6) << setprecision(2) << minimum(ba_ca_c)
       << "\t" << setw(6) << setprecision(2) << minimum(ba_c_n)
       << "\n";
  cout << "Max:"
       << "\t" << setw(6) << setprecision(4) << maximum(bl_n_ca) 
       << "\t" << setw(6) << setprecision(4) << maximum(bl_ca_c)
       << "\t" << setw(6) << setprecision(4) << maximum(bl_c_n)
       << "\t" 
       << "\t" << setw(6) << setprecision(2) << maximum(ba_n_ca) 
       << "\t" << setw(6) << setprecision(2) << maximum(ba_ca_c)
       << "\t" << setw(6) << setprecision(2) << maximum(ba_c_n)
       << "\n";
  cout << "-------------------------------------------------------\n";
  cout << "Avg:" 
       << "\t" << setw(6) << setprecision(4) << average(bl_n_ca) 
       << "\t" << setw(6) << setprecision(4) << average(bl_ca_c)
       << "\t" << setw(6) << setprecision(4) << average(bl_c_n)
       << "\t" 
       << "\t" << setw(6) << setprecision(2) << average(ba_n_ca) 
       << "\t" << setw(6) << setprecision(2) << average(ba_ca_c)
       << "\t" << setw(6) << setprecision(2) << average(ba_c_n)
       << "\n";
  cout << "SD:" 
       << "\t" << setw(6) << setprecision(4) << standardDeviation(bl_n_ca) 
       << "\t" << setw(6) << setprecision(4) << standardDeviation(bl_ca_c)
       << "\t" << setw(6) << setprecision(4) << standardDeviation(bl_c_n)
       << "\t" 
       << "\t" << setw(6) << setprecision(2) << standardDeviation(ba_n_ca) 
       << "\t" << setw(6) << setprecision(2) << standardDeviation(ba_ca_c)
       << "\t" << setw(6) << setprecision(2) << standardDeviation(ba_c_n)
       << "\n";
  cout << "-------------------------------------------------------\n";

  cout.precision(oldPrec);
  cout.flags(oldFlags);

  return 0;
}
