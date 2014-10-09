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
 @Description This program give information about the torsion angle  of a given protein structural model. 
 */
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <AminoAcid.h>
#include <String2Number.h>
using namespace Victor::Biopool; 
using namespace Victor;

void sShowHelp(){
  cout << "pdb2tor - Pdb Torsion Angles\n\n"
       << "This program give information about the torsion angle" 
       << " of a given protein structural model. \n"
       << "Contact the author for further details.\n"
       << "Developed by Silvio Tosatto <silvio.tosatto@unipd.it> \n\n"
       << "Input options: \n"
       << "\t-i <filename> \t\t Input single pdb file.\n"
       << "\t-c <id>  \t ID of chain to load from PDB file.\n"
       << "\t[--nmr] \t\t Calculate average over NMR ensemble ( only for NMR model ).\n"
       << "\t\t\t\t DEFAULT: The first chain and/or the first model.\n\n"
       << "\t-I <filelist> \t\t Input filelist (a file containing the list of PDB file)." <<"\n"
       << "\t\t\t\t For major information upon the format of the input filelist and the option\n"
       << "\t\t\t\t type: pdb2tor --filelist \n"
       << "Output option: \n"
       << "\t -P \t\t\t Give per residue phi and psi angles.\n"
       << "\t -O \t\t\t Give per residue phi, psi and omega angles.\n"
       << "\t -C \t\t\t Give per residue phi, psi, omega and chi angles.\n"
       << "\t -A \t\t\t Give per residue phi, psi, omega, chi, pre-psi and pre-psi angle.\n"
       << "\t -r \t\t\t Output in format for data TorsionPotential.\n"
       << "\n";
}
void sSpecification()
{
  cout <<"\n\nINPUT FILELIST:\n"
       <<"The file list must contain the list of pdb file to examinate or their path.\n"
       <<"The file must be organized in column: the first containing the pdb file path\n"
       <<"and the second the Chain ( optional ).\n"
       <<"If you want to consider the chain specified in the second column you have to\n"
       <<"use the option --complete. If not specified the second column will be ignored\n"
       <<"and the programm will consider the first chain.\n"
       <<"If you use a filelist of nmr structure use the option --nmr for consider all models.\n\n"
       <<"FILELIST FORMAT EXAMPLE:\n\n"
       <<"PATH_OF_PDB_STRUCTURE1\tCHAIN\n"
       <<"PATH_OF_PDB_STRUCTURE2\tCHAIN\n"
       <<"....\n\n";
}

int main(int nArgs, char* argv[])
{
  if (getArg( "h", nArgs, argv)) {
      sShowHelp();
      return 1;
    };
  if (getArg( "-filelist", nArgs, argv)) {
      sSpecification();
      return 1;
    }
  
  string inputFile, inputFilelist, chainID;
  vector<char> allCh; 
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "I", inputFilelist, nArgs, argv, "!");
  getArg( "c", chainID, nArgs, argv, " ");
  bool nmr = getArg( "-nmr", nArgs, argv);
  bool complete = getArg( "-complete", nArgs, argv);
  bool p = getArg( "P", nArgs, argv);
  bool o = getArg( "O", nArgs, argv);
  bool c = getArg( "C", nArgs, argv);
  bool a = getArg( "A", nArgs, argv);
  bool verbose = getArg( "r", nArgs, argv);

  //Control for ChaiId option
  if ((chainID != " ") && (complete))
    ERROR("You are using the 'chainID' option and the 'complete' option at the same time.", exception);

  //Control for input file
  if ((inputFile == "!") && (inputFilelist == "!")) {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  
  if ((inputFile != "!") && (inputFilelist != "!")) {
      cout << "Please choose between filelist and file mode. Aborting. "
	   << "(-h for help)" << endl;
      return -2;     
    }
  if (( !a ) && ( !c ) && ( !o ) && ( !p ) && ( !verbose ))
    ERROR("Choose a valid option for output.", exception);

  ifstream inFile(inputFilelist.c_str());
  
  if ((!inFile) && (inputFilelist != "!"))
    ERROR("File not found.", exception);

  int totalanalized = 0;
  int totalmodel = 0;

  while ((inFile) || (inputFile != "!")) {
        if (inputFilelist != "!"){
            inFile >> inputFile;

            if (!inFile)
              break;
            if (complete) {
                int check = checkForBlankLine(inFile, true);
                cout<<"check value"<<check;
                if ( check == 1 )
                  chainID = " ";
                else
                  inFile >> chainID;
                skipToNewLine(inFile);
              }
	}
      
        if ( inputFile == "!" ){
            cout << "Missing file specification. Aborting. (-h for help)" << endl;
            return -1;
          }
      
        ifstream inFile2(inputFile.c_str());
        if (!inFile2)
          ERROR("File not found.", exception);

        PdbLoader pl(inFile2);
        pl.setNoHAtoms();
        pl.setNoVerbose();

        allCh = pl.getAllChains(); 


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
        else{
            chainID=allCh[0];
        }

        
        
        pl.setPermissive();
        Protein prot;
        prot.load(pl);

        unsigned int max;

        if (!pl.isValid()){
            if (!verbose)
              cout << "Warning: Invalid PDB file found "<<inputFile<<".\n";
            inputFile = "!";
            continue;
          }

        if (!nmr)
          max = 1;
        else{
            max = pl.getMaxModels();
            if ( max == 0)
              {
                max = 1;
                if ( !verbose )
                  cout <<"Warning: the file "<<inputFile<<" probably is not an nmr structure.\n";
              }
          }

        totalanalized += 1;

        for (unsigned int j = 1; j <= max; j++){
            pl.setModel(j);
            Spacer *sp;
            sp=prot.getSpacer(chainID[0]);
            if (!pl.isValid())  {
                if ( !verbose )
                  cout << "Warning: Invalid PDB file found:"<<inputFile<<".\n";
                if ( j == max )
                  {
                    inputFile = "!";
                    totalanalized -= 1;
                  }
                continue;
              }
            totalmodel += 1;
            int n , i;
            n = sp->sizeAmino()-1;
            cout.setf(ios::fixed, ios::floatfield);

            if (!verbose) {
                cout <<"PDB file: "<<inputFile<<"\tChain: "<<chainID<<"\n"; 
                cout <<"Type\tNumber\t";
                if ( p )
                  cout<<" phi\t\tpsi\n";
                if ( o )
                  cout<<" phi\t\tpsi\t\tomega\n";
                if ( c )
                  cout<<" phi\t\tpsi\t\tomega\t\tchi1\t\tchi2\n";
                if ( a )
                  cout<<"pre-phi\t\tpre-psi\t\tphi\t\tpsi\t\tomega\t\tchi1\t\tchi2\n";

                for ( i = 1 ; i < n; i++ ){
                    if ( p ) {
                        cout <<sp->getAmino(i).getType1L()<<"\t"<<i+1<<"\t"
                             <<sp->getAmino(i).getPhi()<<"\t"
                             <<sp->getAmino(i).getPsi()<<"\n";
                      }
                    if ( o )   {
                        cout <<sp->getAmino(i).getType1L()<<"\t"<<i+1<<"\t"
                             <<sp->getAmino(i).getPhi()<<"\t"
                             <<sp->getAmino(i).getPsi()<<"\t"
                             <<sp->getAmino(i).getOmega()<<"\n";
                      }
                    if ( c )   {
                        cout <<sp->getAmino(i).getType1L()<<"\t"<<i+1<<"\t"
                             <<sp->getAmino(i).getPhi()<<"\t"
                             <<sp->getAmino(i).getPsi()<<"\t"
                             <<sp->getAmino(i).getOmega()<<"\t"; 
                        if (sp->getAmino(i).getMaxChi() != 0 ){
                            if ( sp->getAmino(i).getMaxChi() == 1 )
                              cout <<sp->getAmino(i).getChi(0)<<"\n";
                            if ( sp->getAmino(i).getMaxChi() != 1 )
                              {
                                cout<<sp->getAmino(i).getChi(0)<<"\t"
                                    <<sp->getAmino(i).getChi(1)<<"\n";
                              }
                          }
                        else
                          cout <<"\n";
                      }

                    if ( a )  {
                        if ( i == 1 ){

                            cout <<sp->getAmino(i).getType1L()<<"\t"<<i+1<<"\t";
                            cout <<"NOT\t\tNOT\t\t"
                                 <<sp->getAmino(i).getPhi()<<"\t"
                                 <<sp->getAmino(i).getPsi()<<"\t"
                                 <<sp->getAmino(i).getOmega()<<"\t"; 
                            if (sp->getAmino(i).getMaxChi() != 0 )   {
                                if ( sp->getAmino(i).getMaxChi() == 1 )
                                  cout <<sp->getAmino(i).getChi(0)<<"\n";
                                if ( sp->getAmino(i).getMaxChi() != 1 ){
                                    cout<<sp->getAmino(i).getChi(0)<<"\t"
                                        <<sp->getAmino(i).getChi(1)<<"\n";
                                  }
                              }
                            else
                              cout <<"\n";
                          }
                        else{
                            cout <<sp->getAmino(i).getType1L()<<"\t"<<i+1<<"\t";
                            cout <<sp->getAmino(i-1).getPhi()<<"\t"
                                 <<sp->getAmino(i-1).getPsi()<<"\t"
                                 <<sp->getAmino(i).getPhi()<<"\t"
                                 <<sp->getAmino(i).getPsi()<<"\t"
                                 <<sp->getAmino(i).getOmega()<<"\t"; 
                            if (sp->getAmino(i).getMaxChi() != 0 )   {
                                if ( sp->getAmino(i).getMaxChi() == 1 )
                                  cout <<sp->getAmino(i).getChi(0)<<"\n";
                                if ( sp->getAmino(i).getMaxChi() != 1 )	{
                                    cout<<sp->getAmino(i).getChi(0)<<"\t"
                                        <<sp->getAmino(i).getChi(1)<<"\n";
                                  }
                              }
                            else
                              cout <<"\n";
                          }
                      }
                  }
              }
            else if ( verbose )  {
                 for ( i = 1 ; i < n; i++ ) {
                     cout <<"\t";
                     cout << sp->getAmino(i).getPhi() << "\t" << sp->getAmino(i).getPsi() << "\t"
                          << sp->getAmino(i).getType() << "\t" 
                          << sp->getAmino(i-1).getPhi() <<"\t" 
                          << sp->getAmino(i-1).getPsi()<<"\t"
                          << sp->getAmino(i).getOmega() <<"\t" <<  sp->getAmino(i).getMaxChi()<<"\t";
                    for (unsigned int k = 0; k < sp->getAmino(i).getMaxChi(); k++)
                      cout << sp->getAmino(i).getChi(k) << "\t";

                     cout << endl;
                   }
              }
          }
        // reset variable to trigger break condition:
        inputFile = "!";
    }
  if ( !verbose )  {
      cout <<"Total file analized:"<<"\t"<<totalanalized<<"\n";
      if ( nmr )
	cout <<"Total nmr model analized:\t"<<totalmodel<<"\n";
    }
}
