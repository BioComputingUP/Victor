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
 @Description FRST - Function of Rapdf, Solvation and Torsion potentials
 *This program calculates a pseudo-energy to evaluate the quality "
	*of a given protein structural model, as expressed in a single 
 *(real) number. 
 */
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <SolvationPotential.h>
#include <EffectiveSolvationPotential.h>
#include <RapdfPotential.h>
#include <TorsionPotential.h>
#include <PhiPsi.h>
#include <PhiPsiOmegaChi1Chi2PreAngle.h>

using namespace Victor::Biopool;
using namespace Victor::Energy;
using namespace Victor;

const long double W_RAPDF = 2.5;
const long double W_SOLV = 500.0;
const long double W_HYDB = -50.0;
const long double W_TORS = 350.0;

const long double MAX_AA = -100.0;
const long double MIN_AA = 50.0;

void sShowHelp(){
	cout << "FRST - Function of Rapdf, Solvation and Torsion potentials\n\n"
	<< "This program calculates a pseudo-energy to evaluate the quality "
	<< "of a given protein structural model, as expressed in a single "
	<< "(real) number. Contact the author for further details.\n"
	<< "Developed by Silvio Tosatto <silvio.tosatto@unipd.it> \n"
	<< "   Options: \n"
	<< "\t-I <filelist> \t\t Input *filelist* file for PDB templates\n"
	<< "\t-i <filename> \t\t Input file for PDB template\n"
        << "\t-c <id>  \t\t ID of chain to load from PDB file\n"
	<< "\t[-T] \t\t\t New (Mk2) torsion angle potential \n"
	<< "\t[-S] \t\t\t New (Mk2) solvation potential \n"
	<< "\t[-s <resol>] \t\t\t Solvation potential resolution (default = 1) \n"
	<< "\t             \t\t\t Attention: must be a clean divisor of 30. \n"
	<< "\t[-v] \t\t\t Verbose mode\n"
	<< "\t[--casp] \t\t\t CASP QA mode 1\n"
	<< "\t[-C <# residues>] \t\t\t CASP QA mode2 , n. residues to assume\n"
	<< "\t[-p] \t\t\t Per residue energy calculation\n"
	
	<< "\n"
	<< "\tRemember to set VICTOR_ROOT enviroment variable" <<endl<<endl;
}

double sHydrogen(Spacer& sp){
	int count = 0;
	
	for (int i = 0; i <(int)sp.sizeAmino(); i++)
		for (int j = 0; j < (int)sp.sizeAmino(); j++)
			if (abs(i-j) > 1){
				double dist = sp.getAmino(i)[N].distance(sp.getAmino(j)[O]);
				double dist2 = sp.getAmino(i)[N].distance(sp.getAmino(j)[C]);
				double dist3 = sp.getAmino(i)[CA].distance(sp.getAmino(j)[O]);
				
				if ((dist <= 4.0) && (dist >= 2.0) && (dist < dist2) && (dist < dist3)){
					count++;
					break;
				}
			}
				return count;
}


int main(int nArgs, char* argv[]){ 
        char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
	
	if (getArg( "h", nArgs, argv)) {
		sShowHelp();
		return 1;
        };
	string victorRoot =getenv("VICTOR_ROOT");
        cout<<victorRoot;
	string inputFile, inputFilelist, solvOut, rapdfOut, torOut, chainID;
	bool newTorsion = getArg( "T", nArgs, argv);
	bool newSolvation = getArg( "S", nArgs, argv);
	bool verbose = getArg( "v", nArgs, argv);
	
	bool casp = getArg( "-casp", nArgs, argv);
	unsigned int caspAll;
	getArg( "C", caspAll, nArgs, argv, 9999);
	if (caspAll < 9999)
		casp = true;
	
	bool perResidue = getArg( "p", nArgs, argv);
	getArg( "i", inputFile, nArgs, argv, "!");
	getArg( "I", inputFilelist, nArgs, argv, "!");
	getArg( "c", chainID, nArgs, argv, "!");
	
	unsigned int solvRes;
	getArg( "s", solvRes, nArgs, argv, 1);
	
	if ((inputFile == "!") && (inputFilelist == "!"))    {
		cout << "Missing file specification. Aborting. (-h for help)" << endl;
		return -1;
        }
	
	if ((inputFile != "!") && (inputFilelist != "!")){
		cout << "Please choose between filelist and file mode. Aborting. "
	    << "(-h for help)" << endl;
		return -2;     
        }
	
	Potential* solv = NULL;
	if (!newSolvation)
		solv = new SolvationPotential(solvRes);
	else
		solv = new EffectiveSolvationPotential;
	RapdfPotential rapdf;
	TorsionPotential* tors = NULL;
	if (!newTorsion)
		tors =  new PhiPsi(10);//default ARCSTEP = 10
	else
		tors =  new PhiPsiOmegaChi1Chi2PreAngle(20);//ARCSTEP 20, ARCSTEP2 40 
		
	ifstream inFile(inputFilelist.c_str());
	if ((!inFile) && (inputFilelist != "!"))
		ERROR("File not found.", exception);
		
	while ((inFile) || (inputFile != "!"))		{
		if (inputFilelist != "!") {
			inFile >> inputFile;
			if (!inFile)
				break;
		}
		ifstream inFile2(inputFile.c_str());
		if (!inFile2)
			ERROR("File not found.", exception);
		PdbLoader pl(inFile2);
		Spacer *sp; 
                pl.setNoHAtoms();
                vector<char> allCh;
                allCh = pl.getAllChains(); 
                for (unsigned int i = 0; i < allCh.size(); i++)
                  cout << "\t," << allCh[i] << ",";
                cout << "\n";

                /*check on validity of chain: 
                if user select a chain then check validity
                 else select firs valid one by default*/
                if (chainID != "!") {
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
                    chainID[0]=allCh[0];
                      cout << "Using chain " << chainID <<"\n";
                }
                Protein prot;
                prot.load(pl);

                sp=prot.getSpacer(chainID[0]);
		if (!pl.isValid()){
                    cerr << "Warning: Invalid PDB file found.\n";
                    inputFile = "!";
                    continue;
                }
			
		long double res = rapdf.calculateEnergy(*sp);
		long double sres = solv->calculateEnergy(*sp);
		long double hres = - 1.0 *  sHydrogen(*sp);
		long double tres = tors->calculateEnergy(*sp);
		long double sum = W_RAPDF * res + W_SOLV * sres + W_HYDB * hres 
				+ W_TORS * tres;
		cout.setf(ios::fixed, ios::floatfield);
		double divisor = (perResidue == true) ? sp->sizeAmino() : 1;
		if ((verbose) || (casp)){
                    cout << inputFile << "\t";
		}
                if (casp) {
                    double num = sp->sizeAmino();
                    if (caspAll < 9999)
                            num = sp->sizeAmino() + ( 0.5 * fabs(caspAll-sp->sizeAmino()));
                            double gdt = (sum - (MIN_AA * num)) / ((MAX_AA - MIN_AA) * num); 
                            if (gdt < 0.1)
                                    gdt = 0.1;
                            else if (gdt > 1.0)
                                    gdt = 1.0;
                            cout << setw(9) << setprecision(4) << (gdt / divisor);
                            if (caspAll < 9999){
				for (unsigned int index = 1; index <= caspAll; index++){
					if (sp->isGap(index)){
						cout << " " << setw(4) << setprecision(1) << "X";
					}
					else{
                                            unsigned int i = sp->getIndexFromPdbNumber(index);	
                                            // calculate partial energy:	
                                            long double res = 0.0;
                                            long double sres = 0.0;
                                            long double tres = 0.0;
                                            double num = 0.0;
                                            for (unsigned int j = 0; j < 20; j++){
						double d = 20.0 - j;	
						if (i >= j){
						    res += d * rapdf.calculateEnergy(sp->getAmino(i-j), *sp);
                                                    sres += d * solv->calculateEnergy(sp->getAmino(i-j), *sp);
                                                    tres += d * tors->calculateEnergy(sp->getAmino(i-j), *sp);
                                                    num += d * 1.0;
						}
						if (i+j < sp->sizeAmino()){
                                                    res += d * rapdf.calculateEnergy(sp->getAmino(i+j), *sp);
                                                    sres += d * solv->calculateEnergy(sp->getAmino(i+j),   *sp);
                                                    tres += d * tors->calculateEnergy(sp->getAmino(i+j), *sp);
                                                    num += d * 1.0;
						}
                                            }
                                            long double sum = W_RAPDF * res + W_SOLV * sres + W_TORS * tres;
                                            double caErr = (1.0 - sqrt( (sum - ( num * 1.5 * MIN_AA)) / ( num * 1.5 * (MAX_AA - ( 1.0 * MIN_AA)) ) ) ) * 20.0 + 1.0; 
                                            if (caErr < 1.0)
						caErr = 1.0;
                                            else if (caErr > 20.0)
						caErr = 20.0;
                                            cout << " " << setw(4) << setprecision(1) << caErr;
                                        }
					if (((index+1) % 50) == 0)
						cout << "\n";
                                    }
				}
				
		}
		else 
                    cout << setw(9) << setprecision(4) << (sum / divisor);
		if (verbose){
			cout << "\t" << setw(9) << setprecision(4) << (res / divisor) << "\t" 
			<< setw(9) << setprecision(4) << (sres / divisor) << "\t" 
			<< setw(9) << setprecision(4) << (hres / divisor) << "\t" 
			<< setw(9) << setprecision(4) << (tres / divisor);
		}
		cout << endl;
		// reset variable to trigger break condition:
                inputFile = "!";
	}
	delete tors;
}
