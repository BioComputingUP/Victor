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

// --*- C++ -*------x------------------------------------------------------------------------------
//
//
// Description:     This program calculates phylogenetic trees starting from FASTA file.
//					Choose from 2 different methos Neighbor joining and UPGMA.
//					This program calculates ClustalW allignment and use NJ for create guide trees.
// -----------------x------------------------------------------------------------------------------

#include <Alignment.h>
#include <AlignmentData.h>
#include <SequenceData.h>
#include <SubMatrix.h>
#include <NWAlign.h>
#include <AGPFunction.h>
#include <ScoringS2S.h>

#include <GetArg.h>

#include <NewickTree.h>
#include <PhyloSupport.h>
#include <ClustalW.h>

#include <string>
#include <vector>
#include <sstream>
#include <ctime>

using namespace Victor::Align2;
using namespace Victor::Phylo;
using namespace Victor;
using namespace std;



/// Show command line options and help text.

void
sShowHelp() {
    cout << "\nP Phylogenetic Tree GENERATOR" //This program read FASTA file and calculate phylogenetic tree. Use NWALIGN for align.
            << "\nThis program read FASTA file and calculate phylogenetic tree. Use NWALIGN for align."
            << "\nOptions:"
            << "\n"
            << "\n * [--in <name>]     \t Name of input FASTA file"
            << "\n   [--out <name>]    \t Name of output file (default = to screen)"
            << "\n   [-m <name>]       \t Name of substitution matrix file (default = blosum62.dat)"
			<< "\n"
            << "\n   [-d <double>]     \t Min. master vs all seq. id. filter threshold (suggested = 0.00)"
            << "\n   [-D <double>]     \t Min. all vs all seq. id. filter threshold (suggested = 0.00)"
            << "\n   [-u <double>]     \t Max. master vs all seq. id. filter threshold (suggested = 1.00)"
            << "\n   [-U <double>]     \t Max. all vs all seq. id. filter threshold (suggested = 1.00)"
			<< "\n"
            << "\n   [--ws <0|1|2|3>]  \t Weighting scheme for profiles (default = 0, i.e. no weighting scheme)"
            << "\n                     \t --ws=0: No weighting scheme (default)."
            << "\n                     \t --ws=1: Calculate a frequency profile or PSSM using Henikoff weighting scheme."
            << "\n                     \t --ws=2: Calculate a frequency profile or PSSM using PSIC weighting scheme."
            << "\n                     \t --ws=3: Calculate a frequency profile or PSSM using SeqDivergence weighting scheme."
			<< "\n"
			<< "\n   [-o <double>]     \t Open gap penalty (default = 12.00)"
            << "\n   [-e <double>]     \t Extension gap penalty (default = 3.00)"
            << "\n"
            << "\n   [--function <0|1|2>] \t CLuster function (default = 0, i.e. UPGMA)"
            << "\n                     \t --gf=0: UPGMA (default)."
            << "\n                     \t --gf=1: NJ."
			<< "\n                     \t --gf=2: ClustalW use NJ."
            << "\n   [-t <int>] 	   \t For decide dimension of token, very important for ClustalW(10%size of first seq)"
            << "\n   [--cSeq <double>] \t Coefficient for sequence alignment (default = 1.0)"
			<< "\n"
			<< "\n   [--ktuples]     	\t use ktuples method  formula for calculate distance of pairwise sequence"
            << "\n   [--verbose]       	\t Verbose mode"
            << "\n" << endl;
}


/// Show command line options and help text.

int main(int argc, char **argv) {

    string inputFileName, outputFileName, matrixFileName, matrixStrFileName;
    double openGapPenalty, extensionGapPenalty;
    unsigned int  function;
    double cSeq;
    double downs, downa, ups, upa;
    bool verbose=false;
    bool ktuples=false;
    unsigned int weightingScheme;
    int tokenSize;

    struct tm* newtime;

    time_t t;

    // --------------------------------------------------
    // 0. Treat options
    // --------------------------------------------------
    if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    }


    getArg("-in", inputFileName, argc, argv, "cyc1_input.fasta");
	getArg("-out", outputFileName, argc, argv, "!");

	getArg("m", matrixFileName, argc, argv, "blosum62.dat");

	getArg("-function", function, argc, argv, 0);

    getArg("-ws", weightingScheme, argc, argv, 0);

	getArg("o", openGapPenalty, argc, argv, 12.00);
	getArg("e", extensionGapPenalty, argc, argv, 3.00);

    getArg("d", downs, argc, argv, 999.9);
    getArg("D", downa, argc, argv, 999.9);
    getArg("u", ups, argc, argv, 999.9);
    getArg("U", upa, argc, argv, 999.9);

    getArg("t", tokenSize, argc, argv, -1);

	getArg("-cSeq", cSeq, argc, argv, 0.80);

	ktuples = getArg("-ktuples", argc, argv);
	verbose = getArg("-verbose", argc, argv);

    // --------------------------------------------------
    // 1. Load FASTA
    // --------------------------------------------------
	Alignment aliSec;
		string path = getenv("VICTOR_ROOT");
		if (path.length() < 3)
			cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;
		string dataPath = path + "data/";
		string finalPath=dataPath+inputFileName;

		cout<<"path "<<finalPath<<endl;
		ifstream secFile(finalPath.c_str());
				if (!secFile)
					ERROR("Error opening secondary structure FASTA file.", exception);
		aliSec.loadFasta(secFile);
		if (aliSec.size() < 1)
			ERROR("Secondary structure FASTA file must contain two sequences.", exception);

	// --------------------------------------------------
	// 2. Setup Param
	// --------------------------------------------------
		PhyloSupport::extensionGapPenalty=extensionGapPenalty;
		PhyloSupport::openGapPenalty=openGapPenalty;
		PhyloSupport::cSeq=cSeq;
		PhyloSupport::downs=downs;
		PhyloSupport::downa=downa;
		PhyloSupport::ups=downa;
		PhyloSupport::upa=downa;
		PhyloSupport::weightingScheme=weightingScheme;
		PhyloSupport::tokenSize=tokenSize;
		PhyloSupport::verbose=verbose;
		PhyloSupport::matrix=matrixFileName;

    // --------------------------------------------------
    // 3. Load data
    // --------------------------------------------------

	NewickTree tree;
	ClustalW cw;
	string out="Error Tree";
	if(function==0){
		cout<<"##########################################Starting PhyloTreeUPGMA#####################################"<<endl;
		time(&t);
		newtime = localtime(&t);
		cout <<"localTime:"<< newtime->tm_hour << " hour " << newtime->tm_min<<" min "<<newtime->tm_sec<<" sec "<< endl;
		tree.upgma(aliSec,ktuples,verbose);
		cout<<"##########################################Tree create with PhyloTreeUPGMA#############################"<<endl;
		time(&t);
		newtime = localtime(&t);
		cout <<"localTime:"<< newtime->tm_hour << " hour " << newtime->tm_min<<" min "<<newtime->tm_sec<<" sec "<< endl;
		out=tree.printNewickTree();
	}
	else if(function==1){
		cout<<"##########################################Starting PhyloTreeNJ########################################"<<endl;
		time(&t);
		newtime = localtime(&t);
		cout <<"localTime:"<< newtime->tm_hour << " hour " << newtime->tm_min<<" min "<<newtime->tm_sec<<" sec "<< endl;
		tree.neighborJoining(aliSec,ktuples,verbose);
		cout<<"##########################################Tree create with PhyloTreeNJ################################"<<endl;
		time(&t);
		newtime = localtime(&t);
		cout <<"localTime:"<< newtime->tm_hour << " hour " << newtime->tm_min<<" min "<<newtime->tm_sec<<" sec "<< endl;
		out=tree.printNewickTree();
	}
	else if(function==2){
		cout<<"##########################################Starting PhyloTreeNJ########################################"<<endl;
		time(&t);
		newtime = localtime(&t);
		cout <<"localTime:"<< newtime->tm_hour << " hour " << newtime->tm_min<<" min "<<newtime->tm_sec<<" sec "<< endl;
		tree.neighborJoining(aliSec,ktuples,verbose);
		cout<<"##########################################Guide Tree create with PhyloTreeNJ################################"<<endl;
		time(&t);
		newtime = localtime(&t);
		cout <<"localTime:"<< newtime->tm_hour << " hour " << newtime->tm_min<<" min "<<newtime->tm_sec<<" sec "<< endl;
		out=tree.printNewickTree();
		cout<<"##########################################Starting ClustalW###############################"<<endl;
		time(&t);
		newtime = localtime(&t);
		cout <<"localTime:"<< newtime->tm_hour << " hour " << newtime->tm_min<<" min "<<newtime->tm_sec<<" sec "<< endl;
		cout<<"number of leaf from main "<<tree.getNumberOfLeaf()<<endl;
		cw= ClustalW(tree);
		time(&t);
		newtime = localtime(&t);
		cout <<"localTime:"<< newtime->tm_hour << " hour " << newtime->tm_min<<" min "<<newtime->tm_sec<<" sec "<< endl;

	}
	else
	{
		cout<<"invalid cluster, select [--cluster <0|1>]";
	}

    // --------------------------------------------------
    // 2. Output Tree
    // --------------------------------------------------
	if(outputFileName=="!"){
		cout<<"Print the Tree "<<endl;
		cout<<out<<endl;
	}
	else{
		cout<<"Print the Tree in the file "<<outputFileName<<endl;
		//write tree in file
		ofstream outFile(outputFileName.c_str());
		if(!outFile) {
			cout<<"Error During creation File out.clustalw";
		}
		outFile<<out<<endl;

		outFile.close();

	}

	cout<<"------------------------End-------------------------"<<endl;

	return 0;

}
