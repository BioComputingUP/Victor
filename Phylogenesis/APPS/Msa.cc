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

#include <GetArg.h>
#include <SubMatrix.h>
#include <Align.h>
#include <SequenceData.h>
#include <ScoringScheme.h>
#include <GapFunction.h>
#include <AGPFunction.h>
#include <DistanceMatrix.h>
#include <UPGMA.h>
#include <NeighborJoining.h>
#include <NewickPhylogeneticTreeExporter.h>
#include <ClustalW.h>
#include <iostream>
#include <vector>

using namespace Victor::Align2;
using namespace Victor::Phylogenesis;
using namespace std;

void sShowHelp() {
	cout<<"\nMULTILPLE SEQUENCE ALIGNMENT GENERATOR"
			<< "\nThis program is an implementation of ClustalW (Thompson, Higgins & Gibson, 1994),"
			<< "\none variant of the progressive method for multiple sequence alignment."
			<< "\nGuide tree is builded with NJ (or UPGMA) and provided as output in the Newick tree format."
			<< "\nInput must be provided as FASTA format, it is possible to provide also the scoring matrix as input "
			<< "\nor it can be dynamically calculated using global, local or freeshift alignments."
			<<" \n"
			<< "\nOptions:"
			<< "\n * [--in <name>]     \t Name of input FASTA file"
			<< "\n   [--verbose]       \t Verbose mode"
			<< "\n"
			<< "\nDistance matrix:"
			<< "\n   [--global]        \t Needleman-Wunsch-Gotoh global alignment (default)"
			<< "\n   [--local]         \t Smith-Waterman local alignment"
			<< "\n   [--freeshift]     \t Free-shift alignment"
			<< "\n   [--matrix <name>] \t Custom distance matrix"
			<< "\n   [-m <name>]       \t Name of substitution matrix file (default = blosum30.dat)"
            << "\n   [-o <double>]     \t Open gap penalty (default = 10.00)"
            << "\n   [-e <double>]     \t Extension gap penalty (default = 0.10)"
			<< "\n"
			<< "\nClustering method:"
			<< "\n   [--upgma]         \t Unweighted Pair Group Method with Arithmetic Mean (default = Neighbor joining)"
			<< "\n"
			<< "\nClustalW options:"
			<< "\n   [--cwo <double>]  \t Open gap penalty (default = 10.00)"
			<< "\n   [--cwe <double>]  \t Extension gap penalty (default = 0.20)"
			<< "\n   [--cutoff <double>] \t Cutoff value for divergent sequences (default = 0.4)"
			<< "\n   [--pam ] 	\t PAM substitution matrix series (default = Blosum)"
			<< "\n"
			<< "\nOutput format:"
			<< "\n   [--outNewick <name>]\t Name of the Newick tree format output file (default = to screen)"
			<< "\n   [--outMsa <name>] \t Name of the MSA output FASTA file (default = to screen)"
			<< "\n   [--clustal]	\t MSA output in clustal format (default = FASTA)"
            << "\n\n";
}

int main( int argc, char* argv[] ){

	cout<<"VICTOR - Multiple sequence alignment generator\n\n";

	//
	// Init options
	//

	string inputFileName, matrixFileName, cwMatrixFileName, scoringMatrixFileName, newickOutputFileName, msaOutputFileName;
	double openGapPenalty, extensionGapPenalty, cwOpenGapPenalty, cwExtensionGapPenalty, cutoff;
	bool verbose, upgma, pam, clustal;

	if (getArg("h", argc, argv) || argc == 1) {
		sShowHelp();
		return 1;
	}

	getArg("-in", inputFileName, argc, argv, "!");
	getArg("-matrix", scoringMatrixFileName, argc, argv, "!");
	getArg("-outNewick", newickOutputFileName, argc, argv, "!");
	getArg("-outMsa", msaOutputFileName, argc, argv, "!");

	AlignAlgorithm alignAlgorithm = global;
	if(getArg("-local", argc, argv))
		alignAlgorithm = local;
	if(getArg("-freeshift", argc, argv))
		alignAlgorithm = freeshift;

	getArg("m", matrixFileName, argc, argv, "blosum30.dat");
    getArg("o", openGapPenalty, argc, argv, 10.0);
    getArg("e", extensionGapPenalty, argc, argv, 0.1);

    getArg("-cwo", cwOpenGapPenalty, argc, argv, 10.00);
    getArg("-cwe", cwExtensionGapPenalty, argc, argv, 0.2);
    getArg("-cutoff", cutoff, argc, argv, 0.4);

    upgma = getArg("-upgma", argc, argv);
    pam =  getArg("-pam", argc, argv);
    clustal = getArg("-clustal", argc, argv);
    verbose = getArg("-verbose", argc, argv);

    //
    // Load data
    //

    string path = getenv("VICTOR_ROOT");
    if (path.length() < 3)
        cout << "Warning: environment variable VICTOR_ROOT is not set." << endl;
    string dataPath = path + "data/";

    vector<Sequence*> sequences;
    if (inputFileName != "!") {
        ifstream inputFile(inputFileName.c_str());
        if (!inputFile)
            ERROR("Error opening input FASTA file.", exception);
        sequences = Sequence::LoadFasta(inputFile);
        if (sequences.size() < 1)
            ERROR("Input FASTA file must contain at least two sequences.", exception);
    } else
        ERROR("ClustalW needs input FASTA file.", exception);

    matrixFileName = dataPath + matrixFileName;
    ifstream matrixFile(matrixFileName.c_str());
    if (!matrixFile)
        ERROR("Error opening substitution matrix file.", exception);
    SubMatrix* sub = new SubMatrix(matrixFile);

    //
	// Distance matrix
    //

    DistanceMatrix* dm;
	if(scoringMatrixFileName == "!"){
		dm = new DistanceMatrix(sequences, alignAlgorithm, sub, openGapPenalty, extensionGapPenalty);
	}
	else{
		dm = new DistanceMatrix(scoringMatrixFileName);
	}

    //
    // Clustering
    //

    HierarchicalClusterBuilder* builder;
    if(upgma){
    	builder = new UPGMA(sequences, dm);
    }
    else{
    	builder = new NeighborJoining(sequences, dm);
    }

	//
	// ClustalW
	//

	ClustalW* cw = new ClustalW(dm, builder, sequences, cutoff, cwOpenGapPenalty, cwExtensionGapPenalty, pam, verbose);
	MultipleAlignment* result = cw->performMultipleSequenceAlignment();

    //
    // Newick
    //

    PhylogeneticTreeExporter* expoter;
	if(newickOutputFileName == "!"){
		expoter = new NewickPhylogeneticTreeExporter();
		cout<<"\nGuide tree:\n";
	}
	else{
		expoter = new NewickPhylogeneticTreeExporter(newickOutputFileName);
	}
	expoter->save(cw->getGuideTree());

	//
	// Output MSA
	//

	if(msaOutputFileName != "!"){
		cout<<"\nWriting multiple alignment file "<<msaOutputFileName<<endl;
		std::ofstream  output;
		output.open(msaOutputFileName.c_str());
		clustal ? result->saveClustal(output) : result->saveFasta(output);
		output.close();
		return 0;
	}
	cout<<"\nResult:\n"<<endl;
	clustal ? result->saveClustal(cout) : result->saveFasta(cout);
}


