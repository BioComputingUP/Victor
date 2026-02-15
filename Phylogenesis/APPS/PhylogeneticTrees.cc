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
#include <GapFunction.h>
#include <AGPFunction.h>
#include <DistanceMatrix.h>
#include <UPGMA.h>
#include <NeighborJoining.h>
#include <NewickPhylogeneticTreeExporter.h>
#include <ClustalW.h>
#include <SingleCluster.h>
#include <iostream>
#include <vector>

using namespace Victor;
using namespace Victor::Align2;
using namespace Victor::Phylogenesis;

using namespace std;

void sShowHelp() {
	cout<<"\nPHYLOGENETIC TREES GENERATOR"
			<< "\nThis program calculates phylogenetic trees for the sequences provided as FASTA input."
			<< "\nIt is possible to provide the scoring matrix as input or it can be dynamically calculated using global,"
			<<" \nlocal or freeshift alignments."
			<< "\nSupported algorithms are UPGMA and NJ, the tree is generated in the Newick tree format."
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
			<< "\nOutput format:"
			<< "\n   [--out <name>]    \t Name of the Newick tree format output file (default = to screen)"
			<< "\n   [--rooted] 	   \t Reroot the tree with 'mid-point' method (NJ only) "
			<< "\n\n";
}

int main( int argc, char* argv[] ){

	cout<<"VICTOR - Phylogenetic trees generator\n\n";

	//
	// Init options
	//

	string inputFileName, matrixFileName, scoringMatrixFileName, outputFileName;
	double openGapPenalty, extensionGapPenalty;
	bool verbose, upgma, rooted;

	if (getArg("h", argc, argv) || argc == 1) {
		sShowHelp();
		return 1;
	}

	getArg("-in", inputFileName, argc, argv, "!");
	getArg("-matrix", scoringMatrixFileName, argc, argv, "!");
	getArg("-out", outputFileName, argc, argv, "!");

	AlignAlgorithm alignAlgorithm = global;
	if(getArg("-local", argc, argv))
		alignAlgorithm = local;
	if(getArg("-freeshift", argc, argv))
		alignAlgorithm = freeshift;

	getArg("m", matrixFileName, argc, argv, "blosum30.dat");
	getArg("o", openGapPenalty, argc, argv, 10.0);
	getArg("e", extensionGapPenalty, argc, argv, 0.1);

	upgma = getArg("-upgma", argc, argv);
	verbose = getArg("-verbose", argc, argv);
	rooted = getArg("-rooted", argc, argv);

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
	   ERROR("PhylogeneticTrees needs input FASTA file.", exception);

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
	dm->init();

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

    Cluster* tree = builder->performClustering();

    //
    // Newick export
    //

    std::cout<<"\nExporting phylogenetic tree in the Newick format...\n";
    PhylogeneticTreeExporter* exporter = outputFileName == "!" ? new NewickPhylogeneticTreeExporter() : new NewickPhylogeneticTreeExporter(outputFileName);
    exporter->save(tree);

    if(!upgma && rooted){
    	std::cout<<"\nExporting phylogenetic rooted-tree in the Newick format...\n";
    	Cluster* rootedTree = Cluster::midpointRooting(tree, sequences.size());
    	PhylogeneticTreeExporter* expoter = outputFileName == "!" ? new NewickPhylogeneticTreeExporter() : new NewickPhylogeneticTreeExporter("rooted_"+outputFileName);
    	expoter->save(rootedTree);
    }
}
