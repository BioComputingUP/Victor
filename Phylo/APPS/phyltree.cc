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

#include <phylo.h>
#include <SubMatrix.h>
#include <AGPFunction.h>
#include <Alignment.h>
#include <AlignmentBase.h>
#include <GetArg.h>
#include <iostream>

using namespace Victor::Phylo;
using namespace Victor;



/** @brief Shows an helper. */
static void show_help() {
    cout << "PHYLOGENETIC TREE GENERATOR\n"
         << "This program calculates a phylogenetic tree for a given alignment.\n"
         << "Options:\n"
         << "--in <name>    \t Path to input FASTA file\n"
         << "[--out <name>] \t Path to output Newick file (default: to screen\n"
         << "[-m <name>]    \t Path to substitution matrix file (default: blosum62.dat)\n"
         << "[-o <double>]  \t Open gap penalty (default: 10.0)\n"
         << "[-e <double>]  \t Extension gap penalty (default: 0.1)\n"
         << "[-d <0|1|2>]   \t Distance matrix builder criterion (default: 0)\n"
         << "               \t -d 0: Distance as 1 - percentage of identity\n"
         << "               \t -d 1: Use Levenshtein distance\n"
         << "               \t -d 2: Use Feng-Doolittle distance\n"
         << "[-c <0|1|2>]   \t Clustering algorithm (default:2)\n"
         << "               \t -c 0: UPGMA\n"
         << "               \t -c 1: Fitch-Margoliash\n"
         << "               \t -c 2: Neighbor Joining\n"
         << "[-v]           \t Verbose (default: no)\n"
         << "\n" << endl;
}



/** @brief Builds a phylogenetic tree.
 * @param[in]           argc      ARGument Counter
 * @param[in]           argv      ARGument Vector
 * @return 0 in case of success
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
int main(int argc, char *argv[]) {
    string input_filename, output_filename, matrix_filename;
    unsigned int distance_criterion, clustering;
    double gap_open, gap_extension;
    bool verbose;


    ////////////////////////////////////////////////////////////////////
    // Reads options
    if (getArg("h", argc, argv)) {
        show_help();
        return 0;
    }

    getArg("-in",  input_filename,     argc, argv, "!");
    getArg("-out", output_filename,    argc, argv, "!");
    getArg("m",    matrix_filename,    argc, argv, "blosum62.dat");
    getArg("o",    gap_open,           argc, argv, 10.0);
    getArg("e",    gap_extension,      argc, argv, 0.1);
    getArg("d",    distance_criterion, argc, argv, 0);
    getArg("c",    clustering,         argc, argv, 1);
    verbose = getArg("v", argc, argv);


    ////////////////////////////////////////////////////////////////////
    // Loads and prepares data
    // Alignment
    if (verbose) cout << "Loading alignment data..." << endl;
    Alignment alignment;
    if (input_filename != "!") {
        ifstream fasta_file(input_filename);
        if (!fasta_file) {
            ERROR("Error opening input FASTA file.", exception);
        }
        
        alignment.loadFasta(fasta_file);
        if (alignment.size() < 1) {
            ERROR("Input FASTA file must contain two sequences.", exception);
        }
    } else {
        ERROR("Missing input FASTA file.", exception);
    }

    // Substitution matrix
    if (verbose) cout << "Loading substitution matrix..." << endl;
    SubMatrix substitution_matrix;
    if (matrix_filename != "!") {
        ifstream matrix_file(matrix_filename);
        if (!matrix_file) {
            ERROR("Error opening substitution matrix file.", exception);
        }
        
        substitution_matrix = SubMatrix(matrix_file);
    } else {
        ERROR("Missing substitution matrix file.", exception);
    }

    // Gap function
    AGPFunction gap_function(gap_open, gap_extension);

    // Distance builder
    DistanceMatrixBuilder *build_matrix;
    switch (distance_criterion) {
        case 0:
            build_matrix = new IdentityPercentage(substitution_matrix, gap_function);
            break;
        case 1:
            build_matrix = new LevenshteinDistance();
            break;
        case 2:
            build_matrix = new FengDoolittleDistance(substitution_matrix, gap_function);
            break;
        default:
            ERROR("Invalid distance criterion", exception);
            exit(-1);
    }

    // Clustering algorithm
    ClusteringAlgorithm *build_tree;
    switch (clustering) {
        case 0:
            build_tree = new UPGMA();
            break;
        case 1:
            build_tree = new NJ();
            break;
        default:
            ERROR("Invalid clustering algorithm", exception);
            exit(-1);
    }


    ////////////////////////////////////////////////////////////////////
    // Computes phylogenetic tree
    if (verbose) cout << "Generating distance matrix..." << endl;
    DistanceMatrix distance_matrix = (*build_matrix)(alignment);
    if (verbose) cout << "Generating phylogenetic tree..." << endl;
    RootedTree tree = (*build_tree)(distance_matrix).asRootedTree();


    ////////////////////////////////////////////////////////////////////
    // Saves output
    if (verbose) cout << "Saving tree..." << endl;
    if (output_filename != "!") {
        ofstream newick_file(output_filename);
        tree.saveNewick(newick_file);
    } else {
        tree.saveNewick(cout);
    }
    if (verbose) cout << "done." << endl;

    return 0;
}