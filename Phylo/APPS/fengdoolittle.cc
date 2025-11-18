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
#include <AlignmentBase.h>
#include <GetArg.h>
#include <iostream>

using namespace Victor::Phylo;
using namespace Victor;



/** @brief Shows a helper. */
static void show_help() {
    cout << "FENG DOOLITTLE - MULTIPLE SEQUENCE ALIGNMENT TOOL\n"
         << "TThis program calculates a multiple sequence alignment.\n"
         << "Options:\n"
         << "--in <name>    \t Name of input FASTA file\n"
         << "[--out <name>] \t Name of output alignment file (default: to screen)\n"
         << "[-o <double>]  \t Open gap penalty (default: 10.0)\n"
         << "[-e <double>]  \t Extension gap penalty (default: 0.1)\n"
         << "[-v]           \t Verbose (default: no)\n"
         << "\n" << endl;
}



/** @brief Feng-Doolittle multiple sequence alignment program.
 * @param[in]           argc      ARGument Counter
 * @param[in]           argv      ARGument Vector
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 * @return 0 in case of success
 */
int main(int argc, char *argv[]) {
    string input_filename, output_filename;
    double gap_open, gap_extension;
    bool verbose;


    ////////////////////////////////////////////////////////////////////
    // Reads options
    if (getArg("h", argc, argv)) {
        show_help();
        return 0;
    }

    getArg("-in",  input_filename,  argc, argv, "!");
    getArg("-out", output_filename, argc, argv, "!");
    getArg("o",    gap_open,        argc, argv, 10.0);
    getArg("e",    gap_extension,   argc, argv, 0.1);
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


    ////////////////////////////////////////////////////////////////////
    // Creates FengDoolittle object
    if (verbose) cout << "Configuring Feng-Doolittle algorithm..." << endl;
    FengDoolittle feng_doolittle(gap_open, gap_extension);
    if (verbose) cout << "Generating multiple sequence alignment..." << endl;
    MultipleAlignment MSA = feng_doolittle(alignment);


    ////////////////////////////////////////////////////////////////////
    // Saves output
    if (verbose) cout << "Saving multiple sequence alignment..." << endl;
    if (output_filename != "!") {
        ofstream clustalw_file(output_filename);
        MSA.saveClustalW(clustalw_file);
    } else {
        MSA.saveClustalW(cout);
    }
    if (verbose) cout << "done." << endl;


    return 0;
}