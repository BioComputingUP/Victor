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
 */
#include <Protein.h>
#include <PdbLoader.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <GetArg.h>

using namespace Biopool;

void sShowHelp() {
    cout << "Pdb 2 Seq $Revision: 2.0 $ -- converts a PDB file into SEQ\n"
            << "(torsion angles) protein structure backbone torsion angles\n"
            << " Options: \n"
            << "\t-i <filename> \t Input PDB file\n"
            << "\t-o <filename> \t Output to file (default stdout)\n"
            << "\t-c <id>       \t Chain identifier to read\n"
            << "\t--all         \t All chains\n"
            << "\t-m <number>   \t Model number to read (NMR only, default is first model)\n"
            << "\t--chi         \t Write Chi angles (default false)\n"
            << "\t-v            \t verbose output\n\n"
            << "\tIf both -c and --all are missing, only the first chain is processed.\n\n";

}

int main(int argc, char* argv[]) {

    if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    }

    string inputFile, outputFile, chainID;
    unsigned int modelNum;
    bool chi, all;

    getArg("i", inputFile, argc, argv, "!");
    getArg("o", outputFile, argc, argv, "!");
    getArg("c", chainID, argc, argv, "!");
    getArg("m", modelNum, argc, argv, 999);
    all = getArg("-all", argc, argv);
    chi = getArg("-chi", argc, argv);

    // Check input file
    if (inputFile == "!") {
        cout << "Missing input file specification. Aborting. (-h for help)" << endl;
        return -1;
    }
    ifstream inFile(inputFile.c_str());
    if (!inFile)
        ERROR("Input file not found.", exception);


    PdbLoader pl(inFile);

    // Set PdbLoader variables
    pl.setModel(modelNum);
    pl.setNoHAtoms();
    pl.setNoHetAtoms();
    pl.setNoSecondary();
    if (!getArg("v", argc, argv)) {
        pl.setNoVerbose();
    }


    // Check chain args
    if ((chainID != "!") && all) {
        ERROR("You can use --all or -c, not both", error);
    }
    // User selected chain
    if (chainID != "!") {
        if (chainID.size() > 1)
            ERROR("You can choose only 1 chain", error);
        pl.setChain(chainID[0]);
    }// All chains
    else if (all) {
        pl.setAllChains();
    }// First chain
    else {
        pl.setChain(pl.getAllChains()[0]);
    }

    // Load the protein object
    Protein prot;
    prot.load(pl);

    // Open the proper output stream (file or stdout)
    std::ostream* os = &cout;
    std::ofstream fout;
    if (outputFile != "!") {
        fout.open(outputFile.c_str());
        if (!fout) {
            ERROR("Could not open file for writing.", exception);
        } else {
            os = &fout;
        }
    }


    Spacer* sp;
    for (unsigned int i = 0; i < prot.sizeProtein(); i++) {

        sp = prot.getSpacer(i);

        // Write the sequence
        SeqSaver ss(*os);
        if (!chi)
            ss.setWriteChi(false);
        sp->save(ss);
    }

    return 0;
}
