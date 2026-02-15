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
#include <Protein.h>
#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <vector3.h>
#include <CIFLoader.h>
#include <CIFSaver.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>
#include <GetArg.h>

using namespace Victor;using namespace Victor::Biopool;

void sShowHelp() {
    cout << "CIF Secondary $Revision: 1.0 $ -- calculate the secondary structure\n"
            << "(torsion angles) protein structure backbone torsion angles\n"
            << "\tOptions: \n"
            << "\t-i <filename> \t Input CIF file\n"
            << "\t-o <filename> \t Output to file(the chain letter is appended)\n"
            << "\t-c <id>       \t Chain identifier to read(default is all chains)\n"
            << "\t-m <number>   \t Model number to read (NMR only, default is first model)\n"
            << "\t-s <1,2,3>    \t SS calculation(default 3): 1 = CIF fields, 2 = torsion angles, 3 = DSSP\n"
            << "\t              \t   1,2:\t H = helix, E = extended(strand,sheet), . = other.\n"
            << "\t              \t     3:\t H = alpha-helix, E = sheet, B = bridge,\n"
            << "\t                       \t G = 3-10-helix, I = pi-helix, T = n-Turn, S = bend\n"
            << "\t--ext         \t Extended output (default false). Write the type of aminoacid:\n"
            << "\t              \t   N = negative charge, \t P = positive charge\n"
            << "\t              \t   h = hydrophilic, \t + = hydrophobic\n"
            << "\t              \t   , = neutral (hydrophobic)\n"
            << "\t-v            \t verbose output\n\n"
            << "\tWhen parsing multiple chains, one file for each chain is created\n\n";
}

void writeOutput(Spacer* sp, bool ext, int ssType, ostream& os) {

    if (ext) {
        for (unsigned int i = 0; i <= sp->sizeAmino(); i++) {
            if ((i + 1) % 10 == 0)
                os << setw(1) << (((i + 1) % 100) / 10);
            else
                os << " ";
            if ((i + 1) % 60 == 0)
                os << "\n";
        }
        os << "\n";

        for (unsigned int i = 0; i < sp->sizeAmino(); i++) {
            os << sp->getAmino(i).getType1L();
            if ((i + 1) % 60 == 0)
                os << "\n";
        }
        os << "\n";
        for (unsigned int i = 0; i < sp->sizeAmino(); i++) {
            switch (sp->getAmino(i).getCode()) {
                case ASP:
                case GLU:
                    os << "P";
                    break;
                case LYS:
                case ARG:
                    os << "N";
                    break;
                case ASN:
                case GLN:
                case SER:
                case THR:
                case HIS:
                    os << "h";
                    break;
                case VAL:
                case LEU:
                case ILE:
                    os << "+";
                    break;
                default:
                    os << ",";
            };
            if ((i + 1) % 60 == 0)
                os << "\n";
        }
        os << "\n";
    }
    if (ssType != 3) {
        for (unsigned int i = 0; i < sp->sizeAmino(); i++) {
            switch (sp->getAmino(i).getState()) {
                case HELIX:
                    os << "H";
                    break;
                case STRAND:
                    os << "E";
                    break;
                default:
                    os << ".";
            };
            if ((i + 1) % 60 == 0)
                os << "\n";
        }
    } else {
        vector<set<char > > ss = sp->getDSSP();
        for (unsigned int i = 0; i < ss.size(); i++) {
            if (!ss[i].empty()) {
                for (set<char> ::iterator it = ss[i].begin(); it != ss[i].end(); ++it) {
                    os << (*it);
                }
            } else {
                os << ".";
            }
            if ((i + 1) % 60 == 0)
                os << "\n";
        }
    }
    os << "\n";
}

int main(int argc, char* argv[]) {

    if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    }

    string inputFile, outputFile, chainID;
    unsigned int modelNum;
    unsigned int ssType;
    bool extendedOutput, all;


    getArg("i", inputFile, argc, argv, "!");
    getArg("o", outputFile, argc, argv, "!");
    getArg("c", chainID, argc, argv, "!");
    getArg("m", modelNum, argc, argv, 999);
    getArg("s", ssType, argc, argv, 3);
    all = getArg("-all", argc, argv);
    extendedOutput = getArg("-ext", argc, argv);




    // Check input file
    if (inputFile == "!") {
        cout << "Missing input file specification. Aborting. (-h for help)" << endl;
        return -1;
    }
    ifstream inFile(inputFile.c_str());
    if (!inFile)
        ERROR("Input file not found.", exception);


    CIFLoader CIFL(inFile);

    // Set CIFLoader variables
    CIFL.setModel(modelNum);
    CIFL.setNoHetAtoms();

    if (!getArg("v", argc, argv)) {
        CIFL.setNoVerbose();
    }

    // Check chain args
    if ((chainID != "!") && all) {
        ERROR("You can use --all or -c, not both", error);
    }
    // User selected chain
    if (chainID != "!") {
        if (chainID.size() > 1)
            ERROR("You can choose only 1 chain", error);
        CIFL.setChain(chainID[0]);
    }// All chains
    else if (all) {
        CIFL.setAllChains();
    }// First chain
    else {
        CIFL.setChain(CIFL.getAllChains()[0]);
    }

    // Load the protein object
    Protein prot;
    prot.load(CIFL);

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
        writeOutput(sp, extendedOutput, ssType, (*os));
    }

    return 0;
}
