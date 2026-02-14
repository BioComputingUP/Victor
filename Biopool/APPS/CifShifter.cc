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

#include <string>
#include <GetArg.h>
#include <CifLoader.h>
#include <CifSaver.h>

using namespace Victor;
using namespace Victor::Biopool;

void sShowHelp() {
    cout << "CIF Shifter\n"
            << "Allows to shift all residues in file by fixed offset.\n"
            << " Options: \n"
            << "\t-i <filename> \t\t Input CIF file\n"
            << "\t-o <filename> \t\t Output CIF file\n"
            << "\t[-p <number>] \t\t Positive *residue* offset\n"
            << "\t[-n <number>] \t\t Negative *residue* offset\n"
            << "\t[-P <number>] \t\t Positive *atom* offset\n"
            << "\t[-N <number>] \t\t Negative *atom* offset\n"
            << "\t[--nohydrogen] \t\t Skip hydrogen atoms\n"
            << "\t[--renum] \t\t Reset residue numbering starting from 1\n"
            << "\n";
}

void sAddLine() {
    cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

void sRenumberAtoms(Spacer& sp) {
    unsigned int counter = 1;
    for (unsigned int i = 0; i < sp.sizeAmino(); i++)
        for (unsigned int j = 0; j < sp.getAmino(i).size(); j++) {
            sp.getAmino(i)[j].setNumber(counter);
            counter++;
        }
}

int main(int argc, char** argv) {

    // --------------------------------------------------
    // 0. treat options
    // --------------------------------------------------

    if (getArg("h", argc, argv)) {
        sShowHelp();
        return 1;
    };

    string inputFile, outputFile;
    int offset, offsetAtom, tmp;
    getArg("i", inputFile, argc, argv, "!");
    getArg("o", outputFile, argc, argv, "!");

    getArg("p", offset, argc, argv, 0);
    getArg("n", tmp, argc, argv, 0);
    offset -= tmp;
    
    getArg("P", offsetAtom, argc, argv, 0);
    getArg("N", tmp, argc, argv, 0);
    offsetAtom -= tmp;

    bool noHydrogen = getArg("-nohydrogen", argc, argv);
    bool renumber = getArg("-renum", argc, argv);

    if ((inputFile == "!") || (outputFile == "!")) {
        cout << "Missing file specification. Aborting. (-h for help)" << endl;
        return -1;
    }

    if ((offset == 0) && (!renumber) && (!noHydrogen)) {
        cout << "Warning: Offset is zero. Mistake? \n";
    }
    // --------------------------------------------------
    // 1. read structure
    // --------------------------------------------------

    ifstream inFile(inputFile.c_str());

    if (!inFile)
        ERROR("File does not exist.\n", exception);

    CifLoader cl(inFile);

    if (noHydrogen)
        cl.setNoHAtoms();

    Protein prot;
    prot.load(cl);
    unsigned int zero = 0;
    Spacer sp = *(prot.getSpacer(zero));
    inFile.close();

    // --------------------------------------------------
    // 1.1 renumber residues (if necessary)
    // --------------------------------------------------

    if (renumber) {
        cout << "Renumbering...\n";
        sp.setStartOffset(1);
        sp.removeAllGaps();
        sRenumberAtoms(sp);
    }

    // --------------------------------------------------
    // 2. shift offset
    // --------------------------------------------------

    sp.setStartOffset(offset + sp.getStartOffset());

    if (offsetAtom != 0)
        sp.setAtomStartOffset(offsetAtom + sp.getAtomStartOffset());

    // --------------------------------------------------
    // 3. write model to disk
    // --------------------------------------------------

    ofstream outFile(outputFile.c_str());
    if (!outFile)
        ERROR("File not found.", exception);
    CifSaver cs(outFile);

    sp.save(cs);
    outFile.close();
    
    return 0;
}

