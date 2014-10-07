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
 * @Description This program generates tables of protein entries for the lobo algorithm.
 */
#include <string>
#include <GetArg.h>
#include <LoopTable.h>
#include <IntCoordConverter.h>

using namespace Biopool;

int main(int nArgs, char* argv[]) {
    if (getArg("h", nArgs, argv))
        cout << "This program generates tables of protein entries for the Lobo "
            << "algorithm.\n"
            << "Usage: \n"
            << " LoopTableTest -A<src1> -B<src2> -O<file> -R<rama> -S<size> "
            << "[-v] [-h] \n"
            << "\t -A<file> \t Name of 1st input file or 1.\n"
            << "\t -B<file> \t Name of 2nd input file or 1.\n"
            << "\t -O<file> \t Name of the output file.\n"
            << "\t -R<file> \t Name of the ramachandran data file.\n"
            << "\t -S<size> \t Size of output file: xxs xs s m l xl. \n"
            << "\t --tor <value>  Set torsion angle flexibility factor.\n"
            << "\t --bl <value>  Set bond length flexibility factor.\n"
            << "\t --ba <value>  Set bond angle flexibility factor.\n"
            << "\t -v \t\t Verbose.\n"
            << "\t -h \t\t Display this help.\n"
            << endl;

    bool verbose = getArg("v", nArgs, argv);

    string inputFile1, inputFile2, outputFile;
    string rama, size;
    getArg("A", inputFile1, nArgs, argv, "!");
    getArg("B", inputFile2, nArgs, argv, "!");
    getArg("O", outputFile, nArgs, argv, "!");
    getArg("R", rama, nArgs, argv, "!");
    getArg("S", size, nArgs, argv, "!");

    float bl, ba, tor;
    getArg("-bl", bl, nArgs, argv, -1);
    getArg("-ba", ba, nArgs, argv, -1);
    getArg("-tor", tor, nArgs, argv, -1);
    if (bl > -1) {
        cout << "bl = " << bl << "\n";
        LoopTableEntry::setBondLengthTol(bl);
    }
    if (ba > -1) {
        cout << "ba = " << ba << "\n";
        LoopTableEntry::setBondAngleTol(ba);
    }

    if (tor > -1) {
        cout << "tor = " << tor << "\n";
        RamachandranData::setAngleTol(tor);
    }

    if ((inputFile1 == "!") || (inputFile2 == "!") || (outputFile == "!")
            || (rama == "!") || (size == "!")) {
        cout << "Missing file specification. Aborting. (-h for help)" << endl;
        return -1;
    }

    if (verbose)
        cout << "Verify: \n" << "Input from: " << inputFile1 << "  and  "
            << inputFile2 << "\n"
            << "Output to: " << outputFile << "\t Size: " << size << endl;

    RamachandranData* ramaData = new RamachandranData;
    ifstream ramaIn(rama.c_str());
    if (!ramaIn)
        ERROR("Rama filename invalid.", exception);
    ramaData->load(ramaIn);

    LoopTable P1, P2, P3;
    P1.setRama(ramaData);
    P2.setRama(ramaData);
    P3.setRama(ramaData);

    if (inputFile1 == "1")
        P1.setToSingleAminoAcid();
    else
        P1.read(inputFile1);

    if (inputFile2 == "1")
        P2.setToSingleAminoAcid();
    else
        P2.read(inputFile2);

    if (size == "xxs")
        P3.concatenate(P1, P2, 1, 1);
    else if (size == "xs")
        P3.concatenate(P1, P2, 128, 1);
    else if (size == "s")
        P3.concatenate(P1, P2, 16384, 1);
    else if (size == "m")
        P3.concatenate(P1, P2, 131072, 1);
    else if (size == "l")
        P3.concatenate(P1, P2, 1048576, 1);
    else if (size == "xl")
        P3.concatenate(P1, P2, 2097152, 1);
    else
        ERROR("Error selecting table size.", exception);

    if (verbose)
        cout << "writing...\n";

    P3.write(outputFile);

    if (verbose)
        P3.printTable(1);

    return 0;
}
