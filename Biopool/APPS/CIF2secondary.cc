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

#include <string>
#include <GetArg.h>
#include <Spacer.h>
#include <CIFLoader.h>
#include <IoTools.h>

using namespace Victor;using namespace Victor::Biopool;

void sShowHelp() {
    cout << "CIF 2 Secondary Structure converter\n"
            << "\t H = helix, \t E = extended (strand, sheet), \t . = other.\n"
            << "   Options: \n"
            << "\t-i <filename> \t\t Input file for CIF structure\n"
            << "\n";
}

int main(int nArgs, char* argv[]) {
    if (getArg("h", nArgs, argv)) {
        sShowHelp();
        return 1;
    };
    vector<char> allCh;
    string chainID = "!";
    string inputFile;
    getArg("i", inputFile, nArgs, argv, "!");

    if (inputFile == "!") {
        cout << "Missing file specification. Aborting. (-h for help)" << endl;
        return -1;
    }

    Spacer *sp;
    ifstream inFile(inputFile.c_str());
    if (!inFile)
        ERROR("File not found.", exception);

    CIFLoader il(inFile);
    il.setNoHAtoms();
    // sp.load(il); 
    allCh = il.getAllChains();
    for (unsigned int i = 0; i < allCh.size(); i++)
        cout << "\t," << allCh[i] << ",";
    cout << "\n";

    /*check on validity of chain: 
    if user select a chain then check validity
     else select first valid one by default*/
    if (chainID != "!") {
        bool validChain = false;
        for (unsigned int i = 0; i < allCh.size(); i++) {
            if (allCh[i] == chainID[0]) {
                il.setChain(chainID[0]);
                cout << "Loading chain " << chainID << "\n";
                validChain = true;
                break;
            }
        }
        if (!validChain) {
            cout << "Chain " << chainID << " is not available\n";
            return -1;
        }

    } else {
        chainID[0] = allCh[0];
        cout << "Using chain " << chainID << "\n";
    }
    Protein prot;
    prot.load(il);
    sp = prot.getSpacer(chainID[0]);

    allCh = il.getAllChains();
    cout << ">" << inputFile << "\n";
    for (unsigned int i = 0; i < sp->sizeAmino(); i++) {
        switch (sp->getAmino(i).getState()) {
            case HELIX:
                cout << "H";
                break;
            case STRAND:
                cout << "E";
                break;
            default:
                cout << ".";
        };
        if ((i + 1) % 60 == 0)
            cout << "\n";
    }
    cout << "\n";

    return 0;
}
