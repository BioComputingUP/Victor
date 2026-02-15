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
#include <Spacer.h>
#include <CIFLoader.h>
#include <CIFSaver.h>
#include <IoTools.h>

using namespace Victor;using namespace Victor::Biopool;

int main(int nArgs, char* argv[]) {
    if (nArgs != 3) {
        cout << "CIF Editor $Revision: 1.0 $ -- allows sequential manipulation of "
                << "protein structure backbone torsion angles" << endl;
        cout << "  Usage: \t\t CIFEditor <input_filename> <output_filename> \n";
        return 1;
    };

    Spacer sp;
    ifstream inFile(argv[1]);
    if (!inFile)
        ERROR("File not found.", exception);

    CIFLoader il(inFile);
    sp.load(il);

    cout << "Editing " << argv[1] << " output goes to " << argv[2] << "\n";

    int aaid = -1;
    do {
        cout << "Aminoacid# or -1: ";
        cin >> aaid;
        if (aaid <= -1) {
            cout << "Bye.\n";
            return 0;
        };

        if (aaid >= (int) sp.sizeAmino()) {
            cout << "\t Invalid aa#!\n";
        } else {
            double newVal = 999;
            cout << "\t " << aaid << "  " << sp.getAmino(aaid).getType() << "\n";
            cout << "Phi=   " << sp.getAmino(aaid).getPhi()
                    << "\t new phi or 999:   ";
            cin >> newVal;
            if (newVal != 999)
                sp.getAmino(aaid).setPhi(newVal);
            cout << "Psi=   " << sp.getAmino(aaid).getPsi()
                    << "\t new psi or 999:   ";
            cin >> newVal;
            if (newVal != 999)
                sp.getAmino(aaid).setPsi(newVal);
            cout << "Omega= " << sp.getAmino(aaid).getOmega()
                    << "\t new omega or 999: ";
            cin >> newVal;
            if (newVal != 999)
                sp.getAmino(aaid).setOmega(newVal);

            ofstream outFile2(argv[2]);

            if (!outFile2)
                ERROR("Couldn't write file.", exception);

            CIFSaver pss2(outFile2);

            sp.save(pss2);
        };
    } while (aaid != -1);

    return 0;
}
