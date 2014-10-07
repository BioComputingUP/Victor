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

// Includes:

#include <AminoAcidHydrogen.h>
#include <IntCoordConverter.h>
#include <IoTools.h>
#include <vector>



// Global constants, typedefs, etc. (to avoid):
using namespace Biopool;

double BOND_LENGTH_H_TO_ALL = 1.00;
map<AminoAcidCode, vector<vector<string> > > AminoAcidHydrogen::paramH;

/**
 *@Description   Load a file with arguments for building hydrogens thank of
 * the IntCoordConverter class
 *@param   void
 *@return  void
 */
void
AminoAcidHydrogen::loadParam(string inputFile) {
    ifstream input(inputFile.c_str());
    if (!input)
        ERROR("File not found.", exception);

    input.clear(); // reset file to previous content 
    input.seekg(0, ios::beg);

    string tok = "";
    vector<string> tokens; // One parsed line splitted by white spaces
    AminoAcidCode aaCode;

    string line;
    line = readLine(input); // header
    //line = readLine(input);
    // read all lines
    do {
        if (line[0] != '#') {
            for (unsigned int i = 0; i < line.size(); i++) {
                if (line[i] == ' ') {
                    tokens.push_back(tok);
                    tok.clear();
                } else {
                    tok += line[i];
                }
            }
            // last tok
            tokens.push_back(tok);
            tok.clear();

            aaCode = aminoAcidThreeLetterTranslator(tokens.front());

            if (paramH.find(aaCode) == paramH.end()) {
                paramH[aaCode] = vector<vector<string > >();
            }
            paramH[aaCode].push_back(tokens);
            tokens.clear();
        }

        line = readLine(input);
    } while (input);

    /*
    for (unsigned int i=0; i<paramH[aaCode].size();i++){
        for (unsigned int j=0; j<paramH[aaCode][i].size();j++){
            cout << i << " " << j << " " << paramH[aaCode][i][j] << "\n";
        }
    }*/

}

/**
 *@Description   Add Hydrogen to an AminoAcid
 *@param   void
 *@return  void
 */
void
AminoAcidHydrogen::setHydrogen(AminoAcid* aa, bool verbose) {



    if (aa->isMember(H) || aa->isMember(HA)) {
        return;
    }

    AminoAcidCode aaCode = AminoAcidCode(aa->getCode());
    string aaName = aminoAcidThreeLetterTranslator(aaCode);

    if (verbose) {
        cout << aaName << "\n";

    }

    vector<vector<string > > paramList = paramH[aaCode];
    vector<string> args;


    //TRP z H N CA 118.50 C 119.50 -1
    // Add H-N backbone, it takes previus aminoacid reference



    if (((*aa).sizeInBonds() > 0) && ((*aa).getType1L() != 'P')) {
        IntCoordConverter icc;
        Atom atH;
        atH.setType("H");

        AminoAcid before = (*aa).getInBond(0);
        atH.bindIn((*aa)[N]);
        (*aa).addAtom(atH);


        icc.zAtomToCartesian((*aa)[N], BOND_LENGTH_H_TO_ALL, (*aa)[CA], 118.50, before[C], 119.50, -1, (*aa)[H]);

    } else {
        if (verbose) {
            cout << "First residue, undefined H rotation!\n";
        }
    }

    AminoAcid aaCopy = (*aa);

    for (unsigned int i = 0; i < paramList.size(); i++) {

        args.clear();
        args = paramList[i];

        IntCoordConverter icc;

        int chiral;

        AtomCode atHCod = AtomTranslator(args[2]);
        AtomCode atBindCod = AtomTranslator(args[3]); // Binding atom
        AtomCode atAngle1Cod = AtomTranslator(args[4]); // Angle partner 1
        AtomCode atAngle2Cod;

        if (args[1] == "z") {
            atAngle2Cod = AtomTranslator(args[6]); // Angle partner 2
        } else if (args[1] == "t") {
            atAngle2Cod = AtomTranslator(args[5]); // Angle partner 2
        } else {
            if (verbose)
                cout << args[1] << " Wrong H function type!\n";
        }

        if (verbose)
            cout << args[2] << "\n";

        if (args[1] == "z") {

            Atom atH;

            chiral = atoi(args[8].c_str());

            atH.setType(args[2]);

            if (aa->getSideChain().isMember(atBindCod)) {

                atH.bindIn(aa->getSideChain()[atBindCod]);
                // WARNING: It changes coords when it shouldn't
                aa->getSideChain().addAtom(atH);
                // PATCH: to prevent coords shift caused by addAtom (ex. CZ in PHE)
                for (unsigned int i = 0; i < aa->getSideChain().size(); i++) {
                    if (!isHAtom(aa->getSideChain().getAtom(i).getCode())) {
                        if (isKnownAtom(aa->getSideChain().getAtom(i).getCode())) {
                            aa->getSideChain().getAtom(i).setCoords(aaCopy.getSideChain()[aa->getSideChain().getAtom(i).getCode()].getCoords());
                        } else {
                            if (verbose)
                                cout << aa->getSideChain().getAtom(i).getCode() << " " << aa->getSideChain().getAtom(i).getNumber() << " " << AtomTranslator(aa->getSideChain().getAtom(i).getCode()) << "\n";
                        }
                    }
                }

                if (aa->getSideChain().isMember(atAngle1Cod)) {
                    if (aa->getSideChain().isMember(atAngle2Cod)) {
                        icc.zAtomToCartesian(aa->getSideChain()[atBindCod], BOND_LENGTH_H_TO_ALL, aa->getSideChain()[atAngle1Cod], atof(args[5].c_str()), aa->getSideChain()[atAngle2Cod], atof(args[7].c_str()), chiral, aa->getSideChain()[atHCod]);
                    } else if ((*aa).isMember(atAngle2Cod)) {
                        icc.zAtomToCartesian(aa->getSideChain()[atBindCod], BOND_LENGTH_H_TO_ALL, aa->getSideChain()[atAngle1Cod], atof(args[5].c_str()), (*aa)[atAngle2Cod], atof(args[7].c_str()), chiral, aa->getSideChain()[atHCod]);
                    } else {
                        aa->getSideChain().removeAtom(atH);
                        if (verbose)
                            cout << " Wrong partner2 atom!\n";
                    }
                } else if ((*aa).isMember(atAngle1Cod)) {
                    if (aa->getSideChain().isMember(atAngle2Cod)) {
                        icc.zAtomToCartesian(aa->getSideChain()[atBindCod], BOND_LENGTH_H_TO_ALL, (*aa)[atAngle1Cod], atof(args[5].c_str()), aa->getSideChain()[atAngle2Cod], atof(args[7].c_str()), chiral, aa->getSideChain()[atHCod]);
                    } else if ((*aa).isMember(atAngle2Cod)) {
                        icc.zAtomToCartesian(aa->getSideChain()[atBindCod], BOND_LENGTH_H_TO_ALL, (*aa)[atAngle1Cod], atof(args[5].c_str()), (*aa)[atAngle2Cod], atof(args[7].c_str()), chiral, aa->getSideChain()[atHCod]);
                    } else {
                        aa->getSideChain().removeAtom(atH);
                        if (verbose)
                            cout << " Wrong partner2 atom!\n";
                    }
                } else {
                    aa->getSideChain().removeAtom(atH);
                    if (verbose)
                        cout << " Wrong partner1 atom!\n";
                }
            } else if ((*aa).isMember(atBindCod)) {

                atH.bindIn((*aa)[atBindCod]);
                (*aa).addAtom(atH);

                if (aa->getSideChain().isMember(atAngle1Cod)) {

                    if (aa->getSideChain().isMember(atAngle2Cod)) {
                        icc.zAtomToCartesian((*aa)[atBindCod], BOND_LENGTH_H_TO_ALL, aa->getSideChain()[atAngle1Cod], atof(args[5].c_str()), aa->getSideChain()[atAngle2Cod], atof(args[7].c_str()), chiral, (*aa)[atHCod]);
                    } else if ((*aa).isMember(atAngle2Cod)) {
                        icc.zAtomToCartesian((*aa)[atBindCod], BOND_LENGTH_H_TO_ALL, aa->getSideChain()[atAngle1Cod], atof(args[5].c_str()), (*aa)[atAngle2Cod], atof(args[7].c_str()), chiral, (*aa)[atHCod]);
                    } else {
                        (*aa).removeAtom(atH);
                        if (verbose)
                            cout << " Wrong partner2 atom!\n";
                    }
                } else if ((*aa).isMember(atAngle1Cod)) {
                    if (aa->getSideChain().isMember(atAngle2Cod)) {
                        icc.zAtomToCartesian((*aa)[atBindCod], BOND_LENGTH_H_TO_ALL, (*aa)[atAngle1Cod], atof(args[5].c_str()), aa->getSideChain()[atAngle2Cod], atof(args[7].c_str()), chiral, (*aa)[atHCod]);
                    } else if ((*aa).isMember(atAngle2Cod)) {
                        icc.zAtomToCartesian((*aa)[atBindCod], BOND_LENGTH_H_TO_ALL, (*aa)[atAngle1Cod], atof(args[5].c_str()), (*aa)[atAngle2Cod], atof(args[7].c_str()), chiral, (*aa)[atHCod]);
                    } else {
                        (*aa).removeAtom(atH);
                        if (verbose)
                            cout << " Wrong partner2 atom!\n";
                    }
                } else {
                    (*aa).removeAtom(atH);
                    if (verbose)
                        cout << " Wrong partner1 atom!\n";
                }
            } else {
                if (verbose)
                    cout << args[3] << " Wrong binding atom!\n";
            }
        }

    }

    if (verbose)
        cout << aaName << " loaded\n";
}
