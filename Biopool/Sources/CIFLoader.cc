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
#include <CIFLoader.h>
#include <CIFFileManager.h>

using namespace std;
using namespace Victor;
using namespace Victor::Biopool;

CIFLoader::CIFLoader(istream& _input, bool _permissive, bool _noHAtoms,
        bool _noHetAtoms, bool _noSecondary, bool _noConnection, bool _noWater,
        bool _verb, bool _allChains, string _NULL, bool _onlyMetal,
        bool _noNucleotideChains)
: input(_input), permissive(_permissive), valid(true),
noHAtoms(_noHAtoms), noHetAtoms(_noHetAtoms), noSecondary(_noSecondary),
noConnection(_noConnection), noWater(_noWater), verbose(_verb),
allChains(_allChains), chain(' '), model(999), altAtom('A'), helixCode(_NULL),
sheetCode(_NULL), onlyMetalHetAtoms(_onlyMetal), noNucleotideChains(_noNucleotideChains) {
    //PrintOut("CIFLoader()");
    cifFM = new CIFFileManager(_input);
    cifFM->Init();
}

unsigned int CIFLoader::getMaxModels() {

    unsigned int max = 1;

    for (int i = 0; i < cifFM->AtomLoop.RowsCount(); i++) {
        if (cifFM->AtomLoop.GetCellValueTrim(AtomGroupPDB, i) == "ATOM") {
            unsigned int temp = stouiDEF(cifFM->AtomLoop.GetCellValueTrim(AtomPDBModelNumber, i));
            if (temp > max)
                max = temp;
        }

    }

    return max;
}

void CIFLoader::checkAndSetChain() {
    vector<char> chainList = getAllChains();
    if (chain != ' ') {
        bool validChain = false;
        for (unsigned int i = 0; i < chainList.size(); i++)
            if (chain == chainList[i]) {
                validChain = true;
                break;
            }
        if (validChain == false) {
            ERROR("Please check CHAIN id. This is not valid", exception);
        }
    } else
        chain = chainList[0]; //the first valid chain is default choice
}

void CIFLoader::checkModel() {
    if ((model != 999)&&(model > getMaxModels())) {
        ERROR("Please check MODEL number", exception);
    }
}

vector<char> CIFLoader::getAllChains() {

    //PrintOut("CIFLoader::getAllChains()");
    vector<char> res;
    char last = ' ';

    for (int i = 0; i < cifFM->AtomLoop.RowsCount(); i++) {
        int modelNumber = stoiDEF(cifFM->AtomLoop.GetCellValueTrim(AtomPDBModelNumber, i));
        if (modelNumber > 1)
            break;

        if (cifFM->AtomLoop.GetCellValueTrim(AtomGroupPDB, i) == "ATOM") {
            char current = (cifFM->AtomLoop.GetCellValueTrim(AtomChainID, i).c_str())[0];

            if (current != last) {
                last = current;
                res.push_back(current);
            }
        }
    }

    return res;

}

bool CIFLoader::setBonds(Spacer& sp) {
    //cout << sp.getAmino(0).getType1L() << "\n";
    sp.getAmino(0).setBondsFromPdbCode(true);
    for (unsigned int i = 1; i < sp.size(); i++) {
        //cout << sp.getAmino(i).getType1L() << "\n";
        if (!sp.getAmino(i).setBondsFromPdbCode(true, &(sp.getAmino(i - 1))))
            return false;
    }
    return true;
}

bool CIFLoader::inSideChain(const AminoAcid& aa, const Atom& at) {
    if (isBackboneAtom(at.getCode()))
        return false;

    if ((at.getType() == "H") || (at.getType() == "HN")
            || ((at.getType() == "HA") && (!aa.isMember(HA)))
            || (at.getType() == "1HA") || (at.getType() == "1H")
            || (at.getType() == "2H") || (at.getType() == "3H"))
        return false; // special case for GLY H (code HA)

    return true; // rest of aminoacid is its sidechain
}

void CIFLoader::assignSecondary(Spacer& sp) {
    if (helixData.size() + sheetData.size() == 0) {
        sp.setStateFromTorsionAngles();
        return;
    }

    for (unsigned int i = 0; i < helixData.size(); i++) {
        if (helixCode[i] == chain) {
            for (int j = helixData[i].first; j <= const_cast<int&> (helixData[i].second); j++) {
                // important: keep ifs separated to avoid errors
                if (j < sp.maxPdbNumber())
                    if (!sp.isGap(sp.getIndexFromPdbNumber(j)))
                        sp.getAmino(sp.getIndexFromPdbNumber(j)).setState(HELIX);
            }
        }
    }

    for (unsigned int i = 0; i < sheetData.size(); i++)
        if (sheetCode[i] == chain)
            for (int j = sheetData[i].first; j <= const_cast<int&> (sheetData[i].second); j++) {
                // important: keep ifs separated to avoid errors
                if (j < sp.maxPdbNumber())
                    if (!sp.isGap(sp.getIndexFromPdbNumber(j)))
                        sp.getAmino(sp.getIndexFromPdbNumber(j)).setState(STRAND);
            }
}

//setOnlyMetalHetAtoms
void CIFLoader::setOnlyMetalHetAtoms() {
    if (noHetAtoms) {
        ERROR("can't load metal ions if hetAtoms option is disabled", exception);
    }

    onlyMetalHetAtoms = true;
    noWater = true;
}

//setWater
void CIFLoader::setWater() {
    if (noHetAtoms || onlyMetalHetAtoms) {
        ERROR("can't load water if hetAtoms option is disabled\nor onlyMetalHetAtoms is enabled", exception);
    }
    noWater = false;
}

void CIFLoader::loadProtein(Protein& prot) {
    PRINT_NAME;

    vector<char> chainList = getAllChains();

    if (chainList.size() == 0) {
        if (verbose)
            cout << "Warning: Missing chain ID in the CIF, assuming the same chain for the entire file.\n";
        chainList.push_back(char(' '));
    }

    bool loadChain = false;

    helixCode = "";
    sheetCode = "";

    string path = "data/AminoAcidHydrogenData.txt";
    const char* inputFile = getenv("VICTOR_ROOT");
    if (inputFile == NULL)
        ERROR("Environment variable VICTOR_ROOT was not found.", exception);

    AminoAcidHydrogen::loadParam(((string) inputFile + path).c_str());
    
    for (unsigned int i = 0; i < chainList.size(); i++) {
        loadChain = false;
        // Load all chains
        if (allChains) {
            loadChain = true;
        } else {
            // Load only first chain
            if (chain == ' ') {
                loadChain = true;
                chain = '#';
            }// Load only selected chain
            else if (chainList[i] == chain) {
                loadChain = true;
                chain = '#';
            }
        }

        if (loadChain) {

            if (verbose) {
                cout << "\nLoading chain: ->" << chainList[i] << "<-\n";
            }
            setChain(chainList[i]);

            Spacer* sp = new Spacer();
            LigandSet* ls = new LigandSet();
            AminoAcid* aa = new AminoAcid();
            Ligand* lig = new Ligand();

            int aaNum = -100000; // infinite negative
            int oldAaNum = -100000;

            sp->setType(cifFM->GetInLineValue("struct_keywords.pdbx_keywords"));
            parseHelixData();
            parseSheetData();
            parseAtomHetatmData(chainList[i], aaNum, oldAaNum, aa, sp, lig, ls);
            
            // Print some indexes for the debug
            /*
            cout << aa->getType1L() << " offset:" << sp->getStartOffset() << " gaps:"
                << sp->sizeGaps() << " sizeAmino:" <<  sp->sizeAmino() <<  " maxPdbNum:"
                << sp->maxPdbNumber() << " aaNum:" << aaNum
                << " oldAaNum:" << oldAaNum << "\n";
            */

            // last residue/ligand
            // AminoAcid
            if ((aa->size() > 0) && (aa->getType1L() != 'X')) {
                if (sp->sizeAmino() == 0) {
                    sp->setStartOffset(oldAaNum - 1);
                } else {
                    // Add gaps
                    for (int i = sp->maxPdbNumber() + 1; i < oldAaNum; i++) {
                        sp->addGap(i);
                    }
                }
                sp->insertComponent(aa);

            }
            // Ligand
            if (lig->size() > 0) {
                if (onlyMetalHetAtoms) {
                    if (lig->isSimpleMetalIon()) { // skip not metal ions
                        ls->insertComponent(lig);
                    }
                } else {
                    ls->insertComponent(lig);
                }
            }
            if (verbose)
                cout << "Parsing done\n";

            ////////////////////////////////////////////////////////////////////
            // Spacer processing
            if (sp->sizeAmino() > 0) {

                // correct ''fuzzy'' (i.e. incomplete) residues
                for (unsigned int j = 0; j < sp->sizeAmino(); j++) {
                    if ((!sp->getAmino(j).isMember(O)) ||
                            (!sp->getAmino(j).isMember(C)) ||
                            (!sp->getAmino(j).isMember(CA)) ||
                            (!sp->getAmino(j).isMember(N))) {

                        // remove residue
                        sp->deleteComponent(&(sp->getAmino(j)));

                        // Add a gap for removed residues
                        sp->addGap(sp->getStartOffset() + j + 1);

                        if (verbose) {
                            cout << "Warning: Residue number " << sp->getPdbNumberFromIndex(j) << " is incomplete and had to be removed.\n";
                        }
                    }
                }
                if (verbose)
                    cout << "Removed incomplete residues\n";

                // connect aminoacids
                if (!noConnection) {

                    if (!setBonds(*sp)) { // connect atoms...
                        valid = false;
                        if (verbose)
                            cout << "Warning: Fail to connect residues in chain: " << chainList[i] << ".\n";
                    }
                    if (verbose)
                        cout << "Connected residues\n";
                }

                // correct position of leading N atom
                sp->setTrans(sp->getAmino(0)[N].getTrans());
                vgVector3<double> tmp(0.0, 0.0, 0.0);
                sp->getAmino(0)[N].setTrans(tmp);
                sp->getAmino(0).adjustLeadingN();
                if (verbose)
                    cout << "Fixed leading N atom\n";

                // Add H atoms
                if (!noHAtoms) {
                    for (unsigned int j = 0; j < sp->sizeAmino(); j++) {
                        AminoAcidHydrogen::setHydrogen(&(sp->getAmino(j)), false); // second argument is VERBOSE
                    }
                    if (verbose)
                        cout << "H assigned\n";

                    if (!noSecondary) {
                        sp->setDSSP(false); // argument is VERBOSE
                        if (verbose)
                            cout << "DSSP assigned\n";
                    }
                }

                // assign secondary structure from torsion angles
                if (!noSecondary) {
                    assignSecondary(*sp);
                    if (verbose)
                        cout << "Torsional SS assigned\n";
                }

            } else {
                if (verbose)
                    cout << "Warning: No residues in chain: " << chainList[i] << ".\n";
            }

            ////////////////////////////////////////////////////////////////////
            // Load data into protein object
            Polymer* pol = new Polymer();
            pol->insertComponent(sp);
            if (verbose)
                cout << "Loaded AminoAcids: " << sp->size() << "\n";

            if (!(noHetAtoms)) {
                if (ls->sizeLigand() > 0) { //insertion only if LigandSet is not empty
                    pol->insertComponent(ls);
                    if (verbose)
                        cout << "Loaded Ligands: " << ls->size() << "\n";
                } else {
                    if (verbose)
                        cout << "Warning: No ligands in chain: " << chainList[i] << ".\n";
                }
            }

            prot.addChain(chainList[i]);
            prot.insertComponent(pol);

        } // end loadChain
    } // chains iteration
}

void CIFLoader::parseHelixData() {
    // read helix entry

    int start, end;

    for (int i = 0; i < cifFM->HelixLoop.RowsCount(); i++) {
        start = stoiDEF(cifFM->HelixLoop.GetCellValueTrim(HelixStart, i));
        end = stoiDEF(cifFM->HelixLoop.GetCellValueTrim(HelixEnd, i));

        helixData.push_back(pair<const int, int>(start, end));
        helixCode += (cifFM->HelixLoop.GetCellValueTrim(HelixCode, i)).c_str()[0];
    }
}

void CIFLoader::parseSheetData() {
    // read sheet entry

    int start, end;

    for (int i = 0; i < cifFM->SheetLoop.RowsCount(); i++) {
        start = stoiDEF(cifFM->SheetLoop.GetCellValueTrim(SheetStart, i));
        end = stoiDEF(cifFM->SheetLoop.GetCellValueTrim(SheetEnd, i));

        sheetData.push_back(pair<const int, int>(start, end));
        sheetCode += (cifFM->SheetLoop.GetCellValueTrim(SheetCode, i)).c_str()[0];
    }
}

void CIFLoader::parseAtomHetatmData(const char currentChain, int &aaNum, int &oldAaNum, AminoAcid* aa, Spacer* sp, Ligand* lig, LigandSet* ls) {
    // Parse one line of the "ATOM" and "HETATM" fields

    unsigned int readingModel = model;

    for (int i = 0; i < cifFM->AtomLoop.RowsCount(); i++) {

        string tag = cifFM->AtomLoop.GetCellValueTrim(AtomGroupPDB, i);
        char chainID = (cifFM->AtomLoop.GetCellValueTrim(AtomChainID, i).c_str())[0];

        // Get only the first model if not specified
        readingModel = stouiDEF(cifFM->AtomLoop.GetCellValueTrim(AtomPDBModelNumber, i));
        
        if (readingModel > model)
        {
            break;
        }
        
        if (model == 999)
        {
            model = readingModel;
        }

        if (currentChain == chainID) {
            
            if ((model == 999) || (model == readingModel)) {
                aaNum = stoiDEF(cifFM->AtomLoop.GetCellValueTrim(AtomSeqID, i));

                // Insert the Ligand object into LigandSet
                if (aaNum != oldAaNum) {
                    // Print some indexes for the debug
                    /*
                    cout << aa->getType1L() << " offset:" << sp->getStartOffset() << " gaps:"
                         << sp->sizeGaps() << " sizeAmino:" <<  sp->sizeAmino() <<  " maxPdbNum:"
                         << sp->maxPdbNumber() << " aaNum:" << aaNum
                         << " oldAaNum:" << oldAaNum << endl;
                    cout << "aa->size(): " << aa->size() << " - aa->getType1L()" << aa->getClassName() << endl;
                    */

                    if ((aa->size() > 0) && (aa->getType1L() != 'X')) { // Skip the first empty AminoAcid
                        if (sp->sizeAmino() == 0) {
                            sp->setStartOffset(oldAaNum - 1);
                        } else {
                            // Add gaps
                            for (int i = sp->maxPdbNumber() + 1; i < oldAaNum; i++) {
                                sp->addGap(i);
                            }
                        }

                        sp->insertComponent(aa);
                    }

                    // Ligand
                    if (lig->size() > 0) {

                        if (onlyMetalHetAtoms) {
                            if (lig->isSimpleMetalIon()) { // skip not metal ions
                                ls->insertComponent(lig);
                            }
                        } else {
                            ls->insertComponent(lig);
                        }
                    }

                    aa = new AminoAcid();
                    lig = new Ligand();
                }

                oldAaNum = parseCIFline(i, tag, lig, aa);

            } // end model check
        } // end chain check
    }
}

int CIFLoader::parseCIFline(int atomLineIndex, string tag, Ligand* lig, AminoAcid * aa) {
    
    int atNum = stoiDEF(cifFM->AtomLoop.GetCellValueTrim(AtomID, atomLineIndex)); // stoi convert from string to int
    int aaNum = stoiDEF(cifFM->AtomLoop.GetCellValueTrim(AtomSeqID, atomLineIndex));
    char altAaID = (cifFM->AtomLoop.GetCellValueTrim(AtomInsCode, atomLineIndex)).c_str()[0]; // "Code for insertion of residues"
    if (altAaID == '?')
    {
        altAaID = ' ';
    }
    
    vgVector3<double> coord;
    coord.x = stodDEF(cifFM->AtomLoop.GetCellValueTrim(AtomX, atomLineIndex));
    coord.y = stodDEF(cifFM->AtomLoop.GetCellValueTrim(AtomY, atomLineIndex));
    coord.z = stodDEF(cifFM->AtomLoop.GetCellValueTrim(AtomZ, atomLineIndex));

    double bfac = stodDEF(cifFM->AtomLoop.GetCellValueTrim(AtomBFact, atomLineIndex));

    string atType = cifFM->AtomLoop.GetCellValueTrim(AtomName, atomLineIndex);
    string aaType = cifFM->AtomLoop.GetCellValueTrim(AtomResidueName, atomLineIndex);

    // take care of deuterium atoms
    if (atType == "D") {
        cerr << "--> " << atType << "\n";
        atType = "H";
    }

    //cout << "atNum:" << atNum << " - aaNum:" << aaNum << " - altAaID:" << altAaID << " - coord.x:" << coord.x << " - coord.y:" << coord.y << " - coord.z:" << coord.z << " - bfac:" << bfac << " - atType:" << atType << " - aaType:" << aaType<< " - TAG:" << tag << endl;
        
    // Initialize the Atom object
    Atom* at = new Atom();
    at->setNumber(atNum);
    at->setType(atType);
    at->setCoords(coord);
    at->setOccupancy(stodDEF(cifFM->AtomLoop.GetCellValueTrim(AtomOccupancy, atomLineIndex)));
    at->setBFac(bfac);

    //Ligand object(includes DNA / RNA in "ATOM" field)
    if ((tag == "HETATM") || isKnownNucleotide(nucleotideThreeLetterTranslator(aaType))) {

        if (noWater) {
            if (!(aaType == "HOH")) {
                lig->addAtom(*at);
                lig->setType(aaType);
            }
        } else {
            lig->addAtom(*at);
            lig->setType(aaType);
        }
    }// AminoAcid
    else if ((tag == "ATOM")) {

        // skip N-terminal ACE groups
        if (aaType != "ACE") {

            // DEBUG: it would be nice to load also alternative atoms
            // skip alternative atoms,
            if (altAaID != ' ') {
                if (verbose)
                    cout << "Warning: Skipping extraneous amino acid entry " << aaNum << " " << atNum << " " << altAaID << ".\n";
            } else {
                aa->setType(aaType);
                aa->getSideChain().setType(aaType);

                if (!noHAtoms || isHeavyAtom(at->getCode())) {
                    if (!inSideChain(*aa, *at))
                    {
                        aa->addAtom(*at);
                    }
                    else {
                        aa->getSideChain().addAtom(*at);
                    }
                }
            }

        } else {
            if (verbose)
                cout << "Warning: Skipping N-terminal ACE group " << aaNum << " " << atNum << ".\n";
        }
    }
    delete at;
    return aaNum;
}
