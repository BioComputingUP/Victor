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
#include <PdbLoader.h>
#include <IoTools.h>
#include <vector3.h>
#include <AtomCode.h>
#include <AminoAcid.h>
#include <String2Number.h>
#include <Ligand.h>
#include <Nucleotide.h>
#include <AminoAcidHydrogen.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Victor; using namespace Victor::Biopool;

// CONSTRUCTORS/DESTRUCTOR:

// PREDICATES:

/**
 *  Reads in the maximum allowed number of NMR models, zero otherwise
 *@param none
 */

unsigned int
PdbLoader::getMaxModels() {
    input.clear(); // reset file to previous content 
    input.seekg(0);

    unsigned int max = 0;

    while (input) {
        string atomLine = readLine(input);
        string tag = atomLine.substr(0, 6);

        if (tag == "MODEL ")
            max++;
    }
    return max;
}

/**
 *    If user selected a chain, it check validity of this choice,
 *    otherwise it select first available chain.
 *@param   void
 *@return  void
 */
void
PdbLoader::checkAndSetChain() {
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

/**
 *    If user selected a Model, it check validity of this choice,
 *    otherwise it select first available chain.
 *@param   void
 *@return  void
 */
void
PdbLoader::checkModel() {
    if ((model != 999)&&(model > getMaxModels())) {
        ERROR("Please check MODEL number", exception);
    }
}

/**
 *    Reads in the maximum allowed number of NMR models, zero otherwise.
 *    This is a faster version because we only need to read the line with
 *    EXPDTA
 *@param   void
 *@return  unsigned int
 */
unsigned int
PdbLoader::getMaxModelsFast() {
    input.clear();
    input.seekg(0);
    while (input) {
        string line = readLine(input);
        if (line.substr(0, 6) == "EXPDTA") {
            istringstream str(line);
            string header, nmr, structures;
            unsigned int models;
            str >> header >> nmr >> models >> structures;
            return models;
        }
    }

    return 0;
}

/**
 *    Returns all available chain IDs for a PDB file.
 *    
 *@param   void
 *@return  vector of chars
 */
vector<char>
PdbLoader::getAllChains() {
    vector<char> res;
    char lastChain = ' ';
    input.clear(); // reset file to previous content 
    input.seekg(0);
    unsigned int max = 0;

    while (input) {
        string atomLine = readLine(input);

        if (atomLine.substr(0, 5) == "MODEL")
            max++;
        if (max > 1) // only consider first model: others duplicate chainIDs
            break;

        // check for new chains containing amino acids
        if ((atomLine.substr(0, 4) == "ATOM")
                && (atomLine.substr(21, 1)[0] != lastChain))
            //&& (aminoAcidThreeLetterTranslator(atomLine.substr(17,3)) != XXX))
        {
            lastChain = atomLine.substr(21, 1)[0];
            res.push_back(lastChain);
        }
    }

    return res;
}

/**
 *    Private helper function to set bond structure after loading the spacer.
 *@param   Spacer reference
 *@return  bool
 */
bool PdbLoader::setBonds(Spacer& sp) {
    //cout << sp.getAmino(0).getType1L() << "\n";
    sp.getAmino(0).setBondsFromPdbCode(true);
    for (unsigned int i = 1; i < sp.size(); i++) {
        //cout << sp.getAmino(i).getType1L() << "\n";
        if (!sp.getAmino(i).setBondsFromPdbCode(true, &(sp.getAmino(i - 1))))
            return false;
    }
    return true;
}

/**
 *    Private helper function to determine if atom is backbone or sidechain. 
 *@param   Spacer reference
 *@return  bool
 */
bool PdbLoader::inSideChain(const AminoAcid& aa, const Atom& at) {
    if (isBackboneAtom(at.getCode()))
        return false;
    if ((at.getType() == "H") || (at.getType() == "HN")
            || ((at.getType() == "HA") && (!aa.isMember(HA)))
            || (at.getType() == "1HA") || (at.getType() == "1H")
            || (at.getType() == "2H") || (at.getType() == "3H"))
        return false; // special case for GLY H (code HA)
    return true; // rest of aminoacid is its sidechain
}


/**
 *    Try to assigns the secondary structure from the PDB header. If not present
 *  uses Spacer's setStateFromTorsionAngles().
 *@param   Spacer reference
 */


void PdbLoader::assignSecondary(Spacer& sp) {
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

void
PdbLoader::setOnlyMetalHetAtoms() {
    if (noHetAtoms) {
        ERROR("can't load metal ions if hetAtoms option is disabled", exception);
    }

    onlyMetalHetAtoms = true;
    noWater = true;
}
//setWater

void
PdbLoader::setWater() {
    if (noHetAtoms || onlyMetalHetAtoms) {
        ERROR("can't load water if hetAtoms option is disabled\nor onlyMetalHetAtoms is enabled", exception);
    }
    noWater = false;
}

/*
void 
PdbLoader::loadSpacer(Spacer& sp){
    Protein prot;
    PdbLoader::loadProtein(prot);
    sp = prot.getSpacer(0);    
}
 */

/**
 *   Core function for PDB file parsing. 
 * @param prot (Protein&)
 */

void
PdbLoader::loadProtein(Protein& prot) {

    PRINT_NAME;

    vector<char> chainList = getAllChains();

    if (chainList.size() == 0) {
        if (verbose)
            cout << "Warning: Missing chain ID in the PDB, assuming the same chain for the entire file.\n";
        chainList.push_back(char(' '));
    }


    unsigned int readingModel = model;
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
            }                // Load only selected chain
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

            input.clear(); // reset file to previous content 
            input.seekg(0, ios::beg);

            Spacer* sp = new Spacer();
            LigandSet* ls = new LigandSet();

            string atomLine;
            atomLine = readLine(input);

            int aaNum = -100000; // infinite negative
            int oldAaNum = -100000;
            //int lastAa = -10000;

            AminoAcid* aa = new AminoAcid();
            Ligand* lig = new Ligand();


            int start, end;

            string name = "";
            string tag = "";

            // read all lines
            do {

                tag = atomLine.substr(0, 6);

                if ((tag == "HEADER") && (name == "")) {
                    name = atomLine;
                    sp->setType(name);
                } else if (tag == "MODEL ") {
                    readingModel = stouiDEF(atomLine.substr(6, 10));
                    if (readingModel > model)
                        break;
                    // Get only the first model if not specified
                    if (model == 999) {
                        model = readingModel;
                    }
                }

                    // read helix entry
                else if (tag == "HELIX ") {
                    start = stoiDEF(atomLine.substr(21, 4));
                    end = stoiDEF(atomLine.substr(33, 4));

                    helixData.push_back(pair<const int, int>(start, end));
                    helixCode += atomLine.substr(19, 1).c_str()[0];
                }                    // read sheet entry
                else if (tag == "SHEET ") {
                    start = stoiDEF(atomLine.substr(22, 4));
                    end = stoiDEF(atomLine.substr(33, 4));

                    sheetData.push_back(pair<const int, int>(start, end));
                    sheetCode += atomLine.substr(21, 1).c_str()[0];
                }

                    // Parse one line of the "ATOM" and "HETATM" fields
                else if ((tag == "ATOM  ") || (tag == "HETATM")) {

                    char chainID = atomLine.substr(21, 1)[0];
                    if (chainList[i] == chainID) {

                        if ((model == 999) || (model == readingModel)) {
                            aaNum = stoiDEF(atomLine.substr(22, 4));

                            // Insert the Ligand object into LigandSet
                            if (aaNum != oldAaNum) {
                                // Print some indexes for the debug
                                /* 
                                cout << aa->getType1L() << " offset:" << sp->getStartOffset() << " gaps:" 
                                     << sp->sizeGaps() << " sizeAmino:" <<  sp->sizeAmino() <<  " maxPdbNum:" 
                                     << sp->maxPdbNumber() << " aaNum:" << aaNum  
                                     << " oldAaNum:" << oldAaNum << " lastAa:" << lastAa << "\n";
                                 */
                                if ((aa->size() > 0) && (aa->getType1L() != 'X')) { // Skip the first empty AminoAcid
                                    if (sp->sizeAmino() == 0) {
                                        sp->setStartOffset(oldAaNum - 1);
                                    } else {
                                        // Add gaps
                                        //for (int i = lastAa+1; i < oldAaNum; i++){
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

                            oldAaNum = parsePDBline(atomLine, tag, lig, aa);

                        } // end model check
                    } // end chain check
                }
                atomLine = readLine(input);

            } while (input);

            
            /*
            // Print some indexes for the debug
            cout << aa->getType1L() << " offset:" << sp->getStartOffset() << " gaps:" 
                << sp->sizeGaps() << " sizeAmino:" <<  sp->sizeAmino() <<  " maxPdbNum:" 
                << sp->maxPdbNumber() << " aaNum:" << aaNum  
                << " oldAaNum:" << oldAaNum << " lastAa:" << lastAa << "\n";
             */
            
            // last residue/ligand
            // AminoAcid
            if ((aa->size() > 0) && (aa->getType1L() != 'X')) {
                if (sp->sizeAmino() == 0) {
                    sp->setStartOffset(oldAaNum - 1);
                } else {
                    // Add gaps
                    //for (int i = lastAa+1; i < oldAaNum; i++){
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

/**
 *   Parse a single line of a PDB file.
 * @param atomLine (string) the whole PDB line as it is
 * @param tag (string) = the first field (keyword) in a PDB line
 * @param lig (Ligand) pointer
 * @param aa (AminoAcid) pointer
 * @return Residue number read from the PDB line (int)
 */
int
PdbLoader::parsePDBline(string atomLine, string tag, Ligand* lig, AminoAcid* aa) {

    int atNum = stoiDEF(atomLine.substr(6, 5)); // stoi convert from string to int
    //char altAtID = atomLine.substr(16,1)[0];    // "Alternate location indicator"               
    int aaNum = stoiDEF(atomLine.substr(22, 4));
    char altAaID = atomLine.substr(26, 1)[0]; // "Code for insertion of residues"
    vgVector3<double> coord;
    coord.x = stodDEF(atomLine.substr(30, 8));
    coord.y = stodDEF(atomLine.substr(38, 8));
    coord.z = stodDEF(atomLine.substr(46, 8));
    double bfac = 0.0;
    if (atomLine.length() >= 66) {
        if (atomLine.substr(60, 6) != "      ") { // empty bfac
            bfac = stodDEF(atomLine.substr(60, 6));
        }
    }
    string atType = "";
    for (int i = 11; i < 17; i++) {
        if (atomLine[i] != ' ')
            atType.append(atomLine.substr(i, 1));
    }
    string aaType = "";
    for (int i = 17; i < 20; i++) {
        if (atomLine[i] != ' ')
            aaType.append(atomLine.substr(i, 1));
    }
    // take care of deuterium atoms
    if (atType == "D") {
        cerr << "--> " << atType << "\n";
        atType = "H";
    }

    // Initialize the Atom object
    Atom* at = new Atom();
    at->setNumber(atNum);
    at->setType(atType);
    at->setCoords(coord);
    //at->setCoords(stod(atomLine.substr(30,8)),stod(atomLine.substr(38,8)),stod(atomLine.substr(46,8)));
    at->setBFac(bfac);

    // Ligand object (includes DNA/RNA in "ATOM" field)
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
    }        // AminoAcid
    else if ((tag == "ATOM  ")) {

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
                        aa->addAtom(*at);
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
