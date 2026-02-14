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
#include <string>

#include <IoTools.h>
#include <vector3.h>
#include <AtomCode.h>
#include <String2Number.h>
#include <Ligand.h>
#include <Nucleotide.h>
#include <AminoAcidHydrogen.h>

#include "CifLoader.h"

// Global constants, typedefs, etc. (to avoid):

using namespace Victor;
using namespace Victor::Biopool;
using namespace std;

// CONSTRUCTORS/DESTRUCTOR:

CifLoader::CifLoader(istream& _input, ostream& output, bool _permissive, bool _noHAtoms,
        bool _noHetAtoms, bool _noSecondary, bool _noConnection, bool _noWater,
        bool _verb, bool _allChains, string _NULL, bool _onlyMetal,
        bool _noNucleotideChains ) :
input(_input), output(output), permissive(_permissive), valid(true), noHAtoms(_noHAtoms),
noHetAtoms(_noHetAtoms), noSecondary(_noSecondary), noConnection(_noConnection),
noWater(_noWater), verbose(_verb), allChains(_allChains), chain(' '),
model(999), altAtom('A'), helixCode(_NULL),
//sheetCode(_NULL), helixData(), sheetData(), onlyMetalHetAtoms(_onlyMetal), 
sheetCode(_NULL), onlyMetalHetAtoms(_onlyMetal), noNucleotideChains(_noNucleotideChains) {
    cif = new CifStructure(_input, output);
}

CifLoader::~CifLoader() {
    delete cif;
    PRINT_NAME;
}

// PREDICATES:

/**
 * If user selected a Model, it check validity of this choice,
 * otherwise it select first available chain.
 */
void CifLoader::checkModel() {
    if ((model != 999) && (model > getMaxModels())) {
        ERROR("Please check model number", exception);
    }
}

/**
 * If user selected a chain, it check validity of this choice,
 * otherwise it select first available chain.
 */
void CifLoader::checkAndSetChain() {
    vector<char> chainList = getAllChains();

    if (chain != ' ') {
        bool validChain = false;
        for (unsigned int i = 0; i < chainList.size(); i++)
            if (chain == chainList[i]) {
                validChain = true;
                break;
            }
        if (validChain == false) {
            ERROR("Please check chain id. This is not valid", exception);
        }
    } else {
        chain = chainList[0]; //the first valid chain is default choice
    }
}

/**
 * Reads in the maximum allowed number of NMR models, zero otherwise.
 */
unsigned int CifLoader::getMaxModels() {
    input.clear(); // reset file error flags
    input.seekg(0);

    string atomLine = readLine(input);

    unsigned int max = 0;

    // search column's number of the model field in the atom group
    cif->parseGroup("atom", atomLine);
    int col = cif->getGroupColumnNumber("atom", "model");

    if (col != 0) {
        while (input) {
            if (atomLine.substr(0, 4) == "ATOM") {
                max = stoiDEF(cif->getGroupField("atom", atomLine, col));
            }
            atomLine = readLine(input);
        }
    }
    return max;
}

/**
 * Returns all available chain IDs for a PDB file.
 * @return  all available chain IDs
 */
vector<char> CifLoader::getAllChains() {
    vector<char> res;
    char lastChain = ' ';

    input.clear(); // reset file error flags
    input.seekg(0);

    string atomLine = readLine(input);

    unsigned int modelNum = 0;

    cif->parseGroup("atom", atomLine);
    int modelCol = cif->getGroupColumnNumber("atom", "model");
    int chainCol = cif->getGroupColumnNumber("atom", "chain");

    while (input) {
        if (atomLine.substr(0, 4) == "ATOM") {
            modelNum = stoiDEF(cif->getGroupField("atom", atomLine, modelCol));;
            // only consider first model: others duplicate chain IDs
            if (modelNum > 1) {
                break;
            }
            // check for new chains containing amino acids
            char id = (cif->getGroupField("atom", atomLine, chainCol).c_str())[0];
            if (id != lastChain) {
                lastChain = id;
                res.push_back(id);
            }
        }
        atomLine = readLine(input);
    }
    return res;
}

void CifLoader::setOnlyMetalHetAtoms() {
    if (noHetAtoms) {
        ERROR("can't load metal ions if hetAtoms option is disabled", exception);
    }
    onlyMetalHetAtoms = true;
    noWater = true;
}

void CifLoader::setWater() {
    if (noHetAtoms || onlyMetalHetAtoms) {
        ERROR("can't load water if hetAtoms option is disabled\nor onlyMetalHetAtoms is enabled", exception);
    }
    noWater = false;
}

// HELPERS

/**
 * Private helper function to set bond structure after loading the spacer.
 * @param   Spacer reference
 * @return  bool
 */
bool CifLoader::setBonds(Spacer& sp) {
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
 * Private helper function to determine if atom is backbone or sidechain. 
 * @param   Spacer reference
 * @return  bool
 */
bool CifLoader::inSideChain(const AminoAcid& aa, const Atom& at) {
    if (isBackboneAtom(at.getCode())) {
        return false;
    }
    if ((at.getType() == "H") || (at.getType() == "HN")
            || ((at.getType() == "HA") && (!aa.isMember(HA)))
            || (at.getType() == "1HA") || (at.getType() == "1H")
            || (at.getType() == "2H") || (at.getType() == "3H")) {
        return false; // special case for GLY H (code HA)
    }
    return true; // rest of aminoacid is its sidechain
}

/**
 * Try to assigns the secondary structure from the PDB header. If not present
 * uses Spacer's setStateFromTorsionAngles().
 * @param   Spacer reference
 */
void CifLoader::assignSecondary(Spacer& sp) {
    if (helixData.size() + sheetData.size() == 0) {
        sp.setStateFromTorsionAngles();
        return;
    }

    for (unsigned int i = 0; i < helixData.size(); i++) {
        if (helixCode[i] == chain) {
            for (int j = helixData[i].first; j <= const_cast<int&> (helixData[i].second); j++) {
                // important: keep ifs separated to avoid errors
                if (j < sp.maxPdbNumber()) {
                    if (!sp.isGap(sp.getIndexFromPdbNumber(j))) {
                        sp.getAmino(sp.getIndexFromPdbNumber(j)).setState(HELIX);
		    }
		}
            }
        }
    }

    for (unsigned int i = 0; i < sheetData.size(); i++) {
        if (sheetCode[i] == chain) {
            for (int j = sheetData[i].first; j <= const_cast<int&> (sheetData[i].second); j++) {
                // important: keep ifs separated to avoid errors
                if (j < sp.maxPdbNumber()) {
                    if (!sp.isGap(sp.getIndexFromPdbNumber(j))) {
                        sp.getAmino(sp.getIndexFromPdbNumber(j)).setState(STRAND);
		    }
		}
            }
	}
    }
}

/*
void 
CifLoader::loadSpacer(Spacer& sp){
    Protein prot;
    CifLoader::loadProtein(prot);
    sp = prot.getSpacer(0);    
}
 */

/**
 * Parse a single line of a CIF file.
 * @param atomLine the whole CIF line as it is
 * @param tag the first field (keyword) in a CIF line
 * @param lig pointer to a ligan
 * @param aa pointer to an amino acid
 * @return Residue number read from the CIF line.
 */
int
CifLoader::parseCifline(string atomLine, string tag, Ligand* lig, AminoAcid* aa) {
    // get atom id
    int atNum = stoiDEF(cif->getGroupField("atom", atomLine,
            cif->getGroupColumnNumber("atom", "atom id")));
    // get residue number
    int aaNum = stoiDEF(cif->getGroupField("atom", atomLine,
            cif->getGroupColumnNumber("atom", "residue num")));
    // get code for insertion of residues
    char altAaID = cif->getGroupField("atom", atomLine,
            cif->getGroupColumnNumber("atom", "residue ins")).c_str()[0]; 

    // get x, y, z coordinates
    vgVector3<double> coord;
    coord.x = stodDEF(cif->getGroupField("atom", atomLine,
            cif->getGroupColumnNumber("atom", "x")));
    coord.y = stodDEF(cif->getGroupField("atom", atomLine,
            cif->getGroupColumnNumber("atom", "y")));
    coord.z = stodDEF(cif->getGroupField("atom", atomLine,
            cif->getGroupColumnNumber("atom", "z")));

    // get b-factor
    double bfac = 0.0;
    int colBfac = cif->getGroupColumnNumber("atom", "bfac");
    if (colBfac != -1) {
        string sbfac = cif->getGroupField("atom", atomLine, colBfac);
        if (sbfac != "?" && sbfac != ".") {
            bfac = stodDEF(sbfac);
        }
    }

    // get atom name
    string atType = cif->getGroupField("atom", atomLine,
            cif->getGroupColumnNumber("atom", "atom name"));

    // get residue name
    string aaType = cif->getGroupField("atom", atomLine,
            cif->getGroupColumnNumber("atom", "residue name"));

    // take care of deuterium atoms
    if (atType == "D") {
        cerr << "--> " << atType << "\n";
        atType = "H";
    }
    
    // get asym id
    char asymId = cif->getGroupField("atom", atomLine,
	    cif->getGroupColumnNumber("atom", "atom asym")).c_str()[0];

    // get entity id
    string entityId = cif->getGroupField("atom", atomLine,
	    cif->getGroupColumnNumber("atom", "atom entity"));
    
    // get occupancy
    double occ = stodDEF(cif->getGroupField("atom", atomLine,
	    cif->getGroupColumnNumber("atom", "occupancy")));
    
    // get model
    int model = stoiDEF(cif->getGroupField("atom", atomLine,
	    cif->getGroupColumnNumber("atom", "model")));
    
    // Initialize the Atom object
    Atom* at = new Atom();
    at->setNumber(atNum);
    at->setType(atType);
    at->setCoords(coord);
    at->setBFac(bfac);
    at->setAsymId(asymId);
    at->setEntityId(entityId);
    at->setOccupancy(occ);
    at->setModel(model);

    // Ligand object (includes DNA/RNA in "ATOM" field)
    if (tag == "HETATM" ||
            isKnownNucleotide(nucleotideThreeLetterTranslator(aaType))) {
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
    else if (tag == "ATOM  ") {
        // skip N-terminal ACE groups
        if (aaType != "ACE") {
            // DEBUG: it would be nice to load also alternative atoms
            // skip alternative atoms
            if (altAaID != '?' && altAaID != '.') {
                if (verbose)
                    cout << "Warning: Skipping extraneous amino acid entry "
                        << aaNum << " " << atNum << " " << altAaID << ".\n";
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
                cout << "Warning: Skipping N-terminal ACE group "
                    << aaNum << " " << atNum << ".\n";
        }
    }
    delete at;
    return aaNum;
}

/**
 * Core function for CIF file parsing. 
 * @param prot (Protein&)
 */
void CifLoader::loadProtein(Protein& prot) {
    PRINT_NAME;

    vector<char> chainList = getAllChains();

    if (chainList.size() == 0) {
        if (verbose)
            cout << "Warning: Missing chain ID in the CIF,"
		    "assuming the same chain for the entire file.\n";
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
            }
	    // Load only selected chain
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

            input.clear(); // reset file error flags
            input.seekg(0, ios::beg);

            Spacer* sp = new Spacer();
            LigandSet* ls = new LigandSet();

            string atomLine;
            atomLine = readLine(input);
	    //output << "atomLine: " << atomLine << endl;

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
                // read header entry
		if (atomLine.find(cif->getTag("header")) != string::npos
                        && (name == "")) {
                    sp->setType(cif->getInlineField(atomLine));
                }
		// read helix entry
                else if (atomLine.find(cif->getTag("helix")) != string::npos) {
                    cif->parseGroup("helix", atomLine);
		    
		    while (atomLine != "# " && input) {
			start = stoiDEF(cif->getGroupField("helix", atomLine,
				cif->getGroupColumnNumber("helix", "helix start")));
			end = stoiDEF(cif->getGroupField("helix", atomLine, 
				cif->getGroupColumnNumber("helix", "helix end")));
			helixData.push_back(pair<const int, int>(start, end));

			helixCode += cif->getGroupField("helix", atomLine,
				cif->getGroupColumnNumber("helix", "helix chain"));
			
			char s[256];
			input.getline(s, 256);
			atomLine.assign(s);
			//output << "atomLine: " << atomLine << endl;
		    }
                }
		// read sheet entry
                else if (atomLine.find(cif->getTag("sheet range")) != string::npos) {
                    cif->parseGroup("sheet range", atomLine);
		    
		    while (atomLine != "# " && input) {
			start = stoiDEF(cif->getGroupField("sheet range", atomLine,
				cif->getGroupColumnNumber("sheet range", "sheet start")));
			end = stoiDEF(cif->getGroupField("sheet range", atomLine,
				cif->getGroupColumnNumber("sheet range", "sheet end")));
			sheetData.push_back(pair<const int, int>(start, end));

			sheetCode += cif->getGroupField("sheet range", atomLine,
				cif->getGroupColumnNumber("sheet range", "sheet chain"));

			char s[256];
			input.getline(s, 256);
			atomLine.assign(s);
			//output << "atomLine: " << atomLine << endl;
		    }
                }
		// Parse one line of the "ATOM" and "HETATM" fields
                else if (atomLine.substr(0, 6) == "ATOM  " ||
                        atomLine.substr(0, 6) == "HETATM") {
                    tag = atomLine.substr(0, 6);

                    // Control model number
                    readingModel = stouiDEF(cif->getGroupField("atom", atomLine,
			    cif->getGroupColumnNumber("atom", "model")));
                    if (readingModel > model)
                        break;
                    // Get only the first model if not specified
                    if (model == 999) {
                        model = readingModel;
                    }

                    char chainID = cif->getGroupField("atom", atomLine,
			    cif->getGroupColumnNumber("atom", "chain")).c_str()[0];

                    if (chainList[i] == chainID) {
                        if ((model == 999) || (model == readingModel)) {
                            aaNum = stoiDEF(cif->getGroupField("atom", atomLine,
				    cif->getGroupColumnNumber("atom", "residue num")));

                            // Insert the Ligand object into LigandSet
                            if (aaNum != oldAaNum) {
                                // Print some indexes for the debug
                                /* 
                                cout << aa->getType1L() << " offset:" << sp->getStartOffset() << " gaps:" 
                                     << sp->sizeGaps() << " sizeAmino:" <<  sp->sizeAmino() <<  " maxPdbNum:" 
                                     << sp->maxPdbNumber() << " aaNum:" << aaNum  
                                     << " oldAaNum:" << oldAaNum << " lastAa:" << lastAa << "\n";
                                 */
                                // Skip the first empty AminoAcid
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
                                aa = new AminoAcid();
                                lig = new Ligand();
                            }
                            oldAaNum = parseCifline(atomLine, tag, lig, aa);
                        } // end model check
                    } // end chain check
                }
                atomLine = readLine(input);
		//output << "atomLine: " << atomLine << endl;
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
            if (verbose) {
                cout << "Parsing done\n";
            }
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
                            cout << "Warning: Residue number "
                                    << sp->getPdbNumberFromIndex(j)
                                    << " is incomplete and had to be removed.\n";
                        }
                    }
                }
                if (verbose) {
                    cout << "Removed incomplete residues\n";
		}
                // connect aminoacids
                if (!noConnection) {
                    if (!setBonds(*sp)) { // connect atoms...
                        valid = false;
                        if (verbose) {
                            cout << "Warning: Fail to connect residues in chain: "
                                    << chainList[i] << ".\n";
                        }
                    }
                    if (verbose) {
                        cout << "Connected residues\n";
                    }
                }

                // correct position of leading N atom
                sp->setTrans(sp->getAmino(0)[N].getTrans());
                vgVector3<double> tmp(0.0, 0.0, 0.0);
                sp->getAmino(0)[N].setTrans(tmp);
                sp->getAmino(0).adjustLeadingN();
                if (verbose) {
                    cout << "Fixed leading N atom\n";
		}
                // Add H atoms
                if (!noHAtoms) {
                    for (unsigned int j = 0; j < sp->sizeAmino(); j++) {
                        AminoAcidHydrogen::setHydrogen(&(sp->getAmino(j)), false); // second argument is VERBOSE
                    }
                    if (verbose) {
                        cout << "H assigned\n";
                    }
                    if (!noSecondary) {
                        sp->setDSSP(false); // argument is VERBOSE
                        if (verbose) {
                            cout << "DSSP assigned\n";
                        }
                    }
                }

                // assign secondary structure from torsion angles
                if (!noSecondary) {
                    assignSecondary(*sp);
                    if (verbose) {
                        cout << "Torsional SS assigned\n";
                    }
                }
            } else {
                if (verbose) {
                    cout << "Warning: No residues in chain: " << chainList[i] << ".\n";
                }
            }

            ////////////////////////////////////////////////////////////////////
            // Load data into protein object
            Polymer* pol = new Polymer();
            pol->insertComponent(sp);
            if (verbose) {
                cout << "Loaded AminoAcids: " << sp->size() << "\n";
            }
            if (!(noHetAtoms)) {
		//insertion only if LigandSet is not empty
                if (ls->sizeLigand() > 0) { 
                    pol->insertComponent(ls);
                    if (verbose) {
                        cout << "Loaded Ligands: " << ls->size() << "\n";
                    }
                } else {
                    if (verbose) {
                        cout << "Warning: No ligands in chain: " << chainList[i] << ".\n";
                    }
                }
            }

            prot.addChain(chainList[i]);
            prot.insertComponent(pol);

        } // end loadChain
    } // chains iteration

}