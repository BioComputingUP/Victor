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
#include <IoTools.h>
#include <vector3.h>

#include "CifSaver.h"

// Global constants, typedefs, etc. (to avoid):
using namespace Victor;
using namespace Victor::Biopool;
using namespace std;

// CONSTRUCTORS/DESTRUCTOR:

CifSaver::CifSaver(ostream& _output) :
output(_output), writeSeq(true), writeSecStr(true), writeTer(true),
atomOffset(0), aminoOffset(0), ligandOffset(0), chain(' '),
atomGroupPrinted(false) {
    cif = new CifStructure(_output);
}

CifSaver::~CifSaver() {
    delete cif;
    PRINT_NAME;
}

// PREDICATES:

// MODIFIERS:

/**
 * Saves a group in CIF format.
 * @param group reference 
 * @return void
 */
void CifSaver::saveGroup(Group& gr) {
    gr.sync();

    if (!atomGroupPrinted) {
	cif->printGroup("atom");
	atomGroupPrinted = true;
    }

    for (unsigned int i = 0; i < gr.size(); i++) {
	string atName = gr[i].getType();

	// cosmetics: OXT has to be output after
	// the sidechain and therefore goes in saveSpacer
	if (atName == "OXT") {
	    continue;
	}

	// Added variable for correcting atom type H (last column in PDBs)
	char atomOneLetter;
	if (!isdigit(atName[0])) {
	    atomOneLetter = atName[0];
	} else {
	    atomOneLetter = atName[1];
	}

	// Added control for size by Damiano Piovesan
	// example HG12
	if (!isdigit(atName[0]) && (atName.size() < 4)) {
	    atName = ' ' + atName;
	}
	while (atName.size() < 4) {
	    atName += ' ';
	}

	// if fields have default values (0 or X), assigns the CIF unknown value (?)
	// or a possibly correct value
	char asymId = gr[i].getAsymId();
	string entityId = gr[i].getEntityId();
	int model = gr[i].getModel();

	if (asymId == 'X') {
	    asymId = '?';
	}

	if (entityId == "0") {
	    entityId = "?";
	}

	if (model == 0) {
	    model = 1;
	}

	output << setw(7) << left << "ATOM" <<
		setw(6) << gr[i].getNumber() <<
		setw(2) << atomOneLetter <<
		setw(5) << left << atName <<
		setw(2) << "." <<
		setw(4) << gr.getType() <<
		setw(2) << asymId <<
		setw(2) << entityId <<
		setw(4) << aminoOffset <<
		setw(2) << "?" <<
		setw(8) << setprecision(3) << gr[i].getCoords().x <<
		setw(8) << setprecision(3) << gr[i].getCoords().y <<
		setw(8) << setprecision(3) << gr[i].getCoords().z <<
		setw(6) << setprecision(2) << gr[i].getOccupancy() <<
		setw(7) << left << setprecision(2) << gr[i].getBFac() <<
		setw(2) << "?" <<
		setw(2) << "?" <<
		setw(2) << "?" <<
		setw(2) << "?" <<
		setw(2) << "?" <<
		setw(2) << "?" <<
		setw(4) << aminoOffset <<
		setw(4) << gr.getType() <<
		setw(2) << chain <<
		setw(5) << left << atName <<
		setw(2) << model <<
		endl;

	atomOffset = gr[i].getNumber() + 1;
    }

    //aminoOffset++;
}

/**
 * Saves a sidechain in CIF format. 
 * @param sideChain side chain to save
 */
void CifSaver::saveSideChain(SideChain& sc) {
    saveGroup(sc);
}

/**
 * Saves an aminoacid in CIF format.
 * @param AminoAcid aminoacid to save
 */
void CifSaver::saveAminoAcid(AminoAcid& aa) {
    saveGroup(aa);
}

/**
 * Saves a spacer in CIF format. 
 * @param Spacer spacer to save
 */
void CifSaver::saveSpacer(Spacer& sp) {
    PRINT_NAME;

    if (sp.size() > 0) {
	unsigned int oldPrec = output.precision();
	ios::fmtflags oldFlags = output.flags();
	output.setf(ios::fixed, ios::floatfield);

	//method of class Component. It checks how deep is the spacer
	if (sp.getDepth() == 0) {
	    if (writeTer) {
		output << "data_" << sp.getType() << endl;
		output << "# " << endl;
		output << cif->getTag("header") << "   " << sp.getType()
			<< " " << endl;
		output << "# " << endl;
	    }
	    aminoOffset = 0;
	    atomOffset = sp.getAtomStartOffset();
	}

	if (writeSeq)
	    writeSeqRes(sp);
	if (writeSecStr)
	    writeSecondary(sp);

	aminoOffset = sp.getStartOffset();
	atomOffset = sp.getAtomStartOffset();

	//saving is one ammino at a time
	for (unsigned int i = 0; i < sp.sizeAmino(); i++) {
	    aminoOffset++;
	    while ((sp.isGap(aminoOffset)) && (aminoOffset < sp.maxPdbNumber())) {
		aminoOffset++;
	    }
	    //cout << i << " " << aminoOffset << "\n";
	    sp.getAmino(i).save(*this);
	}

	// cosmetics: write OXT after last side chain
	if (sp.getAmino(sp.sizeAmino() - 1).isMember(OXT)) {
	    unsigned int index = sp.sizeAmino() - 1;

	    // if fields have default values (0 or X), assigns the CIF unknown value (?)
	    // or a possibly correct value
	    char asymId = sp.getAmino(index)[OXT].getAsymId();
	    string entityId = sp.getAmino(index)[OXT].getEntityId();
	    int model = sp.getAmino(index)[OXT].getModel();

	    if (asymId == 'X') {
		asymId = '?';
	    }

	    if (entityId == "0") {
		entityId = "?";
	    }

	    if (model == 0) {
		model = 1;
	    }

	    output << setw(7) << left << "ATOM" <<
		    setw(6) << sp.getAmino(index)[OXT].getNumber() <<
		    setw(2) << "O" <<
		    setw(5) << left << "OXT" <<
		    setw(2) << "." <<
		    setw(4) << sp.getAmino(index).getType() <<
		    setw(2) << asymId <<
		    setw(2) << entityId <<
		    setw(4) << aminoOffset <<
		    setw(2) << "?" <<
		    setw(8) << setprecision(3) << sp.getAmino(index)[OXT].getCoords().x <<
		    setw(8) << setprecision(3) << sp.getAmino(index)[OXT].getCoords().y <<
		    setw(8) << setprecision(3) << sp.getAmino(index)[OXT].getCoords().z <<
		    setw(6) << setprecision(2) << sp.getAmino(index)[OXT].getOccupancy() <<
		    setw(7) << left << setprecision(2) << sp.getAmino(index)[OXT].getBFac() <<
		    setw(2) << "?" <<
		    setw(2) << "?" <<
		    setw(2) << "?" <<
		    setw(2) << "?" <<
		    setw(2) << "?" <<
		    setw(2) << "?" <<
		    setw(4) << aminoOffset <<
		    setw(4) << sp.getAmino(index).getType() <<
		    setw(2) << chain <<
		    setw(5) << "OXT" <<
		    setw(2) << model <<
		    endl;
	}

	output.precision(oldPrec);
	output.flags(oldFlags);
	aminoOffset = 0; //necessary if the's more than one spacer
    }

}

/**
 * Saves a Ligand in CIF format. 
 * @param Ligand ligand to save
 */
void CifSaver::saveLigand(Ligand& gr) {
    gr.sync();
    unsigned int oldPrec = output.precision();
    ios::fmtflags oldFlags = output.flags();
    output.setf(ios::fixed, ios::floatfield);

    string aaType = gr.getType();

    string tag = "HETATM";
    if (isKnownNucleotide(nucleotideThreeLetterTranslator(aaType))) {
	tag = "ATOM  ";
    }

    //print all HETATM of a ligand
    for (unsigned int i = 0; i < gr.size(); i++) {
	string atType = gr[i].getType();
	aaType = gr.getType();
	string atTypeShort; //last column in a Pdb File

	if (atType != aaType) {
	    atTypeShort = atType[0];
	} else {
	    atTypeShort = atType;
	}

	// if fields have default values (0 or X), assigns the CIF unknown value (?)
	// or a possibly correct value
	char asymId = gr[i].getAsymId();
	string entityId = gr[i].getEntityId();
	int model = gr[i].getModel();

	if (asymId == 'X') {
	    asymId = '?';
	}

	if (entityId == "0") {
	    entityId = "?";
	}

	if (model == 0) {
	    model = 1;
	}

	output << setw(7) << left << tag <<
		setw(6) << gr[i].getNumber() <<
		setw(2) << atTypeShort <<
		setw(5) << left << atType <<
		setw(2) << "." <<
		setw(4) << aaType <<
		setw(2) << asymId <<
		setw(2) << entityId <<
		setw(4) << ligandOffset <<
		setw(2) << "?" <<
		setw(8) << setprecision(3) << gr[i].getCoords().x <<
		setw(8) << setprecision(3) << gr[i].getCoords().y <<
		setw(8) << setprecision(3) << gr[i].getCoords().z <<
		setw(6) << setprecision(2) << gr[i].getOccupancy() <<
		setw(7) << left << setprecision(2) << gr[i].getBFac() <<
		setw(2) << "?" <<
		setw(2) << "?" <<
		setw(2) << "?" <<
		setw(2) << "?" <<
		setw(2) << "?" <<
		setw(2) << "?" <<
		setw(4) << aminoOffset <<
		setw(4) << aaType <<
		setw(2) << chain <<
		setw(5) << atType <<
		setw(2) << model <<
		endl;
    }
    output << "# " << endl;

    ligandOffset++;
    output.precision(oldPrec);
    output.flags(oldFlags);
}

/**
 * Saves a LigandSet in CIF format. 
 * @param LigandSet set of ligands to save
 */
void CifSaver::saveLigandSet(LigandSet& ls) {
    ligandOffset = ls.getStartOffset(); //set the offset for current LigandSet

    for (unsigned int i = 0; i < ls.sizeLigand(); i++) {
	while ((ls.isGap(ligandOffset))
		&& (ligandOffset < ls.maxPdbNumber()))
	    ligandOffset++;
	ls[i].save(*this);
    }
}

/**
 * Saves a Protein in PDB format. 
 * @param Protein protein to save
 */
void CifSaver::saveProtein(Protein& prot) {
    //if (prot.sizeProtein()==0)
    //        ERROR("Empty Protein",exception);

    Spacer* sp = NULL;
    LigandSet* ls = NULL;

    for (unsigned int i = 0; i < prot.sizeProtein(); i++) {
	setChain(prot.getChainLetter(i)); //set the actual chain's ID
	sp = prot.getSpacer(i);
	saveSpacer(*sp);
    }

    for (unsigned int i = 0; i < prot.sizeProtein(); i++) {
	setChain(prot.getChainLetter(i)); //set the actual chain's ID
	ls = prot.getLigandSet(i);

	if (ls != NULL) {
	    saveLigandSet(*ls);
	}
    }
}

/**
 * Writes the SEQRES entry (CIF format) for a spacer.
 * @param Spacer spacer to write
 */
void CifSaver::writeSeqRes(Spacer& sp) {
    cif->printGroup("entity poly");


    for (unsigned int i = 0; i < sp.sizeAmino(); i++) {
	// if fields have default values (0 or X), assigns the CIF unknown value (?)
	// or a possibly correct value
	string entityId = sp.getAmino(i).getAtom(0).getEntityId();

	if (entityId == "0") {
	    entityId = "?";
	}

	output << setw(2) << left << entityId <<
		setw(4) << i + 1 <<
		setw(4) << sp.getAmino(i).getType() <<
		setw(2) << "n" <<
		endl;
    }
    output << "# " << endl;
}

/**
 *  Writes the secondary information (PDB format) for a spacer, e.g. HELIX,
 *  SHEET, etc.
 * @param sideChain reference 
 * @return void
 */
void CifSaver::writeSecondary(Spacer& sp) {

}