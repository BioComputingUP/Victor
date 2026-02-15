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

/* 
 * File:   CIFSaver.cc
 * Author: Francesco Menniti
 * 
 * Created on 22 dicembre 2015, 15.28
 */
#include "CIFSaver.h"

using namespace Victor;
using namespace Victor::Biopool;
using namespace std;

void CIFSaver::saveGroup(Group& gr) {
    gr.sync();

    if (!atomSiteWrited) {
        output << "loop_" << endl;
        output << "_atom_site.group_PDB " << endl;
        output << "_atom_site.id " << endl;
        output << "_atom_site.type_symbol " << endl;
        output << "_atom_site.label_atom_id " << endl;
        output << "_atom_site.label_alt_id " << endl;
        output << "_atom_site.label_comp_id " << endl;
        output << "_atom_site.label_asym_id " << endl;
        output << "_atom_site.label_entity_id " << endl;
        output << "_atom_site.label_seq_id " << endl;
        output << "_atom_site.pdbx_PDB_ins_code " << endl;
        output << "_atom_site.Cartn_x " << endl;
        output << "_atom_site.Cartn_y " << endl;
        output << "_atom_site.Cartn_z " << endl;
        output << "_atom_site.occupancy " << endl;
        output << "_atom_site.B_iso_or_equiv " << endl;
        output << "_atom_site.Cartn_x_esd " << endl;
        output << "_atom_site.Cartn_y_esd " << endl;
        output << "_atom_site.Cartn_z_esd " << endl;
        output << "_atom_site.occupancy_esd " << endl;
        output << "_atom_site.B_iso_or_equiv_esd " << endl;
        output << "_atom_site.pdbx_formal_charge " << endl;
        output << "_atom_site.auth_seq_id " << endl;
        output << "_atom_site.auth_comp_id " << endl;
        output << "_atom_site.auth_asym_id " << endl;
        output << "_atom_site.auth_atom_id " << endl;
        output << "_atom_site.pdbx_PDB_model_num " << endl;

        atomSiteWrited = true;

    }

    for (unsigned int i = 0; i < gr.size(); i++) {
        string atName = gr[i].getType();

        if (atName == "OXT") // cosmetics: OXT has to be output after 
            continue; // the sidechain and therefore goes in saveSpacer

        // Added variable for correcting atom type H (last column in PDBs)
        char atomOneLetter;
        if (!isdigit(atName[0])) {
            atomOneLetter = atName[0];
        } else {
            atomOneLetter = atName[1];
        }

        // Added control for size by Damiano Piovesan
        // example HG12
        if (!isdigit(atName[0]) && (atName.size() < 4))
            atName = ' ' + atName;
        while (atName.size() < 4)
            atName += ' ';

        output << setw(7) << left << "ATOM" <<
                setw(6) << gr[i].getNumber() <<
                setw(2) << atomOneLetter <<
                setw(5) << left << atName <<
                setw(2) << "." <<
                setw(4) << gr.getType() <<
                setw(2) << "?" <<
                setw(2) << "?" <<
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
                setw(2) << "1" <<
                endl;

        atomOffset = gr[i].getNumber() + 1;
    }

    //aminoOffset++;
}

void CIFSaver::saveSideChain(SideChain& sc) {
    saveGroup(sc);
}

void CIFSaver::saveAminoAcid(AminoAcid& aa) {
    saveGroup(aa);
}

void CIFSaver::saveSpacer(Spacer& sp) {
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
		output << "Header" << "   " << sp.getType()
			<< " " << endl;
		output << "# " << endl;
            }
            aminoOffset = 0;
            atomOffset = sp.getAtomStartOffset();
        }

        if (writeSeq)
        {
            writeSeqRes(sp);
        }
        if (writeSecStr)
        {
            writeSecondary(sp);
        }

        aminoOffset = sp.getStartOffset();
        atomOffset = sp.getAtomStartOffset();

        //saving is one ammino at a time
        for (unsigned int i = 0; i < sp.sizeAmino(); i++) {
            aminoOffset++;
            while ((sp.isGap(aminoOffset)) && (aminoOffset < sp.maxPdbNumber())) {
                aminoOffset++;
            }
            sp.getAmino(i).save(*this);
        }

        // cosmetics: write OXT after last side chain
        if (sp.getAmino(sp.sizeAmino() - 1).isMember(OXT)) {
            unsigned int index = sp.sizeAmino() - 1;
            
                output << setw(7) << left << "ATOM" <<
                setw(6) << sp.getAmino(index)[OXT].getNumber() <<
                setw(2) << "O" <<
                setw(5) << left << "OXT" <<
                setw(2) << "." <<
                setw(4) << sp.getAmino(index)[OXT].getType() <<
                setw(2) << "?" <<
                setw(2) << "?" <<
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
                setw(4) << sp.getAmino(index)[OXT].getType() <<
                setw(2) << chain <<
                setw(5) << left << "OXT" <<
                setw(2) << "1" <<
                endl;
        }

        output.precision(oldPrec);
        output.flags(oldFlags);
        aminoOffset = 0; //necessary if the's more than one spacer
    }

}

void CIFSaver::saveLigand(Ligand& gr) {
    gr.sync();
    unsigned int oldPrec = output.precision();
    ios::fmtflags oldFlags = output.flags();
    output.setf(ios::fixed, ios::floatfield);

    string aaType = gr.getType();

    string tag = "HETATM";
    if (isKnownNucleotide(nucleotideThreeLetterTranslator(aaType))) {
	tag = "ATOM";
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

	output << setw(7) << left << tag <<
		setw(6) << gr[i].getNumber() <<
		setw(2) << atTypeShort <<
		setw(5) << left << atType <<
		setw(2) << "." <<
		setw(4) << aaType <<
		setw(2) << '?' <<
		setw(2) << '?' <<
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
		setw(2) << "1" <<
		endl;
    }
    output << "# " << endl;

    ligandOffset++;
    output.precision(oldPrec);
    output.flags(oldFlags);
}

void CIFSaver::saveLigandSet(LigandSet& ls) {
        
    ligandOffset = ls.getStartOffset(); //set the offset for current LigandSet

    for (unsigned int i = 0; i < ls.sizeLigand(); i++) {
        while ((ls.isGap(ligandOffset))
                && (ligandOffset < ls.maxPdbNumber()))
            ligandOffset++;
        ls[i].save(*this);
    }
}

void CIFSaver::saveProtein(Protein& prot) {
    
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

void CIFSaver::writeSeqRes(Spacer& sp) {

    output << "loop_" << endl;
    output << "_entity_poly_seq.entity_id " << endl;
    output << "_entity_poly_seq.num " << endl;
    output << "_entity_poly_seq.mon_id " << endl;
    output << "_entity_poly_seq.hetero " << endl;

    for (unsigned int i = 0; i < sp.sizeAmino(); i++) {

        //string entityId = sp.getAmino(i).getAtom(0).getEntityId();

        output << setw(2) << left << "?" <<
                setw(4) << i + 1 <<
                setw(4) << sp.getAmino(i).getType() <<
                setw(2) << "n" <<
                endl;
    }

    output << "# " << endl;

}

void CIFSaver::writeSecondary(Spacer& sp) {

}


