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
#include <PdbSaver.h>
#include <IoTools.h>
#include <vector3.h>
// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

/**
 *@Description Saves a group in PDB format.
 *@param group reference 
 *@return void
 */
void PdbSaver::saveGroup(Group& gr) {
    gr.sync();

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

        output << "ATOM" << setw(7) << gr[i].getNumber() << " " << atName
                << " "
                << gr.getType() << " " << chain << setw(4) << aminoOffset << "    "
                << setw(8) << setprecision(3) << gr[i].getCoords().x
                << setw(8) << setprecision(3) << gr[i].getCoords().y
                << setw(8) << setprecision(3) << gr[i].getCoords().z
                << "  1.00" << setw(6) << setprecision(2) << gr[i].getBFac()
                << "           " << atomOneLetter << "\n";

        atomOffset = gr[i].getNumber() + 1;
    }

    //aminoOffset++;
}

/**
 *@Description Saves a sidechain in PDB format. 
 *@param sideChain reference 
 *@return void
 */
void PdbSaver::saveSideChain(SideChain& sc) {
    saveGroup(sc);
}

/**
 *@Description Saves an aminoacid in PDB format.
 *@param AminoAcid reference 
 *@return void
 */
void PdbSaver::saveAminoAcid(AminoAcid& aa) {
    saveGroup(aa);
}

/**
 *@Description Saves a spacer in PDB format. 
 *@param Spacer reference 
 *@return void
 */
void PdbSaver::saveSpacer(Spacer& sp) {
    PRINT_NAME;

    if (sp.size() > 0) {
        unsigned int oldPrec = output.precision();
        ios::fmtflags oldFlags = output.flags();
        output.setf(ios::fixed, ios::floatfield);

        //method of class Component. It checks how deep is the spacer
        if (sp.getDepth() == 0) {
            if (writeTer) {
                output << "HEADER    " << sp.getType() << "\n"
                        << "REMARK    created using Biopool2000 $Revision: 1.6.2.3 $ \n";
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
            output << "ATOM" << setw(7) << sp.getAmino(index)[OXT].getNumber()
                    << "  OXT "
                    << sp.getAmino(index).getType() << " " << chain << setw(4) << aminoOffset
                    << "    " << setw(8) << setprecision(3)
                    << sp.getAmino(index)[OXT].getCoords().x
                    << setw(8) << setprecision(3)
                    << sp.getAmino(index)[OXT].getCoords().y
                    << setw(8) << setprecision(3)
                    << sp.getAmino(index)[OXT].getCoords().z
                    << "  1.00" << setw(6) << setprecision(2)
                    << sp.getAmino(index)[OXT].getBFac() << "           O\n";
        }

        if ((sp.getDepth() == 0) && (writeTer))
            output << "TER    " << setw(4) << atomOffset + 1 << "      "
                << sp.getAmino(sp.sizeAmino() - 1).getType() << "  "
            << setw(4) << aminoOffset << "\n";

        output.precision(oldPrec);
        output.flags(oldFlags);
        aminoOffset = 0; //necessary if the's more than one spacer
        output << "TER\n";
    }

}

/**
 *@Description Saves a Ligand in PDB format. 
 *@param Ligand reference 
 *@return void
 */
void PdbSaver::saveLigand(Ligand& gr) {
    gr.sync();
    unsigned int oldPrec = output.precision();
    ios::fmtflags oldFlags = output.flags();
    output.setf(ios::fixed, ios::floatfield);

    string aaType = gr.getType();


    // DEBUG: write TER for DNA/RNA ligands

    string tag = "HETATM";
    if (isKnownNucleotide(nucleotideThreeLetterTranslator(aaType))) {
        tag = "ATOM  ";
    }

    for (unsigned int i = 0; i < gr.size(); i++) //print all HETATM of a ligand
    {
        string atType = gr[i].getType();
        aaType = gr.getType();
        string atTypeShort; //last column in a Pdb File
        unsigned int atNum = gr[i].getNumber();
        if (atType != aaType) {
            atTypeShort = atType[0];
            atTypeShort = ' ' + atTypeShort;
            atType = ' ' + atType;
        } else {
            atTypeShort = atType;
            aaType = ' ' + aaType;
        }
        while (atType.size() < 4)
            atType = atType + ' ';
        while (aaType.size() < 3)
            aaType = ' ' + aaType;



        output << tag << setw(5) << atNum << " " << setw(4) << atType << " "
                << setw(3) << aaType << " " << chain << setw(4) << ligandOffset << "    "
                << setw(8) << setprecision(3) << gr[i].getCoords().x
                << setw(8) << setprecision(3) << gr[i].getCoords().y
                << setw(8) << setprecision(3) << gr[i].getCoords().z
                << "  1.00" << setw(6) << setprecision(2) << gr[i].getBFac()
                << "          " << atTypeShort << "\n";
    }
    if (tag == "ATOM  ") {
        output << "TER\n";
    }

    ligandOffset++;
    output.precision(oldPrec);
    output.flags(oldFlags);
}

/**
 *@Description Saves a LigandSet in PDB format. 
 *@param LigandSet reference 
 *@return void
 */
void PdbSaver::saveLigandSet(LigandSet& ls) {
    ligandOffset = ls.getStartOffset(); //set the offset for current LigandSet

    for (unsigned int i = 0; i < ls.sizeLigand(); i++) {
        while ((ls.isGap(ligandOffset))
                && (ligandOffset < ls.maxPdbNumber()))
            ligandOffset++;
        ls[i].save(*this);
    }
}

/**
 *@Description Saves a Protein in PDB format. 
 *@param Protein reference 
 *@return void
 */
void PdbSaver::saveProtein(Protein& prot) {
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
 *@Description Writes the SEQRES entry (PDB format) for a spacer.
 *@param Spacer reference 
 *@return void
 */
void PdbSaver::writeSeqRes(Spacer& sp) {
    for (unsigned int i = 0; i < sp.sizeAmino() / 13; i++) {
        output << "SEQRES " << setw(3) << i << "   " << setw(3)
                << sp.sizeAmino() << "   ";
        for (unsigned int j = 0; j < 13; j++)
            output << sp.getAmino((i * 13) + j).getType() << " ";
        output << "\n";
    }
    if (sp.sizeAmino() % 13 > 0) {
        output << "SEQRES " << setw(3) << sp.sizeAmino() / 13 + 1 << "   "
                << setw(3) << sp.sizeAmino() << "   ";
        for (unsigned int j = 13 * (sp.sizeAmino() / 13); j < sp.sizeAmino(); j++)
            output << sp.getAmino(j).getType() << " ";
        output << "\n";
    }
}

/**
 *@Description Writes the secondary information (PDB format) for a spacer, e.g. HELIX,
 *    SHEET, etc.
 *@param sideChain reference 
 *@return void
 */
void PdbSaver::writeSecondary(Spacer& sp) {

}
