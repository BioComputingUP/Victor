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
#include <XyzLoader.h>
#include <IoTools.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Victor; using namespace Victor::Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

/**
 *@Description Loads the information from the input stream into an atom group considering the xyz format. 
 *@param  Group reference(Group&)
 *@return changes are made internally
 */
void XyzLoader::loadGroup(Group& gr) {
    PRINT_NAME;
    if (!input)
        return;

    gr.setType(readLine(input)); // read title line

    int n = -1;
    unsigned int i = 0;
    string typeName;
    vgVector3<double> _coord;
    while (readNumber(input, n)) { // read all entries

        Atom at;
        input >> typeName >> _coord;
        at.setType(typeName);
        gr.addAtom(at);
        unsigned int size = gr.size() - 1;
        while (readOnSameLine(input, i)) // read all connections
            if ((i - 1) <= size)
                gr[size].bindIn(gr[i - 1]);
        gr[size].setCoords(_coord);
    }
}

/**
 *@Description  Loads the side chain from the file, and returns the reference to the corresponding amino acid 
 *@param  SideChain reference(SideChain&), AminoAcid pointer(AminoAcid&)
 *@return void
 */
void XyzLoader::loadSideChain(SideChain& sc, AminoAcid* aaRef) {
    PRINT_NAME;
    if (!input)
        return;

    loadGroup(sc);

    if (aaRef != NULL)
        sc.setBackboneRef(aaRef);
}

/**
 *@Description  Loads an aminoacid in xyz format. 
 *  
 *@param  AminoAcid reference
 *@return void
 */
void XyzLoader::loadAminoAcid(AminoAcid& aa) {
    PRINT_NAME;
    if (!input)
        return;

    loadGroup(aa);

    if (checkForKeyword(input, "sidechain")) {
        loadSideChain(aa.getSideChain(), &aa);
        if (connect)
            if (aa.getSideChain().size())
                for (unsigned int i = 0; i < aa.size(); i++)
                    if (aa[i].getType().c_str()[0] == 'C') { // connect backbone
                        if (!aa[i].isOutBond(aa.getSideChain()[0]))
                            aa[i].bindOut(aa.getSideChain()[0]);
                        break;
                    }
    } else
        DEBUG_MSG("XyzLoader::loadAminoAcid: No sidechain found.");

    aa.adjustLeadingN();
}

/**
 *@Description Loads a spacer in xyz format. 
 *  
 *@param  Spacer reference
 *@return void
 */
void XyzLoader::loadSpacer(Spacer& sp) {
    PRINT_NAME;
    if (!input)
        return;

    sp.setType(readLine(input)); // read title line
    while (checkForKeyword(input, "aminoacid")) {
        AminoAcid* aa = new AminoAcid();
        loadAminoAcid(*aa);
        // connect aa to previous chain segment
        if (connect)
            aa->bindIn((*aa)[N], sp.getAmino(sp.size() - 1),
                sp.getAmino(sp.size() - 1)[C]);
        sp.insertComponent(aa);
    }
    if (sp.sizeAmino())
        sp.getAmino(0).adjustLeadingN();
}

/**
 *@Description Loads a Ligand in xyz format. 
 *  
 *@param  Ligand reference
 *@return void
 */

void XyzLoader::loadLigand(Ligand& l) {
    PRINT_NAME;
    ERROR("Not implemented yet", exception);
}






