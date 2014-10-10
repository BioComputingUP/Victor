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
#include <RelSaver.h>
#include <IoTools.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Victor; using namespace Victor::Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

/**
 *   Saves a group in relative format. 
 * @param gr
 */
void RelSaver::saveGroup(Group& gr) {
    PRINT_NAME;
    gr.sync();
    output << gr.getType() << "\n";
    pSaveAtomVector(gr.giveAtoms());
}

/**
 *   Saves a sidechain in relative format. 
 * @param sc
 */
void RelSaver::saveSideChain(SideChain& sc) {
    PRINT_NAME;
    sc.sync();
    output << sc.getType() << "\n";
    pSaveAtomVector(sc.giveAtoms());
}


/**
 *     Saves an aminoacid in relative format. 
 * @param aa
 */
void RelSaver::saveAminoAcid(AminoAcid& aa) {
    PRINT_NAME;
    aa.sync();
    output << aa.getType() << "\n";
    pSaveAtomVector(aa.giveAtoms());

    // write relative position of atom following C, if it exists:
    // can be implemented more efficiently..
    if (aa.isMember(C) && !aa.isMember(OXT))
        for (unsigned int i = 0; i < aa[C].sizeOutBonds(); i++)
            if (aa[C].getOutBond(i).getCode() == N) {
                output << "  " << setw(4) << aa[C].getOutBond(i).getNumber()
                        << "    OXT  " << aa[C].getOutBond(i).getTrans()
                        << "      " << setw(3) << aa[C].getNumber() << "\n";
                break;
            }
    output << "  sidechain\n  ";
    saveSideChain(aa.getSideChain());
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        RelSaver::saveSpacer
//
//  Description:
//    Saves a spacer in relative format. 
//
// ----------------------------------------------------------------------------

void RelSaver::saveSpacer(Spacer& sp) {
    PRINT_NAME;
    output << sp.getType() << "\n";
    for (unsigned int i = 0; i < sp.size(); i++) {
        output << "aminoacid\n";
        sp[i].save(*this);
    }
}


// HELPER:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        RelSaver::pSaveAtomVector
//
//  Description:
//    Private helper method to save a vector of atoms.
//    Attention: ID realignment is currently commented out.
//
// ----------------------------------------------------------------------------

void RelSaver::pSaveAtomVector(vector<Atom>& va) {
    unsigned old_prec = output.precision();
    ios::fmtflags old_flags = output.flags();
    output.setf(ios::fixed, ios::floatfield);

    for (unsigned int k = 0; k < va.size(); k++) // write all entries
    {
        string atName = va[k].getType();
        if (!isdigit(atName[0]))
            atName = ' ' + atName;
        while (atName.size() < 4)
            atName += ' ';

        output << "  " << setw(4) << va[k].getNumber() << "   " << atName << "  "
                << va[k].getTrans() << "   ";
        for (unsigned int i = 0; i < va[k].sizeInBonds(); i++)
            output << "   " << setw(3) << va[k].getInBond(i).getNumber();
        for (unsigned int i = 0; i < va[k].sizeOutBonds(); i++)

            output << "   " << setw(3) << va[k].getOutBond(i).getNumber();
        output << "\n";
    }

    output.precision(old_prec);
    output.flags(old_flags);
}
