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
#include <SeqSaver.h>
#include <IoTools.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Victor; using namespace Victor::Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqSaver::saveSideChain

//
//  Description:
//    Saves a sidechain in SEQ format. 
//
// ----------------------------------------------------------------------------

void SeqSaver::saveSideChain(SideChain& sc, bool header) {
    PRINT_NAME;
    unsigned old_prec = output.precision();
    ios::fmtflags old_flags = output.flags();
    output.setf(ios::fixed, ios::floatfield);
    if (header)
        output << sc.getType() << "   ";

    if (writeChi)
        for (unsigned int i = 0; i < sc.getMaxChi(); i++) // write torsion angles
            output << "   " << setw(8) << setprecision(3) << sc.getChi(i);
    output << "\n";
    output.precision(old_prec);
    output.flags(old_flags);
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqSaver::saveAminoAcid
//
//  Description:
//    Saves an aminoacid in SEQ format. 
//
// ----------------------------------------------------------------------------

void SeqSaver::saveAminoAcid(AminoAcid& aa) {
    PRINT_NAME;
    unsigned old_prec = output.precision();
    ios::fmtflags old_flags = output.flags();
    output.setf(ios::fixed, ios::floatfield);
    output << aa.getType() << "   " // write torsion angles
            << setw(8) << setprecision(3) << aa.getPhi()
            << "   " << setw(8) << setprecision(3) << aa.getPsi()
            << "   " << setw(8) << setprecision(3) << aa.getOmega();
    saveSideChain(aa.getSideChain(), false);
    output.precision(old_prec);
    output.flags(old_flags);
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqSaver::saveSpacer
//
//  Description:
//    Saves a spacer in SEQ format. 
//
// ----------------------------------------------------------------------------

void SeqSaver::saveSpacer(Spacer& sp) {
    PRINT_NAME;
    output << sp.getType() << "\n";
    for (unsigned int i = 0; i < sp.size(); i++)
        sp[i].save(*this);
}

// -*- C++ -*-----------------------------------------------------------------
//
//  Method:        SeqSaver::saveLigand
//
//  Description:
//    Saves a Ligand in SEQ format. 
//
// ----------------------------------------------------------------------------

void SeqSaver::saveLigand(Ligand& l) {
    PRINT_NAME;
    ERROR("Not implemented yet", exception);
}
