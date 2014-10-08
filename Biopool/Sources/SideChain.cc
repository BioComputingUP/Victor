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
/**
 *@Description:
 *    This class implements a simple side chain.
 */

// Includes:
#include <string>
#include <SideChain.h>
#include <AminoAcid.h>
#include <limits.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:

SideChain::SideChain() : Group(0, 0), maxChi(0), chi(), backboneRef(NULL) {
}

SideChain::SideChain(const SideChain& orig) {
    this->copy(orig);
}

SideChain::~SideChain() {
    PRINT_NAME;
}

// PREDICATES:

vector<double> SideChain::getChi() {
    for (unsigned int i = 0; i < chi.size(); i++)
        if (chi[i] > 990)
            chi[i] = RAD2DEG * calculateChi(i);
    return chi;
}

// MODIFIERS:

void SideChain::setChi(vector<double> cv) {
    PRECOND(cv.size() <= chi.size(), exception);

    for (unsigned int i = 0; i < cv.size(); i++) {
        INVARIANT((cv[i] >= -180) && (cv[i] <= 180), exception);
        convertChi(i, DEG2RAD * cv[i]);
        chi[i] = cv[i];
    }
}

void SideChain::copy(const SideChain& orig) {
    PRINT_NAME;
    Group::copy(orig);
    maxChi = orig.maxChi;
    chi = orig.chi;
    backboneRef = orig.backboneRef; // orig keeps its backbone reference
}

void SideChain::setBackboneRef(AminoAcid* br) {
    backboneRef = br;
    if ((backboneRef == NULL) || (size() == 0))
        return;
    if ((*this)[0].sizeInBonds())
        (*this)[0].unbindIn((*this)[0].getInBond(0));
    if (backboneRef->isMember(CA))
        (*this)[0].bindIn((*backboneRef)[CA]);
}

void SideChain::patchAminoAcidCode() { // determines the 3-letter code from the sidechain
    if (!isMember(CB)) {
        if (size() > 1)
            ERROR("SideChain::patchAminoAcidCode() : Unknown sidechain.",
                exception);
        setType("GLY");
        return;
    }

    if (isMember(SG)) {
        setType("CYS");
        return;
    }

    if (isMember(CG)) { // leu, phe, tyr, trp, met, asp, glu, asn, gln, lys, arg, his, pro
        if (isMember(SD))
            setType("MET");
        else if (isMember(CD)) { // glu, gln, lys, arg, pro
            if (isMember(OE1)) { // glu, gln
                if (isMember(OE2))
                    setType("GLU");
                else
                    setType("GLN");
            } else { // lys, arg, pro
                if (isMember(NE))
                    setType("ARG");
                else if (isMember(CE))
                    setType("LYS");
                else
                    setType("PRO");
            }
        } else { // leu, phe, tyr, trp, asp, asn, his
            if (isMember(OD1)) { // asp, asn
                if (isMember(OD2))
                    setType("ASP");
                else
                    setType("ASN");
            } else { // leu, phe, tyr, trp, his
                if (isMember(ND1))
                    setType("HIS");
                else if (isMember(NE1))
                    setType("TRP");
                else { // leu, phe, tyr
                    if (!isMember(CE1))
                        setType("LEU");
                    else if (isMember(OH))
                        setType("TYR");
                    else
                        setType("PHE");
                }
            }
        }
    } else if (isMember(CG1)) {//  ile, thr, val
        if (isMember(OG1))
            setType("THR");
        else if (isMember(CD1))
            setType("ILE");
        else
            setType("VAL");
    } else {// ser, ala 

        if (isMember(OG))
            setType("SER");
        else
            setType("ALA");
    }

}

bool SideChain::setBondsFromPdbCode(bool connect, bool permissive) {
    if (getType() == "X")
        patchAminoAcidCode();

    string seq = (getType() == "TRP") ? "BG" : "BGDEZH";
    // sequence of PDB codes, TRP is a special case (see below)
    unsigned int pos, pos2, pos3, pos4;

    // special case: no sidechain atoms
    if (size() == 0)
        return true;

    // special case sidechain pseudoatoms
    // special case Glycine
    if ((getType() == "GLY") || (size() == 1)) {
        if (size() && (backboneRef != NULL))
            (*this)[0].bindStructure((*backboneRef)[CA], connect);
        return true;
    }

    if (isMember(CB) && (backboneRef != NULL)) // connect CB with backbone
        (*this)[CB].bindStructure((*backboneRef)[CA], connect);


    for (unsigned int i = 0; i < seq.length(); i++) {
        if (findCode(seq[i], pos, pos2)) {

            // bind H Atoms
            if (pos2 != INT_MAX) {
                findHCode(seq[i], '1', (*this)[pos], connect);
                findHCode(seq[i], '2', (*this)[pos2], connect);
            } else
                findHCode(seq[i], ' ', (*this)[pos], connect);
            if (findCode(seq[i + 1], pos3, pos4)) {// chain continues?

                (*this)[pos3].bindStructure((*this)[pos], connect);
                if (pos4 != INT_MAX) {
                    if (pos2 == INT_MAX)
                        (*this)[pos4].bindStructure((*this)[pos], connect);
                        // chain bifurcation

                    else
                        (*this)[pos4].bindStructure((*this)[pos2], connect);
                    // parallel chains
                }
            }
        } else
            break;
    }

    setMaxChiFromCode();

    // special cases: 
    if ((getType() == "TYR") || (getType() == "PHE")) {
        if (isMember(CZ) and isMember(CE2)) {
            (*this)[CZ].setMaxInBonds(2);
            (*this)[CZ].bindIn((*this)[CE2]);
        }
        return true;
    }

    if (getType() == "PRO") {
        if (isMember(N) and isMember(CD)) {
            (*backboneRef)[N].setMaxInBonds(2);
            (*backboneRef)[N].bindIn((*this)[CD]);
        }
        return true;
    }

    if (getType() == "HIS") {
        if (isMember(NE2) and isMember(CE1)) {
            (*this)[NE2].setMaxInBonds(2);
            (*this)[NE2].bindIn((*this)[CE1]);
        }
        return true;
    }

    if (getType() == "ILE") {
        if (isMember(HD11)and isMember(CD1)) // setting H-bonds
            (*this)[HD11].bindStructure((*this)[CD1], connect);
        if (isMember(HD12)and isMember(CD1))
            (*this)[HD12].bindStructure((*this)[CD1], connect);
        if (isMember(HD13)and isMember(CD1))
            (*this)[HD13].bindStructure((*this)[CD1], connect);
    }

    if (getType() == "TRP") { // TRP is a complex special case that has to be built manually 
        if (isMember(CD1)and isMember(CG)) {
            (*this)[CD1].bindStructure((*this)[CG], connect);
        }
        if (isMember(CD2)and isMember(CG)) {
            (*this)[CD2].bindStructure((*this)[CG], connect);
        }
        if (isMember(NE1)and isMember(CD1)) {
            (*this)[NE1].bindStructure((*this)[CD1], connect);
        }
        if (isMember(CE2)and isMember(CD2)) {
            (*this)[CE2].bindStructure((*this)[CD2], connect);
            (*this)[CE2].setMaxInBonds(2);
        }
        if (isMember(CE2)and isMember(NE1)) {
            (*this)[CE2].bindIn((*this)[NE1]);
        }
        if (isMember(CE3)and isMember(CD2)) {
            (*this)[CE3].bindStructure((*this)[CD2], connect);
        }
        if (isMember(CZ2)and isMember(CE2)) {
            (*this)[CZ2].bindStructure((*this)[CE2], connect);
        }
        if (isMember(CZ3)and isMember(CE3)) {
            (*this)[CZ3].bindStructure((*this)[CE3], connect);
        }
        if (isMember(CH2)and isMember(CZ2)) {
            (*this)[CH2].bindStructure((*this)[CZ2], connect);
        }
        if (isMember(CH2)and isMember(CZ3)) {
            (*this)[CH2].setMaxInBonds(2);
            (*this)[CH2].bindIn((*this)[CZ3]);
        }

        // setting H-bonds
        if (isMember(HD1))
            (*this)[HD1].bindStructure((*this)[CD1], connect);
        if (isMember(HE1))
            (*this)[HE1].bindStructure((*this)[NE1], connect);
        if (isMember(HE3))
            (*this)[HE3].bindStructure((*this)[CE3], connect);
        if (isMember(HZ2))
            (*this)[HZ2].bindStructure((*this)[CZ2], connect);
        if (isMember(HZ3))
            (*this)[HZ3].bindStructure((*this)[CZ3], connect);
        if (isMember(HH2))
            (*this)[HH2].bindStructure((*this)[CH2], connect);
    }

    return true;
}

bool SideChain::findCode(char c, unsigned int &pos, unsigned int &pos2) {
    pos = INT_MAX;
    pos2 = INT_MAX;
    for (unsigned int i = 0; i < (*this).size(); i++) {
        string tmp;
        tmp = (*this)[i].getType();
        unsigned int shift = ((tmp[0] == '1') || (tmp[0] == '2')
                || (tmp[0] == '3')) ? 1 : 0;
        if ((tmp[1 + shift] == c) && (tmp[0 + shift] != 'H')) {
            if (pos != INT_MAX) {
                pos2 = i;
                break;
            } else
                pos = i;
        }
    }
    return ((pos != INT_MAX) ? true : false);
}

void SideChain::findHCode(char c, char c2, Atom& at, bool connect) {
    for (unsigned int i = 0; i < (*this).size(); i++) {
        string tmp;
        tmp = (*this)[i].getType() + ' '; // simplifies comparison with c2
        unsigned int shift = ((tmp[0] == '1') || (tmp[0] == '2')
                || (tmp[0] == '3')) ? 1 : 0;
        if ((tmp[1 + shift] == c) && (tmp[2 + shift] == c2) && (tmp[0 + shift] == 'H'))
            (*this)[i].bindStructure(at, connect);
    }
}


// OPERATORS:

SideChain& SideChain::operator=(const SideChain& orig) {
    PRINT_NAME;
    if (&orig != this)
        copy(orig);
    return *this;
}

// HELPERS:

double SideChain::calculateChi(unsigned int n) {
    PRECOND((n < 5), exception);
    IntCoordConverter icc;
    if (size() > 1) // size() == 1 if sidechain has only pseudoatom
        switch (n) {
            case 0: return icc.getTorsionAngle((*backboneRef)[N],
                        (*backboneRef)[CA],
                        (*this)[CB], firstG());
            case 1: return icc.getTorsionAngle((*backboneRef)[CA], (*this)[CB],
                        firstG(), firstD());
            case 2: return icc.getTorsionAngle((*this)[CB], firstG(),
                        firstD(), firstE());
            case 3: return icc.getTorsionAngle(firstG(), firstD(),
                        firstE(), firstZ());
            case 4: return icc.getTorsionAngle(firstD(), firstE(),
                        firstZ(), firstH());
        }
    return 999.0 * DEG2RAD;
}

void SideChain::convertChi(unsigned int n, double c) {
    PRECOND((n < 5), exception);
    IntCoordConverter icc;
    switch (n) {
        case 0: return icc.setTorsionAngle((*backboneRef)[N],
                    (*backboneRef)[CA],
                    (*this)[CB], firstG(), c);
        case 1: return icc.setTorsionAngle((*backboneRef)[CA], (*this)[CB],
                    firstG(), firstD(), c);
        case 2: return icc.setTorsionAngle((*this)[CB], firstG(),
                    firstD(), firstE(), c);
        case 3: return icc.setTorsionAngle(firstG(), firstD(),
                    firstE(), firstZ(), c);
        case 4: return icc.setTorsionAngle(firstD(), firstE(),
                    firstZ(), firstH(), c);
    }
    sync();
}

bool SideChain::hasB() {
    if (isMember(CB))
        return true;
    return false;
}

Atom& SideChain::firstG() {
    if (isMember(CG))
        return (*this)[CG];
    if (isMember(CG1))
        return (*this)[CG1];
    if (isMember(OG))
        return (*this)[OG];
    if (isMember(OG1))
        return (*this)[OG1];
    if (isMember(SG))
        return (*this)[SG];

    for (unsigned int i = 0; i < size(); i++)
        cout << (*this)[i].getType() << "\t";

    ERROR("No G-Atom found.", exception);

    return (*this)[0];
}

bool SideChain::hasG() {
    if (isMember(CG) || isMember(CG1) || isMember(OG) ||
            isMember(OG1) || isMember(SG))
        return true;
    return false;
}

Atom& SideChain::firstD() {
    if (isMember(CD))
        return (*this)[CD];
    if (isMember(CD1))
        return (*this)[CD1];
    if (isMember(ND1))
        return (*this)[ND1];
    if (isMember(OD))
        return (*this)[OD];
    if (isMember(OD1))
        return (*this)[OD1];
    if (isMember(SD))
        return (*this)[SD];
    ERROR("No D-Atom found.", exception);
    return (*this)[0];
}

bool SideChain::hasD() {
    if (isMember(CD) || isMember(CD1) || isMember(ND1) ||
            isMember(OD) || isMember(OD1) || isMember(SD))
        return true;
    return false;
}

Atom& SideChain::firstE() {
    if (isMember(CE))
        return (*this)[CE];
    if (isMember(NE))
        return (*this)[NE];
    if (isMember(CE1))
        return (*this)[CE1];
    if (isMember(NE1))
        return (*this)[NE1];
    if (isMember(OE1))
        return (*this)[OE1];
    ERROR("No E-Atom found.", exception);
    return (*this)[0];
}

bool SideChain::hasE() {
    if (isMember(CE) || isMember(NE) ||
            isMember(CE1) || isMember(NE1) || isMember(OE1))
        return true;
    return false;
}

Atom& SideChain::firstZ() {
    if (isMember(CZ))
        return (*this)[CZ];
    if (isMember(NZ))
        return (*this)[NZ];
    if (isMember(CZ2))
        return (*this)[CZ2];
    ERROR("No Z-Atom found.", exception);
    return (*this)[0];
}

bool SideChain::hasZ() {
    if (isMember(CZ) || isMember(NZ) || isMember(CZ2))
        return true;
    return false;
}

Atom& SideChain::firstH() {
    if (isMember(OH))
        return (*this)[OH];
    if (isMember(NH1))
        return (*this)[NH1];
    if (isMember(CH2))
        return (*this)[CH2];
    ERROR("No H-Atom found.", exception);
    return (*this)[0];
}

bool SideChain::hasH() {
    if (isMember(OH) || isMember(NH1) || isMember(CH2))
        return true;
    return false;
}

void SideChain::setMaxChiFromCode() {
    unsigned int tMaxChi = 0;
    switch (getCode()) {
        case XXX: case ALA:
        case GLY:
            tMaxChi = 0;
            break;

        case SER: case THR:
        case CYS: case VAL:
            tMaxChi = 1;
            break;

        case ASP:
        case PHE: case HIS:
        case ILE: case LEU:
        case ASN: case PRO:
        case TRP: case TYR:
            tMaxChi = 2;
            break;

        case GLU: case MET:
        case GLN:
            tMaxChi = 3;
            break;

        case LYS:
            tMaxChi = 4;
            break;

        case ARG:
            tMaxChi = 5;
            break;

        case AminoAcid_CODE_SIZE:
            ERROR("AminoAcidTranslator(AminoAcidCode code): unknown code", exception);
    }

    // check for sidechain integrity, reduce torsion anlges if necessary:
    if ((tMaxChi >= 1) && (!hasB() || !hasG()))
        tMaxChi = 0;
    if ((tMaxChi >= 2) && !hasD())
        tMaxChi = 1;
    if ((tMaxChi >= 3) && !hasE())
        tMaxChi = 2;
    if ((tMaxChi >= 4) && !hasZ())
        tMaxChi = 3;
    if ((tMaxChi >= 5) && !hasH())
        tMaxChi = 4;

    setMaxChi(tMaxChi);
}

void SideChain::setMaxChi(unsigned int _maxChi) {
    maxChi = _maxChi;
    chi.resize(_maxChi);
    for (unsigned int i = 0; i < _maxChi; i++)
        chi[i] = 999;
}



