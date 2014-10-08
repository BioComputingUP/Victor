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
/**@mainpage Basic information
 *

 * The Victor library (Virtual Construction Toolkit for Proteins) is an open-source project dedicated to providing tools for analyzing and manipulating protein structures. 
 * This Doxygen documentation will help you to find all the needed detail of the code sourceof the Victor package. 
 *  
 * For a general introduction visit the following wiki http://protein.bio.unipd.it/victor
 *
 *  @namespace Biopool
 *  @brief This library contains methods that will allow you to representate one aminoacidic string in a efficient way. 
 *    
 *  @Description Including the ableness of reading the lineal sequence of aminoacids or to process one PDB structure.
 */




// Includes:
#include <AminoAcid.h>
#include <Debug.h>
#include <limits.h>
#include <float.h>
#include <AtomCode.h>

// Global constants, typedefs, etc. (to avoid):
using namespace Biopool;
double AminoAcid::BOND_ANGLE_N_TO_CB = 110.5;
double AminoAcid::BOND_ANGLE_CB_TO_C = 110.1;
double AminoAcid::BOND_LENGTH_CA_TO_CB = 1.53;

/**
 *@Description Constructor
 *@param none
 */
AminoAcid::AminoAcid() : Group(1, 1), type(XXX), state(COIL), phi(999),
psi(999), omega(999), sideChain(), icc() {
    sideChain.setBackboneRef(this);
}

/**
 *@Description Constructor
 *@param AminoAcid
 */
AminoAcid::AminoAcid(const AminoAcid& orig) {
    PRINT_NAME;
    this->copy(orig);
}

/**
 *@Description DESTRUCTOR
 *@param none
 */
AminoAcid::~AminoAcid() {
    PRINT_NAME;
}

/**
 *@Description Returns the atom corresponding to N, 
 * aminoacids have only a single possible, hard coded, open in-bond
 *@param unsigned int 
 *@return const Atom& 
 */
const Atom& AminoAcid::getOpenInBondRef(unsigned int n) const {
    PRECOND((n == 0) && (sizeOpenInBonds() > 0), exception);
    return (*this)[N];
}

/**
 *@Description Returns the atom corresponding to N, 
 * aminoacids have only a single possible, hard coded, open in-bond
 *@param unsigned int 
 *@return  Atom& 
 */
Atom& AminoAcid::getOpenInBondRef(unsigned int n) {
    PRECOND((n == 0) && (sizeOpenInBonds() > 0), exception);
    return (*this)[N];
}

/**
 *@Description Returns the atom corresponding to C, 
 * aminoacids have only a single possible, hard coded, open in-bond
 *@param unsigned int 
 *@return const Atom& 
 */
const Atom& AminoAcid::getOpenOutBondRef(unsigned int n) const {
    // NB: aminoacids have only a single possible, hard coded, open out-bond
    PRECOND((n == 0) && (sizeOpenOutBonds() > 0), exception);
    return (*this)[C];
}

/**
 *@Description Returns the atom corresponding to C, 
 * aminoacids have only a single possible, hard coded, open in-bond
 *@param unsigned int 
 *@return  Atom& 
 */
Atom& AminoAcid::getOpenOutBondRef(unsigned int n) {
    // NB: aminoacids have only a single possible, hard coded, open out-bond
    PRECOND((n == 0) && (sizeOpenOutBonds() > 0), exception);
    return (*this)[C];
}



// MODIFIERS:

/**
 *@Description Connects the aa to the reverse structure
 *@param AminoAcid* , unsigned int (use 0)
 *@return  void 
 */
void AminoAcid::connectIn(AminoAcid* a, unsigned int offset) {
    if (!(a->isMember(OXT)))
        a->addTerminalOXT();
    IntCoordConverter icc;
    a->setOmega(180.0);
    icc.connectReverseStructure((*this), (*a));
    sync();
    resetBoundaries();
}
// MODIFIERS:

/**
 *@Description Connects the aa to the structure  
 *@param AminoAcid* , unsigned int (use 0)
 *@return  void 
 */
void AminoAcid::connectOut(AminoAcid* a, unsigned int offset) {
    if (!isMember(OXT))
        addTerminalOXT();
    IntCoordConverter icc;
    //    if (offset % 2 == 1)
    setOmega(180.0);
    icc.connectStructure((*a), (*this));
    resetBoundaries();
}

/**
 *@Description unconnect aminoacid from predecessor
 *@param none
 *@return  AminoAcid*
 */
AminoAcid* AminoAcid::unconnectIn() {
    if (!sizeInBonds()) {
        DEBUG_MSG("Cannot unconnect aminoacid without predecessor.");
        return NULL;
    }
    AminoAcid* tmp = &getInBond(0);

    // set OXT for tmp:
    Atom at;
    at.setCode(OXT);
    at.setTrans((*this)[N].getTrans());
    at.bindIn((*tmp)[C]);
    tmp->addAtom(at);

    // unbind the two structures and reset the positions:
    unbindIn(*tmp);
    (*this)[N].unbindIn((*tmp)[C]);
    setTrans((*tmp)[C].getCoords());
    tmp->setOmega(999);

    resetBoundaries();
    tmp->resetBoundaries();
    return tmp;
}

/**
 *@Description unconnect aminoacid from the follower  
 *@param none
 *@return  AminoAcid*
 */
AminoAcid* AminoAcid::unconnectOut() {
    if (!sizeOutBonds()) {
        DEBUG_MSG("Cannot unconnect aminoacid without follower.");
        return NULL;
    }
    AminoAcid* tmp = &getOutBond(0);

    // set OXT for this:
    Atom at;
    at.setCode(OXT);
    at.setTrans((*tmp)[N].getTrans());
    at.bindIn((*this)[C]);
    addAtom(at);
    // unbind the two structures and reset the positions:
    unbindOut(*tmp);
    (*this)[C].unbindOut((*tmp)[N]);
    tmp->setTrans((*this)[C].getCoords());
    setOmega(999);

    resetBoundaries();
    tmp->resetBoundaries();
    return tmp;
}

/**
 *@Description Sets phi angle
 *@param double a (angle value)
 *@return  void
 */
void AminoAcid::setPhi(double a) {
    if (a >= 990) {
        phi = a;
        return;
    }

    if (a < -180)
        a += 360;
    else if (a > 180)
        a -= 360;

    PRECOND((a >= -180) && (a <= 180), exception);
    if (sizeInBonds())
        icc.setTorsionAngle(getInBond(0)[C], (*this)[N], (*this)[CA],
            (*this)[C], DEG2RAD * a);
    phi = a;
    sync();
}

/**
 *@Description Sets Psi angle
 *@param double a (angle value)
 *@return  void
 */
void AminoAcid::setPsi(double a) {
    if (a >= 990) {
        psi = a;
        return;
    }

    if (a < -180)
        a += 360;
    else if (a > 180)
        a -= 360;

    PRECOND((a >= -180) && (a <= 180), exception);
    if (sizeOutBonds())
        icc.setTorsionAngle((*this)[N], (*this)[CA], (*this)[C],
            getOutBond(0)[N], DEG2RAD * a);
    psi = a;
    sync();
}

/**
 *@Description Sets Omega angle
 *@param double a (angle value)
 *@return  void
 */
void AminoAcid::setOmega(double a) {
    if (a >= 990) {
        omega = a;
        return;
    }

    if (a < -180)
        a += 360;
    else if (a > 180)
        a -= 360;

    PRECOND((a >= -180) && (a <= 180), exception);
    if (sizeOutBonds())

        icc.setTorsionAngle((*this)[CA], (*this)[C], getOutBond(0)[N],
            getOutBond(0)[CA], DEG2RAD * a);
    omega = a;
    sync();
}

/**
 *@Description Sets side chain
 *@param SideChain& 
 *@return  void
 */
void AminoAcid::setSideChain(SideChain& sc) {
    if (sc.getType() == "GLY")
        return;

    if (sideChain.size() == 0)
        patchBetaPosition();

    if (sideChain.size() == 0)
        ERROR("Sidechain undefined.", exception);

    vgVector3<double> tmpV = sideChain[0].getTrans();
    vgMatrix3<double> res(1);

    sideChain = sc;
    sideChain.setBackboneRef(this);

    alignVectors(tmpV, sideChain[0].getTrans(), res);

    sideChain[0].addRot(res);
    sync();
}

/**
 *@Description Construct side chain
 *@param SideChain& sc, double chi1, double chi2, 
                              double chi3, double chi4, double chi5 (values for chi, 
 *                              as much as needed, max 5)
 *@return  void
 */
void AminoAcid::constructSideChain(SideChain& sc, double chi1, double chi2,
        double chi3, double chi4, double chi5) {
    setSideChain(sc);
    unsigned int maxChi = getMaxChi();
    switch (maxChi) {
        case 1:
        {
            //control temporarily deactivated to allow computation in some erroneus cases 
            if (chi2 != 999) {
                ERROR("Trying to set more side chain torsion angles than available in current aminoacid", exception);
            }
            setChi(0, chi1);
            break;
        }
        case 2:
        {
            if (chi3 != 999) {
                ERROR("Trying to set more side chain torsion angles than available in current aminoacid", exception);
            }
            setChi(0, chi1);
            setChi(1, chi2);
            break;
        }
        case 3:
        {
            if (chi4 != 999) {
                ERROR("Trying to set more side chain torsion angles than available in current aminoacid", exception);
            }
            setChi(0, chi1);
            setChi(1, chi2);
            setChi(2, chi3);
            break;
        }
        case 4:
        {
            if (chi5 != 999) {
                ERROR("Trying to set more side chain torsion angles than available in current aminoacid", exception);
            }
            setChi(0, chi1);
            setChi(1, chi2);
            setChi(2, chi3);
            setChi(3, chi4);
            break;
        }
        case 5:
        {
            setChi(0, chi1);
            setChi(1, chi2);
            setChi(2, chi3);
            setChi(3, chi4);
            setChi(4, chi5);
            break;
        }
    }
    sync();
}

/**
 *@Description Construct side chain
 *@param SideChain& sc, vector<double> chi (values for chi, 
 *                              as much as needed, max 5)
 *@return  void
 */
void AminoAcid::constructSideChain(SideChain& sc, vector<double> chi) {
    setSideChain(sc);

    if (chi.size() > getSideChain().getMaxChi())
        ERROR("Trying to set more side chain torsion angles than available in current aminoacid", exception);

    setChi(chi);
    sync();
}

/**
 *@Description Sets the state from torsion angles from definition was taken from McGuffin et al.,
  // Bioinformatics (17):63-72 (2001)
 *@param none
 *@return  void
 */
void AminoAcid::setStateFromTorsionAngles() {

    double phi = getPhi();
    double psi = getPsi();

    if (((phi >= -80) && (phi <= -40) && (psi >= -65) && (psi <= -5))
            || ((phi >= -110) && (phi <= -40) && (psi >= -74) && (psi <= -0)))
        setState(HELIX);
    else if (((phi >= -180) && (phi <= -60) && (psi >= 60) && (psi <= 180))
            || ((phi >= -180) && (phi <= -60) && (psi >= -180) && (psi <= -140)))
        setState(STRAND);

}

/**
 *@Description adjusts the translation of the N atom, set trans for N relative to CA & C,
 *         adjust this' translation, the phi angle in leading structures is always undefined
 *@param none
 *@return  void
 */
void AminoAcid::adjustLeadingN() {
    if ((*this)[N].getTrans().length() != 0)
        return;
    const double BOND_ANGLE_N_TO_CA = DEG2RAD * 116.5;
    const double BOND_LENGTH_N_TO_CA = 1.45;
    vgVector3<double> normal = (*this)[CA].getTrans().cross(
            (*this)[C].getTrans());
    vgMatrix3<double> res = vgMatrix3<double>::createRotationMatrix(
            normal, (DEG2RAD * 180.0) - BOND_ANGLE_N_TO_CA);
    (*this)[N].setTrans(BOND_LENGTH_N_TO_CA
            * (res * ((*this)[CA].getTrans())).normalize());
    addTrans(-(*this)[N].getTrans());
}

/**
 *@Description adds an OXT atom to the C terminus
 *@param none
 *@return  void
 */
void AminoAcid::addTerminalOXT() { // adds an OXT atom to the C terminus
    if (isMember(OXT))
        return;

    if ((*this)[C].sizeOutBonds() > 1) {
        DEBUG_MSG("Cannot add terminal OXT on bound aminoacids.");
        return;
    }

    const double BOND_ANGLE_CA_TO_OXT = DEG2RAD * 117.0;
    const double BOND_ANGLE_O_TO_N = 122.0;
    const double BOND_LENGTH_C_TO_OXT = 1.35;
    Atom at;
    at.setType("OXT");
    at.setBFac(15.0);
    at.bindIn((*this)[C]);
    addAtom(at);
    IntCoordConverter icc;
    icc.zAtomToCartesian((*this)[C], BOND_LENGTH_C_TO_OXT, (*this)[CA],
            BOND_ANGLE_CA_TO_OXT, (*this)[O],
            BOND_ANGLE_O_TO_N, 1, (*this)[OXT]);
    sync();
}

/**
 *@Description  adds an O atom, if missing
 *@param none
 *@return  void
 */
void AminoAcid::addMissingO() {
    if (isMember(O))
        return;

    if (sizeOutBonds() < 1)
        return;

    const double BOND_ANGLE_CA_TO_O = 120.80;
    const double BOND_ANGLE_N_TO_O = 123.0;
    const double BOND_LENGTH_C_TO_O = 1.231;

    Atom at;
    at.setType("O");
    at.bindIn((*this)[C]);
    addAtom(at);


    IntCoordConverter icc;
    icc.zAtomToCartesian((*this)[C], BOND_LENGTH_C_TO_O, (*this)[CA],
            BOND_ANGLE_CA_TO_O, getOutBond(0)[N],
            BOND_ANGLE_N_TO_O, 1, (*this)[O]);
    sync();
}

/**
 *@Description  removeHAtomsfromLeadingNH3,
 * check for NH3+ and, if so, remove superfluous H (i.e. 2H, 3H) atoms
 *@param none
 *@return  void
 */
void AminoAcid::removeHAtomsfromLeadingNH3() {
    if (!isMember(H))
        return;

    bool first = false;
    for (unsigned int i = 0; i < (*this)[N].sizeOutBonds(); i++) {
        if ((*this)[N].getOutBond(i).getCode() == H) {
            if (!first)
                first = true;
            else { // remove H atom
                Atom tmp = (*this)[N].getOutBond(i);
                (*this)[N].unbindOut(tmp);
                (*this).removeAtom(tmp);
            }
        }
    }
}

/**
 *@Description  adds a CB atom to the sidechain, if necessary
 *@param unsigned int n (value of number to set)
 *@return  void
 */
void AminoAcid::patchBetaPosition(unsigned int n) {
    if (type == GLY)
        return;

    Atom at;
    at.setType("CB");

    if (n > 0)
        at.setNumber(n);

    getSideChain().addAtom(at);

    (*this)[CA].bindOut(getSideChain()[CB]);

    IntCoordConverter icc;
    icc.zAtomToCartesian((*this)[CA], BOND_LENGTH_CA_TO_CB,
            (*this)[N], BOND_ANGLE_N_TO_CB,
            (*this)[C], BOND_ANGLE_CB_TO_C,
            1, (*this).getSideChain()[CB]);
    sync();
}

/**
 *@Description  Copies an aa
 *@param const AminoAcid& (copy from the orig)
 *@return  void
 */
void AminoAcid::copy(const AminoAcid& orig) {
    PRINT_NAME;
    Group::copy(orig);

    type = orig.type;
    phi = orig.phi;
    psi = orig.psi;
    omega = orig.omega;

    sideChain = orig.sideChain;

    if ((orig.getSideChain().size() > 0)
            && (orig[CA].isBond(orig.getSideChain()[0])))
        sideChain.setBackboneRef(this);

    // fix Proline CD to N bond:
    if ((getType() == "PRO")) {
        (*this)[N].setMaxInBonds(2);
        if (getSideChain().isMember(CD) && orig.getSideChain().isMember(CD) &&
                (orig[N].isBond(orig[CD])))
            getSideChain()[CD].bindOut((*this)[N]);
    }

    // set absolute position to orig's:
    if (orig[0].sizeInBonds()) {
        setTrans(const_cast<AminoAcid&> (orig)[0].getInBond(0).getCoords());
    }
}

/**
 *@Description  Syncronize and sets boundaries
 *@param none
 *@return  void
 */
void AminoAcid::sync() {
    Group::sync();
    sideChain.sync();
    resetBoundaries();
}

/**
 *@Description  Clone the aa
 *@param none
 *@return  Component* 
 */
Component* AminoAcid::clone() {
    AminoAcid* tmp = new AminoAcid;
    tmp->copy(*this);
    return tmp;
}

// OPERATORS:

/**
 *@Description  Operator =, assign the aa
 *@param AminoAcid reference
 *@return  AminoAcid
 */
AminoAcid& AminoAcid::operator=(const AminoAcid& orig) {
    PRINT_NAME;
    if (&orig != this)
        copy(orig);
    return *this;
}

// HELPERS:

/**
 *@Description  Reset the boundaries 
 *@param none
 *@return  void
 */
void AminoAcid::resetBoundaries() {
    Group::resetBoundaries();
    sideChain.resetBoundaries();

    for (unsigned int i = 0; i < 3; i++) {
        if (sideChain.getLowerBound()[i] < lowerBound[i])
            lowerBound[i] = sideChain.getLowerBound()[i];
        if (sideChain.getUpperBound()[i] > upperBound[i])
            upperBound[i] = sideChain.getUpperBound()[i];
    }
}

/**
 *@Description  set bond length and angles to cristallographic values
 *@param none
 *@return  void
 */
void AminoAcid::setDefault() {
    const double BOND_ANGLE_CA_TO_O = 120.80;
    const double BOND_LENGTH_C_TO_O = 1.231;
    const double BOND_LENGTH_CALPHA_TO_CPRIME = 1.52;
    const double BOND_LENGTH_N_TO_CALPHA = 1.458;
    const double BOND_ANGLE_AT_CALPHA_TO_CPRIME = 111.6;

    if ((*this)[N].isBond((*this)[CA]) && (*this)[CA].isBond((*this)[C]) && (*this)[C].isBond((*this)[O])) {
        IntCoordConverter icc;

        icc.setBondLength((*this)[N], (*this)[CA], BOND_LENGTH_N_TO_CALPHA);
        icc.setBondLength((*this)[CA], (*this)[C], BOND_LENGTH_CALPHA_TO_CPRIME);
        icc.setBondLength((*this)[C], (*this)[O], BOND_LENGTH_C_TO_O);

        icc.setBondAngle((*this)[N], (*this)[CA], (*this)[C], DEG2RAD * BOND_ANGLE_AT_CALPHA_TO_CPRIME);
        icc.setBondAngle((*this)[CA], (*this)[C], (*this)[O], DEG2RAD * BOND_ANGLE_CA_TO_O);
    } else {
        cout << "not in bond\n";
    }
    sync();
}

/**
 *@Description  Set the bonds based on lengths and angles given as params
 *@param double NToCaLen, double CaToCLen, double CToOLen, double atCaToCAng, double CaToOAng
 *@return  void
 */
void
AminoAcid::setBonds(double NToCaLen, double CaToCLen, double CToOLen, double atCaToCAng, double CaToOAng) {

    icc.setBondLength((*this)[N], (*this)[CA], NToCaLen);
    icc.setBondLength((*this)[CA], (*this)[C], CaToCLen);
    icc.setBondLength((*this)[C], (*this)[O], CToOLen);
    icc.setBondAngle((*this)[N], (*this)[CA], (*this)[C], atCaToCAng);
    icc.setBondAngle((*this)[CA], (*this)[C], (*this)[O], CaToOAng);
    sync();
}

/**
 * @Description  sets the bond structure for an aminoacid, and connect it (if connect is true)
 * @param connect, if True connect the AminoAcid with previous AminoAcid (set "Trans" and "Coords")
 * @param prev, the previous AminoAcid 
 * @param permissive, allows to operate on an incomplete AminoAcid
 * @return True, if the AminoAcid has been connected with the previous AminoAcid
 */
bool
AminoAcid::setBondsFromPdbCode(bool connect, AminoAcid* prev, bool permissive) {

    if ((!isMember(C)) || (!isMember(CA)) || (!isMember(N))) {
        if (permissive) {
            return false;
        } else
            ERROR("Cannot set bonds for residues with missing atoms.", exception);
    }

    if (getType() == "X")
        patchAminoAcidCode();

    if (prev != NULL) {
        if (connect)
            (*this)[N].setTrans((*this)[N].getCoords() - (*prev)[C].getCoords());
        bindIn((*this)[N], (*prev), (*prev)[C]);
    }
    (*this)[CA].bindStructure((*this)[N], connect);
    if (isMember(H)) //HN
        (*this)[H].bindStructure((*this)[N], connect);

    /* OBSOLETE, removed by Damiano Piovesan 2014  
    else if (isMember(H))                // HN sometimes called H instead
      {
        for (unsigned int i = 0; i < sizeBackbone(); i++)
          if( (*this)[i].getCode() == H)
            (*this)[i].bindStructure((*this)[N], connect);
      }
     */

    (*this)[C].bindStructure((*this)[CA], connect);
    if (isMember(HA))
        (*this)[HA].bindStructure((*this)[CA], connect);

    /* OBSOLETE, removed by Damiano Piovesan 2014  
    else if (isMember(HA1))
      (*this)[HA1].bindStructure((*this)[CA], connect);
     */


    if (isMember(O))
        (*this)[O].bindStructure((*this)[C], connect);

    if (isMember(OXT))
        (*this)[OXT].bindStructure((*this)[C], connect);

    return sideChain.setBondsFromPdbCode(connect, permissive);
}
