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
 *    This class implements a "Spacer" or protein fragment. 
 *  ***Attention***: The current implementation allows for "1 to 1" spacers,
 *    ie. spacers composed of a single aminoacid chain. Both spacers composed 
 *    of more than one aminoacid chain (eg. 2 beta sheet parts).
 *
 */
// Includes:
#include <Spacer.h>
#include <Debug.h>
#include <IntCoordConverter.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <vector3.h>
#include <matrix3.h>
#include <iostream>
#include <PdbLoader.h>
#include <cassert>
using namespace std;

using namespace Biopool;

// Global constants, typedefs, etc. (to avoid):
double Spacer::BOND_ANGLE_AT_CPRIME_TO_N = 116.4;
double Spacer::BOND_LENGTH_CPRIME_TO_N = 1.33;
double Spacer::BOND_ANGLE_CA_TO_O = 120.80;
double Spacer::BOND_LENGTH_C_TO_O = 1.231;
double Spacer::BOND_LENGTH_CALPHA_TO_CPRIME = 1.52;
double Spacer::BOND_LENGTH_N_TO_CALPHA = 1.458;
double Spacer::BOND_ANGLE_AT_CALPHA_TO_CPRIME = 111.6;
double Spacer::BOND_ANGLE_O_TO_N = 123.2;
double Spacer::BOND_ANGLE_AT_N_TO_CALPHA = 121.9;

// CONSTRUCTORS/DESTRUCTOR:

/**
 *@Description Basic constructor
 */
Spacer::Spacer() : Polymer(1, 1), startOffset(0), startAtomOffset(0), gaps(),
subSpacerList() {
    PRINT_NAME;
}

/**
 *@Description constructor based in another object
 *@param Reference to the original object to copy(const Spacer&)
 */
Spacer::Spacer(const Spacer& orig) : subSpacerList() {
    PRINT_NAME;
    this->copy(orig);
}

/**
 *@Description Basic destructor
 */
Spacer::~Spacer() {
    PRINT_NAME;
}

// PREDICATES:

/**
 *@Description Predicate used to get the internal array number corresponding to 
 * a given PDB aminoacid number (Residue sequence number, columns 23-26 of a pdb file).
 *@param  pdb number (int) 
 *@return  the corresponding internal array number(int)
 */
int Spacer::getIndexFromPdbNumber(int index) {
    int tmpIndex = index;

    for (unsigned int i = 0; i < gaps.size(); i++) {
        if (index < gaps[i])
            break;
        tmpIndex--;
    }

    tmpIndex -= startOffset;

    if ((tmpIndex < 0) || (static_cast<unsigned int> (tmpIndex) >= sizeAmino())) {
        cerr << "index= " << tmpIndex << "\t startOffset= " << startOffset
                << "\t sizeAmino()= " << sizeAmino() << "\n";
        ERROR("Invalid PDB number requested.", exception);
    }

    return tmpIndex;
}

/**
 *@Description Predicate used to get the PDB aminoacid number (Residue sequence number, columns 23-26 of a pdb file)  corresponding to 
 * a given internal array number
 *@param internal array number (int) 
 *@return  the corresponding PDB aminoacid number(int)
 */
int Spacer::getPdbNumberFromIndex(int index) {
    if (index > static_cast<int> (sizeAmino())) {
        cerr << "index= " << index << "\t sizeAmino()= " << sizeAmino() << "\n";
        ERROR("Invalid PDB number requested.", exception);
    }

    int ind = startOffset + index;

    for (unsigned int i = 0; i < gaps.size(); i++) {
        if (ind < gaps[i])
            break;
        ind++;
    }

    return ind;
}

/**
 *@Description  Predicate used to determine if there is a gap in a given position/PDB aminoacid number 
 *@param  given position to evaluate(int) 
 *@return  result from the verification (bool)
 */
bool Spacer::isGap(int index) {

    if (index < startOffset)
        return true;
    for (unsigned int i = 0; i < gaps.size(); i++) {
        if (index == gaps[i])
            return true;
    }

    if (index >= static_cast<int> (sizeAmino() + startOffset + gaps.size()))
        return true;
    else
        return false;
}

/**
 *@Description Predicate used to a list of the positions/pdb amino acid numbers where a gap is present
 *@param none
 *@return  prints the corresponding values
 */
void Spacer::printGaps() {
    cout << "Gaps listing:\t";
    for (unsigned int i = 0; i < gaps.size(); i++)
        cout << gaps[i] << "\t";
    cout << "\n";
}

/**
 *@Description used to obtain the corresponding reference to the InBond of the n amino acid  
 *@param index of the atom (unsigned int )
 *@return  atom reference(const Atom&)
 */

const Atom& Spacer::getOpenInBondRef(unsigned int n) const {
    PRECOND((n == 0) && (sizeOpenInBonds() > 0), exception);
    return getAmino(0)[N];
}

/**
 *@Description used to obtain the corresponding reference to the InBond of the n amino acid 
 *@param index of the atom (unsigned int )
 *@return  atom reference(const Atom&)
 */
Atom& Spacer::getOpenInBondRef(unsigned int n) {
    PRECOND((n == 0) && (sizeOpenInBonds() > 0), exception);
    return getAmino(0)[N];
}

/**
 *@Description used to obtain the corresponding reference to the OutBond of the n amino acid  
 *@param index of the atom (unsigned int )
 *@return  atom reference(const Atom&)
 */
const Atom& Spacer::getOpenOutBondRef(unsigned int n) const {
    PRECOND((n == 0) && (sizeOpenOutBonds() > 0), exception);
    return getAmino(sizeAmino() - 1)[C];
}

/**
 *@Description used to obtain the corresponding reference to the OutBond of the n amino acid  
 *@param index of the atom (unsigned int )
 *@return  atom reference(const Atom&)
 */
Atom& Spacer::getOpenOutBondRef(unsigned int n) {
    PRECOND((n == 0) && (sizeOpenOutBonds() > 0), exception);
    return getAmino(sizeAmino() - 1)[C];
}

/**
 *@Description Returns the count of all aminoacids in the spacer  
 *@param void
 *@return  Quantity of amino acids in the spacer(unsigned int)
 */
const unsigned int Spacer::sizeAmino() const {
    unsigned int count = 0;
    unsigned int subSpacerSize = subSpacerList.size();
    for (unsigned int loop = 0; loop < subSpacerSize; loop++)
        count += subSpacerList[loop].second - subSpacerList[loop].first;
    return ( size() + count);
}

/**
 *@Description Predicate used to get the range for the Helix amino acids
 *@param none
 *@return  vector containing the ranges(vector < pair < unsigned int, unsigned int> >)
 */
vector<pair<unsigned int, unsigned int> > Spacer::getHelixData() {
    vector<pair<unsigned int, unsigned int> > helixData;
    unsigned int i = 0;
    while (i < sizeAmino())
        if (getAmino(i).getState() == HELIX) {
            unsigned int start = i;
            unsigned int end = i;

            while ((i < sizeAmino()) && (getAmino(i).getState() == HELIX)) {
                end++;
                i++;
            }

            helixData.push_back(pair<unsigned int, unsigned int>(start, end));
        } else
            i++;

    return helixData;
}

/**
 *@Description Predicate used to get the range for the Strand amino acids
 *@param none
 *@return  vector containing the ranges(vector < pair < unsigned int, unsigned int> >)
 */
vector<pair<unsigned int, unsigned int> > Spacer::getStrandData() {
    vector<pair<unsigned int, unsigned int> > strandData;
    unsigned int i = 0;
    while (i < sizeAmino())
        if (getAmino(i).getState() == STRAND) {
            unsigned int start = i;
            unsigned int end = i;

            while ((i < sizeAmino()) && (getAmino(i).getState() == STRAND)) {
                end++;
                i++;
            }

            strandData.push_back(pair<unsigned int, unsigned int>(start, end));
        } else
            i++;

    return strandData;
}



// MODIFIERS:

/**
 *@Description Set the start offset, considers the first amino acid number in the pdb file
 *@param  value for the offset(int)
 *@return  changes are made internally(void)
 */
void Spacer::setStartOffset(int _offset) {
    int diff = _offset - startOffset;

    for (unsigned int i = 0; i < gaps.size(); i++)
        gaps[i] += diff;

    startOffset = _offset;
}

/**
 *@Description Modifier used to add a gap into the list of gaps.
 *@param  amino acid number from on the pdb file (int) 
 *@return  changes are made internally (void)
 */
void Spacer::addGap(int index) {
    if (index < startOffset) {
        startOffset = index;
        return;
    }

    if ((gaps.size() == 0) || (index > gaps[gaps.size() - 1]))
        gaps.push_back(index);
    else
        for (unsigned int i = gaps.size() - 1; i >= 0; i--)
            if (index == gaps[i])
                break;
            else
                if (index < gaps[i]) {
                gaps.insert(gaps.begin() + i, 1, index);
                break;
            }
}

/**
 *@Description Modifier used to remove all gaps from the gap list.
 *@param none
 *@return changes are made internally (void)
 * */
void Spacer::removeAllGaps() {
    while (gaps.size() > 0) {
        gaps.erase(gaps.begin());
    }
}

/**
 *@Description Modifier used to remove a specific gap from the gaps list created.
 *@param  amino acid number(int)
 *@return  changes are made internally(void)
 */
void Spacer::removeGap(int index) {
    if (index < startOffset)
        ERROR("Cannot remove a leading gap.", exception);

    for (unsigned int i = 0; i < gaps.size(); i++)
        if (index == gaps[i]) {
            gaps.erase(gaps.begin() + i);
            break;
        }
}

/**
 *@Description This function is used to insert an amino after a GAP. The general
 *       idea is to develop a mixed function between inserting a new amino and connecting
 *       aminoacids in the spacer). 
 *@param  
 *  The difference is
 *       placed into the parameters beginHole & endHole. They represent those
 *       two indices such that (beginHole+1..endHole-1) is an hole. We need
 *       these parameters for deciding how to change the internal numbering
 *@param string code, 	unsigned int n, 	double phi, 	double psi, 	double omega,	int beginHole, 
      int endHole, vgVector3 <double> ca1,  vgVector3 <double> ca2, string target,  Spacer *refSpacer, 
      Spacer* pOriginalSpacer, double NToCaLen, double CaToCLen,  double CToOLen, double atCaToCAng, 	double CaToOAng,
      double CtoNLen, double atCToNAng, double OToNAng, double atNToCaAng
 *@return  void
 */
void Spacer::insertAminoAfterWithGaps(
        string code,
        unsigned int n,
        double phi,
        double psi,
        double omega,
        int beginHole,
        int endHole,
        vgVector3 <double> ca1,
        vgVector3 <double> ca2,
        string target,
        Spacer *refSpacer,
        Spacer* pOriginalSpacer,
        double NToCaLen,
        double CaToCLen,
        double CToOLen,
        double atCaToCAng,
        double CaToOAng,
        double CtoNLen,
        double atCToNAng,
        double OToNAng,
        double atNToCaAng) {
    // -----------------------------------------------------------------------
    // This is only a control code.
    // -----------------------------------------------------------------------
    assert(pOriginalSpacer != 0);
    assert(refSpacer != 0);
    assert(beginHole <= endHole);
    vector< pair<int, int> > arrAllHoles = pOriginalSpacer->getHoles();
    bool bHoleExists = false;
    for (vector< pair<int, int> >::iterator it = arrAllHoles.begin(); it != arrAllHoles.end(); ++it)
        if (it->first == beginHole && it->second == endHole) {
            //cout <<"The hole [" <<beginHole <<"," <<endHole <<"] exists" <<endl;
            bHoleExists = true;
            break;
        }
    if (!bHoleExists)
        cout << "The hole [" << beginHole << "," << endHole << "] does not exist: possible internal numbering error" << endl;
    // +----------------------------------------------------------------------+
    // |GAP management:                                                       |
    // |                                                                      |
    // +----------------------------------------------------------------------+
    for (int i = beginHole + 1; i <= endHole - 1; ++i)
        addGap(i);

    // +----------------------------------------------------------------------+
    // |Now we insert the amino:                                              |
    // |                                                                      |
    // +----------------------------------------------------------------------+

    IntCoordConverter icc;
    AminoAcid* aa = new AminoAcid;


    for (unsigned int i = 0; i < refSpacer->size(); i++)
        if (refSpacer->getAmino(i).getType1L() == threeLetter2OneLetter(code)) {
            *aa = refSpacer->getAmino(i);
            break;
        }

    if (aa == 0)
        ERROR("Internal error insertAminoAfterWithGaps()", exception);
    aa->bindIn((*aa)[N], dynamic_cast<AminoAcid&> (*components[n]), dynamic_cast<AminoAcid&> (*components[n])[C]);
    // We set the internal numbering of each of the new amminoacids.
    unsigned int index = dynamic_cast<AminoAcid&> (*components[n])[O].getNumber() + dynamic_cast<AminoAcid&> (*components[n]).getSideChain().size();
    (*aa)[N].setNumber(index + 1);
    (*aa)[CA].setNumber(index + 2);
    (*aa)[C].setNumber(index + 3);
    (*aa)[O].setNumber(index + 4);
    //	//  1. TRANSLATION to origin of N-atom

    vgVector3<double> orig = dynamic_cast<AminoAcid&> (*components[n])[N].getCoords();

    aa->addTrans(orig);
    aa->sync();
    //	//  2. FIRST ROTATION	
    //	// vector that identifies final aminoacid orientation
    vgVector3<double> v12 = ca2 - ca1;
    //	// vector that identifies initial aminoacid orientation
    vgVector3<double> N_C = (*aa)[C].getCoords() - (*aa)[N].getCoords();
    //	// CA-N vector
    vgVector3<double> CA_N = (*aa)[N].getCoords() - (*aa)[CA].getCoords();
    //	// CA-C vector
    vgVector3<double> CA_C = (*aa)[C].getCoords() - (*aa)[CA].getCoords();
    //	//normal to plane identified by CA-N and CA-C
    vgVector3<double> CAN_CAC_ax = (CA_N.normalize()).cross(CA_C.normalize());
    //	//angle between N_C and v12 
    double NC_v12_ang = icc.getAngle(N_C.normalize(), v12.normalize());
    //	//rotate
    aa->addRot(vgMatrix3<double>::createRotationMatrix(CAN_CAC_ax, -NC_v12_ang + 0));
    aa->sync();
    //	//  3.  SECOND ROTATION 	
    //	//aminoacid orientation updating
    N_C = (*aa)[C].getCoords() - (*aa)[N].getCoords();
    //	//new angle between NC and v12
    NC_v12_ang = icc.getAngle(N_C.normalize(), v12.normalize());
    //	//rotation axis: normal to plane identified by v12 and current N_C
    vgVector3<double> NC_v12_ax = (N_C.normalize()).cross(v12.normalize());
    //	//rotate
    aa->addRot(vgMatrix3<double>::createRotationMatrix(NC_v12_ax, -NC_v12_ang + 0));
    aa->sync();
    //	//  4.  TRANSLATION back to final position
    //	// Translate Calpha to final position
    vgVector3<double> caTransBack = ca1 - (*aa)[CA].getCoords();
    //	//translate
    aa->addTrans(caTransBack);
    aa->sync();
    components.push_back(aa);
}

/**
 *@Description Insert a new aminoacid of type 'code' after position 'n'.
 *                 Atoms are set according to lenght and angle cristallographic
 *                 values. Torsion angles are given in input (DEFAULT: -62�, -41�)
 *  WARNING: This function is only implemented for 'flat' spacers. 
 *           Ie. no subspacers can be handled.
 *  NB: This function is invoked by the Nazgul module to temporarily (!!!) 
 *      fill gaps.
 *@param string code, unsigned int n, double phi, double psi, double omega, 
                         double NToCaLen, double CaToCLen, double CToOLen, double atCaToCAng, double CaToOAng, 
                         double CtoNLen, double atCToNAng, double OToNAng, double atNToCaAng
 *@return  void
 */
void Spacer::insertAminoAfter(string code, unsigned int n, double phi, double psi, double omega,
        double NToCaLen, double CaToCLen, double CToOLen, double atCaToCAng, double CaToOAng,
        double CtoNLen, double atCToNAng, double OToNAng, double atNToCaAng) {

    if (n == 9999) {
        n = sizeAmino() - 1;
    }
    if (!isGap(getPdbNumberFromIndex(n + 1))) {
        ERROR("Spacer::InsertAminoAfter() is designed to insert aminoacids to fill gaps or at the end of spacer",
                exception);
    }
    // cout <<"====<" <<endl;

    if (subSpacerList.size() > 0)
        ERROR("Spacer::insertAminoAfter() is only implemented for 'flat' spacers.",
            exception);

    if (n > components.size())
        ERROR("Spacer::insertAminoAfter() : Argument out of scope.", exception);

    IntCoordConverter icc;
    AminoAcid* aa = new AminoAcid;

    aa->setType(code);
    aa->getSideChain().setType(code);
    aa->getSideChain().setBackboneRef(aa);

    // index of last atom of previous aminoacid

    unsigned int index = dynamic_cast<AminoAcid&> (*components[n])[O].getNumber()
            + dynamic_cast<AminoAcid&> (*components[n]).getSideChain().size();

    Atom at, at2, at3, at4;

    // set N atom
    at.setType("N");
    at.setNumber(index + 1);
    aa->addAtom(at);

    // bind N-term to C-term of previous aminoacid
    aa->bindIn((*aa)[N], dynamic_cast<AminoAcid&> (*components[n]),
            dynamic_cast<AminoAcid&> (*components[n])[C]);

    Atom &ref1 = dynamic_cast<AminoAcid&> (*components[n])[C];
    Atom &ref2 = dynamic_cast<AminoAcid&> (*components[n])[CA];
    Atom &ref3 = dynamic_cast<AminoAcid&> (*components[n])[O];
    Atom &ref4 = (*aa)[N];
    icc.zAtomToCartesian(ref1,
            CtoNLen,
            ref2,
            atCToNAng,
            ref3,
            OToNAng, 1, ref4);

    // fix position relative to previous:
    (*aa)[N].setTrans((*aa)[N].getCoords()
            - dynamic_cast<AminoAcid&> (*components[n])[C].getCoords());

    // set CA atom
    at2.setType("CA");
    at2.setNumber(index + 2);
    aa->addAtom(at2);

    (*aa)[N].bindOut((*aa)[CA]);

    icc.zAtomToCartesian(dynamic_cast<AminoAcid&> (*aa)[N],
            NToCaLen,
            dynamic_cast<AminoAcid&> (*components[n])[C],
            atNToCaAng,
            dynamic_cast<AminoAcid&> (*components[n])[CA], omega,
            0, (*aa)[CA]);

    // set C atom
    at3.setType("C");
    at3.setNumber(index + 3);
    aa->addAtom(at3);

    (*aa)[CA].bindOut((*aa)[C]);

    icc.zAtomToCartesian((*aa)[CA], CaToCLen,
            (*aa)[N], atCaToCAng,
            dynamic_cast<AminoAcid&> (*components[n])[C], phi,
            0, (*aa)[C]);

    // set O atom
    at4.setType("O");
    at4.setNumber(index + 4);
    aa->addAtom(at4);

    (*aa)[C].bindOut((*aa)[O]);

    icc.zAtomToCartesian((*aa)[C], CToOLen,
            (*aa)[CA], CaToOAng,
            (*aa)[N], psi, 0, (*aa)[O]);

    icc.setTorsionAngle(dynamic_cast<AminoAcid&> (*components[n])[N],
            dynamic_cast<AminoAcid&> (*components[n])[CA],
            dynamic_cast<AminoAcid&> (*components[n])[C],
            (*aa)[N], DEG2RAD * psi);
    aa->patchBetaPosition(index + 5);

    int nOldSizeAmino = sizeAmino();
    if (((int) n < nOldSizeAmino - 1) && (dynamic_cast<AminoAcid&> (*components[n]).sizeOutBonds() > 0)) {
        cout << "Hello world" << endl;
        dynamic_cast<AminoAcid&> (*components[n]).unbindOut(
                dynamic_cast<AminoAcid&> (*components[n]).getOutBond(0));
    }

    static int XXX = 0;
    components.insert(components.begin() + n + 1, aa);
    if ((nOldSizeAmino >= 2) && ((int) n < nOldSizeAmino - 2)) {
        ++XXX;
        cout << "((" << XXX << "))";
        aa->bindOut((*aa)[C], dynamic_cast<AminoAcid&> (*components[n + 2]),
                dynamic_cast<AminoAcid&> (*components[n + 2])[N]);
    }

    aa->sync();
}

/**
 *@Description  This function is only implemented for 'flat' spacers. 
 *           Ie. no subspacers can be handled.
 *  NB: This function is invoked by the Nazgul module to temporarily (!!!) 
 *      fill gaps.
 *  
 *  WARNING: this method still requires to be corrected: problems are found
 *           in positioning the O atom.
 *@param string code, unsigned int p, double phi, double psi, double omega,
                          double NToCaLen, double CaToCLen, double CToOLen, double atCaToCAng, double CaToOAng, 
                          double CtoNLen, double atCToNAng, double OToNAng, double atNToCaAng
 *@return  void
 */
void Spacer::insertAminoBefore(string code, unsigned int p, double phi, double psi, double omega,
        double NToCaLen, double CaToCLen, double CToOLen, double atCaToCAng, double CaToOAng,
        double CtoNLen, double atCToNAng, double OToNAng, double atNToCaAng) {

    if (p == 0)
        addGap(getPdbNumberFromIndex(p - 1));

    //this method modifies startOffset if it is > 0
    if (!isGap(getPdbNumberFromIndex(p - 1))) {
        ERROR("Spacer::InsertAminoBefore() is designed to insert aminoacids to fill gaps or at the beginning of spacer",
                exception);
    }

    if (getPdbNumberFromIndex(p - 1) <= startOffset)
        startOffset--;
    else
        removeGap(getPdbNumberFromIndex(p - 1));

    for (unsigned int a = 0; a < sizeAmino(); a++)
        for (unsigned int atm = 0; atm < dynamic_cast<AminoAcid&> (*components[a]).size(); atm++)
            dynamic_cast<AminoAcid&> (*components[a])[atm].
                setNumber(dynamic_cast<AminoAcid&> (*components[a])[atm].getNumber() + 4);


    if (subSpacerList.size() > 0)
        ERROR("Spacer::insertAminoAfter() is only implemented for 'flat' spacers.",
            exception);

    if (p > components.size())
        ERROR("Spacer::insertAminoAfter() : Argument out of scope.", exception);

    IntCoordConverter icc;
    AminoAcid* aa = new AminoAcid;

    aa->setType(code);
    aa->getSideChain().setType(code);
    aa->getSideChain().setBackboneRef(aa);

    // index of first atom of previous aminoacid
    unsigned int index = dynamic_cast<AminoAcid&> (*components[p])[N].getNumber();

    Atom at, at2, at3, at4;

    // set C atom
    at.setType("C");
    at.setNumber(index - 1);
    aa->addAtom(at);

    // bind C term to N term of subsequent aminoacid
    aa->bindOut((*aa)[C], (*aa), dynamic_cast<AminoAcid&> (*components[p])[N]);

    dynamic_cast<AminoAcid&> (*components[p])[N].bindIn((*aa)[C]);

    // set position according to existent atoms using phi torsion angle
    icc.zAtomToCartesian(dynamic_cast<AminoAcid&> (*components[p])[N],
            /*BOND_LENGTH_CPRIME_TO_N*/CtoNLen,
            dynamic_cast<AminoAcid&> (*components[p])[CA],
            /*BOND_ANGLE_AT_N_TO_CALPHA*/atNToCaAng,
            dynamic_cast<AminoAcid&> (*components[p])[C],
            phi, 0, (*aa)[C]);
    aa->sync();

    // fix position relative to next
    dynamic_cast<AminoAcid&> (*components[p])[N].setTrans(
            dynamic_cast<AminoAcid&> (*components[p])[N].getCoords()
            - (*aa)[C].getCoords());

    // set CA atom
    at3.setType("CA");
    at3.setNumber(index - 2);
    aa->addAtom(at3);

    (*aa)[C].bindIn((*aa)[CA]);
    // set position according to existent atoms using omega torsion angle 
    icc.zAtomToCartesian((*aa)[C],
            CaToCLen,
            dynamic_cast<AminoAcid&> (*components[p])[N],
            atCToNAng,
            dynamic_cast<AminoAcid&> (*components[p])[CA],
            omega, 0, (*aa)[CA]);
    aa->sync();

    // set N atom
    at4.setType("N");
    at4.setNumber(index - 4);

    aa->addAtom(at4);

    (*aa)[CA].bindIn((*aa)[N]);
    // set position according to existent atoms using psi torsion angle
    icc.zAtomToCartesian((*aa)[CA], NToCaLen,
            (*aa)[C], atCaToCAng,
            dynamic_cast<AminoAcid&> (*components[p])[N],
            psi, 0, (*aa)[N]);
    aa->sync();


    at2.setType("O");
    at2.setNumber(index - 3);

    aa->addAtom(at2);

    (*aa)[C].bindOut((*aa)[O]);
    icc.zAtomToCartesian((*aa)[C], CToOLen, (*aa)[CA],
            CaToOAng, dynamic_cast<AminoAcid&> (*components[p])[N],
            OToNAng, 1, (*aa)[O]);

    aa->sync();

    aa->patchBetaPosition();

    components.insert(components.begin() + p, aa);
}



/**
 *@Description synchronize coords with structure
 *@param none
 *@return  changes are made internally(void)
 */
// 

void Spacer::sync() {
    if (!modified)
        return;

    for (unsigned int i = 0; i < sizeAmino(); i++)
        getAmino(i).sync();

    modified = false;

    resetBoundaries();
}

/**
 *@Description Reset the amino acids boundaries
 *@param none
 *@return changes are made internally(void)
 */
void Spacer::resetBoundaries() {
    vgVector3<double> tmpV(DBL_MAX - 1, DBL_MAX - 1, DBL_MAX - 1);
    lowerBound = tmpV;
    upperBound = -tmpV;
    unsigned int sizeA = sizeAmino();
    for (unsigned int i = 0; i < sizeA; i++)
        for (unsigned int j = 0; j < 3; j++) {
            if (getAmino(i).getLowerBound()[j] < lowerBound[j])
                lowerBound[j] = getAmino(i).getLowerBound()[j];
            if (getAmino(i).getUpperBound()[j] > upperBound[j])
                upperBound[j] = getAmino(i).getUpperBound()[j];
        }
}

/**
 *@Description Returns the first and last aminoacid index from the sub-spacer in a <pair>.
 *                 The difference (last - first) is the count of amino acids
 *                 in the whole sub-spacer structure.
 *@param corresponding index(unsigned int) 
 *@return the start and end indexes in the range(pair<unsigned int, unsigned int> )
 */

pair<unsigned int, unsigned int> Spacer::getSubSpacerListEntry(unsigned int n) {
    if (n >= subSpacerList.size())
        ERROR("Index is out of range", exception);
    return subSpacerList[n];
}

/**
 *@Description Returns the first and last aminoacid index from the sub-spacer in a <pair>.
 *                 The difference (last - first) is the count of amino acids
 *                 in the whole sub-spacer structure.
 *@param corresponding index(unsigned int) 
 *@return the start and end indexes in the range(pair<unsigned int, unsigned int> )
 */
pair<unsigned int, unsigned int> Spacer::getSubSpacerListEntry(unsigned int n) const {
    if (n >= subSpacerList.size())
        ERROR("Index is out of range", exception);
    return subSpacerList[n];
}

/**
 *@Description copies the given the spacer
 *@param reference to the spacer to copy(spacer&)
 *@return  changes are made internally(void)
 */
void Spacer::copy(const Spacer& orig) {
    PRINT_NAME;

    startOffset = orig.startOffset;
    startAtomOffset = orig.startAtomOffset;
    gaps = orig.gaps;
    subSpacerList.clear();

    for (unsigned int n = 0; n < orig.sizeSpacer(); n++)
        subSpacerList.push_back(orig.getSubSpacerListEntry(n));
    Polymer::copy(orig);

    if ((sizeAmino() > 1) && (orig.getAmino(0)[C].isBond(orig.getAmino(1)[N])))
        for (unsigned int i = 0; i < sizeAmino() - 1; i++) {
            if (getAmino(i + 1).getCode() == PRO) { // fix proline ring:
                if (getAmino(i + 1).isMember(CD))
                    getAmino(i + 1)[CD].unbindOut(getAmino(i + 1)[CA]);
                getAmino(i).bindOut(getAmino(i)[C], getAmino(i + 1),
                        getAmino(i + 1)[N]);
                if (getAmino(i + 1).isMember(CD))
                    getAmino(i + 1)[CD].bindOut(getAmino(i + 1)[CA]);
            } else
                getAmino(i).bindOut(getAmino(i)[C], getAmino(i + 1), getAmino(i + 1)[N]);
        }
}

/**
 *@Description clone the spacer
 *@param void
 *@return  pointer to the new component  (Component*)
 */
Component* Spacer::clone() {
    Spacer* tmp = new Spacer;
    tmp->copy(*this);
    return tmp;
}

// OPERATORS:

/**
 *@Description Assign the spacer to another spacer.
 *@param reference for the original spacer(const Spacer&)
 *@return  reference to the new spacer(Spacer&)
 */
Spacer& Spacer::operator=(const Spacer& orig) {
    PRINT_NAME;
    if (&orig != this)
        copy(orig);
    return *this;
}


// HELPERS:

/**
 *@Description Insert a component(only an Amino acid or Sapacer class type are allowed) 
 * in the spacer at the back side.  It doesn't connect the components.
 *@param pointer to the component to insert
 *@return  changes are made internally(void)
 */

void Spacer::insertComponent(Component* c) {
    pair<unsigned int, unsigned int> tmp;
    unsigned int count = 0;
    unsigned int count_this = sizeAmino();
    if (c->hasSuperior())
        ERROR("Component does have a superior", exception);
    if (c->getClassName() == "AminoAcid") {
        count = 1;
    } else {
        if (c->getClassName() == "Spacer") {
            count = (dynamic_cast<Spacer*> (c))->sizeAmino();
            tmp.first = count_this;
            tmp.second = count_this + count - 1;
            subSpacerList.push_back(tmp);
        } else
            ERROR("The component is neither an aminoacid nor a spacer.", exception);
    }
    Polymer::insertComponent(c);
    setModified();
    if (hasSuperior())
        (dynamic_cast<Spacer&> (getSuperior()))
        .modifySubSpacerList(this, static_cast<int> (count));
}

/**
 *@Description If something has changed in the spacer structure,
 *                 this method updates the spacer given as a parameter .
 * 
 *@param pointer to the spacer to modify(Spacer*),  number of aminoacids which has inserted
//                 (positive value) or deleted (negativ value)(int).  
 *@return  all changes are made internally(void)
 */

void Spacer::modifySubSpacerList(Spacer* s, int count) {
    unsigned int index = 0;
    while (&getSpacer(index) != s)
        index++;
    subSpacerList[index].second = subSpacerList[index].second + count;
    index++;
    for (unsigned int loop = index; loop < subSpacerList.size(); loop++) {
        subSpacerList[loop].first = subSpacerList[loop].first + count;
        subSpacerList[loop].second = subSpacerList[loop].second + count;
    }
    if (hasSuperior()) {
        (dynamic_cast<Spacer&> (getSuperior()))
                .modifySubSpacerList(this, static_cast<int> (count));
    } else {
        return;
    }
}

/**
 *@Description Returns the spacer of the index 
 *@param unsigned int 
 *@return  spacer reference
 */

Spacer& Spacer::getSpacer(unsigned int n) {
    PRECOND(n < subSpacerList.size(), exception);
    unsigned int idxOffset = 0;
    // count of aminoacids before the first spacer 
    idxOffset = subSpacerList[0].first;
    // count all aminoacids between the spacers
    for (unsigned int loop = 1; loop < n; loop++)
        idxOffset += subSpacerList[loop].first - subSpacerList[loop - 1].second;
    if (components[idxOffset]->getClassName() == "Spacer") {
        return (dynamic_cast<Spacer&> (*components[idxOffset]));
    } else
        ERROR("Component is not a spacer", exception);
}

/**
 *@Description Returns the spacer of the corresponding index 
 *@param unsigned int 
 *@return  spacer reference
 */

const Spacer& Spacer::getSpacer(unsigned int n) const {
    PRECOND(n < subSpacerList.size(), exception);
    unsigned int idxOffset = 0;
    // count of aminoacids before the first spacer
    idxOffset = subSpacerList[0].first;
    // count all aminoacids between the spacers
    for (unsigned int loop = 1; loop < n; loop++)
        idxOffset += subSpacerList[loop].first - subSpacerList[loop - 1].second;
    if (components[idxOffset]->getClassName() == "Spacer") {
        return (dynamic_cast<Spacer&> (*components[idxOffset]));
    } else
        ERROR("Component is not a spacer", exception);

}

/**
 *@Description Returns the aminoacid of the index 
 *@param unsigned int 
 *@return aminoacid reference
 */

AminoAcid& Spacer::getAmino(unsigned int n) { // needs to be recoded & optimized
    PRECOND(n < this->sizeAmino(), exception);
    unsigned int idxOffset = 0;
    unsigned int lastMax = 0;
    unsigned int subSpacerSize = subSpacerList.size();
    unsigned int betweenOffset = 0;
    for (unsigned int loop = 0; loop < subSpacerSize; loop++) {
        betweenOffset += subSpacerList[loop].first - lastMax;
        // count aminoacids between the spacers and the spacer itself
        if ((n >= subSpacerList[loop].first) &&
                (n <= subSpacerList[loop].second)) {
            return ( dynamic_cast<AminoAcid&>
                    (dynamic_cast<Spacer*>
                    (components[betweenOffset])
                    ->getAmino(n - subSpacerList[loop].first)));
        }
        if ((n >= lastMax) && (n < subSpacerList[loop].first)) {
            if ((*components[n - idxOffset])
                    .getClassName() != "AminoAcid")
                ERROR("Component is not an aminoacid !", exception);
            return (dynamic_cast<AminoAcid&>
                    (*components[n - idxOffset])); // return aminoacid
        }
        lastMax = subSpacerList[loop].second;
        idxOffset += lastMax - subSpacerList[loop].first;
    }
    if ((*components[n - idxOffset])
            .getClassName() != "AminoAcid")
        ERROR("The component is not an aminoacid !", exception);
    return (dynamic_cast<AminoAcid&>
            (*components[n - idxOffset])); // return aminoacid
}

/**
 *@Description Returns the aminoacid of index
 *@param unsigned int 
 *@return  aminoacid reference
 */
const AminoAcid& Spacer::getAmino(unsigned int n) const { // needs to be recoded & optimized

    PRECOND(n < this->sizeAmino(), exception);
    unsigned int idxOffset = 0;
    unsigned int lastMax = 0;
    unsigned int subSpacerSize = subSpacerList.size();
    unsigned int betweenOffset = 0;
    for (unsigned int loop = 0; loop < subSpacerSize; loop++) {
        betweenOffset += subSpacerList[loop].first - lastMax;
        // count aminoacids between the spacers and the spacer itself
        if ((n >= subSpacerList[loop].first) &&
                (n <= subSpacerList[loop].second)) {
            return ( dynamic_cast<AminoAcid&>
                    (dynamic_cast<Spacer*>
                    (components[betweenOffset])
                    ->getAmino(n - subSpacerList[loop].first)));
        }
        if ((n >= lastMax) && (n < subSpacerList[loop].first)) {
            if ((*components[n - idxOffset])
                    .getClassName() != "AminoAcid")
                ERROR("Component is not an aminoacid !", exception);
            return (dynamic_cast<AminoAcid&>
                    (*components[n - idxOffset])); // return aminoacid
        }
        lastMax = subSpacerList[loop].second;
        idxOffset += lastMax - subSpacerList[loop].first;
    }
    if ((*components[n - idxOffset])
            .getClassName() != "AminoAcid")
        ERROR("The component is not an aminoacid !", exception);
    return (dynamic_cast<AminoAcid&>
            (*components[n - idxOffset])); // return aminoacid
}

/**
 *@Description Merge the components of the sub-spacer <s> at the place of <s> in the spacer <this>..
 *@param spacer reference 
 *@return  void
 */

void Spacer::mergeSpacer(Spacer* s) {
    PRECOND(&(s->getSuperior()) == this, exception);
    unsigned int idx = 0;
    vector<Component*> tmpComponents;
    tmpComponents.clear();
    for (unsigned int n = 0; n < components.size(); n++)
        if (s == components[n]) {
            idx = n;
            n = components.size();
        }
    for (unsigned int n = 0; n < idx; n++)
        tmpComponents.push_back(components[n]);
    for (unsigned int n = 0; n < s->components.size(); n++) {
        tmpComponents.push_back(s->components[n]);
        if (tmpComponents[tmpComponents.size() - 1]->getClassName() == "Spacer")
            tmpComponents[tmpComponents.size() - 1]->setSuperior(this);
    }
    for (unsigned int n = idx + 1; n < components.size(); n++)
        tmpComponents.push_back(components[n]);
    components.clear();
    for (unsigned int n = 0; n < tmpComponents.size(); n++)
        components.push_back(tmpComponents[n]);
    updateSubSpacerList();
    s->components.clear();
    s->setSuperior(static_cast<Component*> (NULL));
    deleteComponent(s);
}

/**
 *@Description  Create a new sub-spacer at the place of <start> and put
 *                 all components between <start> and <end> in the new 
 *                 sub-spacer. The method returns the pointer at the 
 *                 new sub-spacer.
 *@param unsigned int start, unsigned int end 
 *@return  spacer reference
 */

Spacer* Spacer::splitSpacer(unsigned int start, unsigned int end) {
    PRECOND((start < components.size()) && (end < components.size()), exception);
    if (start > end) { // swap
        unsigned int tmp = end;
        end = start;
        start = tmp;
    }
    Spacer* newSpacer = new Spacer;
    // connect components to the new spacer
    for (unsigned int n = 0; n <= end - start; n++) {
        newSpacer->components.push_back(components[start + n]);
        newSpacer->components[n]->setSuperior(newSpacer);
    }
    newSpacer->updateSubSpacerList();
    unsigned int m = 0;
    for (unsigned int n = end + 1; n < components.size(); n++)
        components[start + 1 + m++] = components[n];
    newSpacer->setSuperior(this);
    components[start] = newSpacer; // insert the new spacer
    for (unsigned int n = 0; n < end - start; n++)
        components.pop_back();
    updateSubSpacerList();
    return newSpacer;
}

/**
 *@Description insert a new sub spacer inside an empty spacer structure
 *                 WARNING: this method works only with superior spacers
 *                 which are made of other subspacers only (child components
 *                 must be all spacers, the presence of aminoacids is not 
 *                 considered)
 *                 (This assumption is sufficient for the Thor module)
 *@param spacer reference 
 *@return  void
 */

///insert a new sub spacer inside an empty spacer structure
///ATTENTION: this method works only with superior spacers which are made of other subspacers 
///only (child components must be all spacers, the presence of aminoacids is not considered)
///(This assumption is sufficient for the Thor module)

void Spacer::insertFirstSpacer(Spacer* sp) {
    if (sizeSpacer() != 0 || sizeAmino() != 0)
        ERROR("Spacer::insertFirstSpacer() can only be invoked on an empty spacer", exception);

    components.push_back(sp);
    sp->setSuperior(this);
    updateSubSpacerList();
}

/**
 *@Description insert a new sub spacer inside a spacer structure.
 *                  N.B. spacers are not positioned, they maintain original coordinates
 *                 WARNING: this method works only with superior spacers
 *                 which are made of other subspacers only (child components
 *                 must be all spacers, the presence of aminoacids is not 
 *                 considered)
 *                 (This assumption is sufficient for the Thor module)
 *@param Spacer reference, unsigned int 
 *@return void
 */
///insert a new sub spacer inside a spacer structure.
///ATTENTION: this method works only with superior spacers which are made of other subspacers 
///only (child components must be all spacers, the presence of single aminoacids is not considered)
///(This assumption is sufficient for the Thor module)

void Spacer::insertSubSpacerAfter(Spacer* sp, unsigned int pos) {
    int lastPdbPred = getSpacer(pos).getPdbNumberFromIndex(getSpacer(pos).sizeAmino() - 1);

    //if inserted spacer is not the last one...
    if (pos != sizeSpacer() - 1) {
        int firstPdbSucc =
                getSpacer(pos + 1).getPdbNumberFromIndex(getSpacer(pos + 1).sizeAmino() - 1);

        if (sp->getPdbNumberFromIndex(0) <= lastPdbPred ||
                sp->getPdbNumberFromIndex(sp->sizeAmino() - 1) >= firstPdbSucc)
            ERROR("Spacer::insertSubSpacerAfter: attempt to insert a new spacer whith aminoacid superposition", exception);

        //if a spacer exists before and after the one to be inserted, gaps are present between spacer 
        //in position pos and pos+1: they have to be removed
        for (unsigned int i = 0; i < sp->sizeAmino(); i++)
            removeGap(sp->getPdbNumberFromIndex(i));

        //The following two controls are maybe unnecessary

        //if first pdb number of inserted spacer and the last one of its predecessor
        //are not contiguous, gaps have to be inserted
        if (lastPdbPred + 1 != sp->getPdbNumberFromIndex(0))
            for (int i = lastPdbPred + 1; i < sp->getPdbNumberFromIndex(0); i++)
                addGap(i);

        else {
            //N term to C term of previous aminoacid

            (sp->getAmino(0)[N]).bindIn(getSpacer(pos).getAmino(getSpacer(pos).sizeAmino() - 1)[C]);

            //bind aminoacids
            (sp->getAmino(0)).bindIn(sp->getAmino(0)[N],
                    getSpacer(pos).getAmino(getSpacer(pos).sizeAmino() - 1),
                    getSpacer(pos).getAmino(getSpacer(pos).sizeAmino() - 1)[C]);
        }

        //if last pdb number of inserted spacer and the first one of its successor
        //are not contiguous, gaps have to be inserted
        if (firstPdbSucc != sp->getPdbNumberFromIndex(sp->sizeAmino() - 1) + 1)
            for (int i = sp->getPdbNumberFromIndex(sp->sizeAmino() - 1) + 1; i < firstPdbSucc; i++)
                addGap(i);

        else {
            // bind C term to N term of subsequent aminoacid
            (sp->getAmino(sp->sizeAmino() - 1)[C]).bindOut(getSpacer(pos).getAmino(0)[C]);

            //bind aminoacids
            (sp->getAmino(sp->sizeAmino() - 1)).bindOut(sp->getAmino(sp->sizeAmino() - 1)[C],
                    getSpacer(pos).getAmino(0),
                    getSpacer(pos).getAmino(0)[N]);
        }
    } else {
        //spacer inserted is last one

        //if first pdb number of inserted spacer and the last one of its predecessor
        //are not contiguous, gaps have to be inserted
        if (lastPdbPred + 1 != sp->getPdbNumberFromIndex(0))
            for (int i = lastPdbPred + 1; i < sp->getPdbNumberFromIndex(0); i++)
                addGap(i);
        else {
            (sp->getAmino(0)[N]).bindIn(getSpacer(pos).getAmino(getSpacer(pos).sizeAmino() - 1)[C]);

            //bind aminoacids
            (sp->getAmino(0)).bindIn(sp->getAmino(0)[N],
                    getSpacer(pos).getAmino(getSpacer(pos).sizeAmino() - 1),
                    getSpacer(pos).getAmino(getSpacer(pos).sizeAmino() - 1)[C]);
        }
    }

    components.insert(components.begin() + pos + 1, sp);
    sp->setSuperior(this);
    updateSubSpacerList();
}

/**
 *@Description Rearrange the spacer (and all sub-spacers) as a flat tree. After this there are no more sub-spacers..
 *@param unsigned int 
 *@return  int
 */

void Spacer::makeFlat() {
    while (sizeSpacer() > 0)
        mergeSpacer(&getSpacer(0));
}

/**
 *@Description Group the spacer (and his sub-spacers) like the state codes
 *                 of the aminoacids. After this the aminoacid structure is
 *                 grouped in HELIX-spacer,STRANDS-spacer,COIL-spacer and 
 *                 TURN-spacer.
 *@param void
 *@return  void
 */

void Spacer::groupLikeStateCode() {
    if (sizeAmino() < 3)
        return;
    unsigned int start = 0;
    unsigned int end = 0;
    StateCode prev = getAmino(0).getState();
    makeFlat();
    for (unsigned int n = 1; n < sizeAmino() - 1; n++) {
        if (getAmino(n).getState() == prev) {
            end++;
        } else {
            if (((end - start) > 0)
                    && ((prev == StateCode(HELIX)) || (prev == StateCode(STRAND)
                    || (prev == StateCode(COIL)) || (prev == StateCode(TURN)))))
                splitSpacer(start, end);
            end -= (end - start) - 1;
            start = end;
            prev = getAmino(n).getState();
        }
        if (getAmino(n).getState() == prev)
            if (((end - start) > 0)
                    && ((prev == StateCode(HELIX)) || (prev == StateCode(STRAND)
                    || (prev == StateCode(COIL)) || (prev == StateCode(TURN))))) {
                cout << "Last splitSpacer " << start << "  " << end + 1 << endl;
                splitSpacer(start, ++end);
            }
    }
}

/**
 *@Description Update the internal structure of a spacer after a change
 *                 in the structure of a spacer.
 *@param void 
 *@return  void
 */

void Spacer::updateSubSpacerList() {
    subSpacerList.clear();
    pair < unsigned int, unsigned int > tmp;
    unsigned int countAA = 0;
    for (unsigned int n = 0; n < components.size(); n++) {
        if (components[n]->getClassName() == "Spacer") {
            if (sizeAmino() == 0)
                ERROR("An empty spacer was found as sub-spacer", exception);
            tmp.first = countAA;
            countAA += (static_cast<Spacer*> (components[n]))->sizeAmino();
            tmp.second = countAA - 1;
            subSpacerList.push_back(tmp);
        } else {
            countAA++;
        }
    }
}

/**
 *@Description Set the state (HELIX or STRAND ...) from the string sec
 *                 which may be generated for example by DSSP.
 *@param string sequence 
 *@return  void
 */
// sets states from secondary string 

void Spacer::setStateFromSecondary(string sec) {
    unsigned int len = sec.length();

    cout << "len= " << sec.length() << "\t sizeAmino= " << sizeAmino() << "\n";

    if (len < sizeAmino())
        cout << "Warning: Spacer::setStateFromSecondary() : string is too short.\n";
    if (len > sizeAmino()) {
        cout << "Warning: Spacer::setStateFromSecondary()"
                << " : string is too long.\n";
        len = sizeAmino();
    }

    for (unsigned int n = 0; n < len; n++) {
        StateCode tmp;
        if (sec[n] == 'H')
            tmp = HELIX;
        else if (sec[n] == 'E')
            tmp = STRAND;
        else
            tmp = COIL;
        getAmino(n).setState(tmp);
    }
}




/**
 *@Description Set the state (HELIX or STRAND ...) from the torsion angles of the aminoacids.
 *@param void
 *@return  void
 */
// sets states from its phi/psi angles

void Spacer::setStateFromTorsionAngles() {
    // NB: This definition was taken from McGuffin et al.,
    // Bioinformatics (17):63-72 (2001)

    for (unsigned int n = 0; n < sizeAmino(); n++)
        getAmino(n).setStateFromTorsionAngles();

    unsigned int start = 0;
    StateCode state = COIL;

    for (unsigned int i = 0; i < sizeAmino(); i++) {
        if (getAmino(i).getState() != COIL) {
            if (state == COIL)
                start = i;
            else if (state != getAmino(i).getState()) {
                if ((((i - start) < 3) && (state == STRAND))
                        || (((i - start) < 4) && (state == HELIX)))
                    for (unsigned int j = start; j < i; j++)
                        getAmino(j).setState(COIL); // too short, remove states
                start = i;
            }
        } else {
            if ((((i - start) < 3) && (state == STRAND))
                    || (((i - start) < 4) && (state == HELIX)))
                for (unsigned int j = start; j < i; j++)
                    getAmino(j).setState(COIL); // too short, remove states
        }
        state = getAmino(i).getState();
    }

    if (sizeAmino() > 0) {
        getAmino(0).setState(COIL); // first & last residue are
        getAmino(sizeAmino() - 1).setState(COIL); // *always* COIL
    }
}

/**
 *@Description Prints the subspacer list.
 *@param void
 *@return  void
 */
void Spacer::printSubSpacerList() {
    cout << "-------------------------------------------------------" << endl;
    cout << "***************** THE SUBSPACERLIST *******************" << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << "SIZE AMINOACIDS.:" << sizeAmino() << endl;
    cout << "SIZE SPACERS....:" << sizeSpacer() << endl;
    for (unsigned int n = 0; n < subSpacerList.size(); n++) {
        cout << "<" << subSpacerList[n].first << "," << subSpacerList[n].second << "> " << endl;
    }
}

/**
 *@Description Prints the components.
 *@param void
 *@return  void
 */
void Spacer::printComponents() {
    cout << "-------------------------------------------------------" << endl;
    cout << "***************** THE COMPONENTS ** *******************" << endl;
    cout << "-------------------------------------------------------" << endl;
    cout << "COUNT COMPONENTS :" << components.size() << endl;
    for (unsigned int n = 0; n < components.size(); n++) {
        if (components[n]->getClassName() == "Spacer") {
            cout << components[n]->getClassName() << " SIZE:" << (static_cast<Spacer*> (components[n]))->sizeAmino() << endl;
        } else {
            cout << components[n]->getClassName() << " adress:" << components[n] << endl;
        }
    }
}

/**
 *@Description This function returs a pair array that define holes into the PDB
 *       numbering.
 * NOTE: The function returns exactly where an hole starts and ends.
 *@param void
 *@return  vector of int pairs
 */
vector< pair<int, int> > Spacer::getHoles() {
    vector< pair<int, int> > arr;
    int prev = getPdbNumberFromIndex(0);
    // We test 1 and not 0 because a PDB file starts from 1 and not from 0.
    if (prev != 1)
        arr.push_back(pair<int, int>(0, prev));
    for (unsigned int i = 1; i < sizeAmino(); ++i) {
        int curr = getPdbNumberFromIndex(i);
        if (curr != prev + 1)
            arr.push_back(pair<int, int>(prev, curr));
        prev = curr;
    }

    return arr;
}

/**
 *@Description This helper function, is used to implement an optimization to mjollnir.
 *       What it does is simple and it works only for NMR PDB files. We load
 *       it and we read the number of different models stored. After it, we
 *       load every single model and once we did it, we calculate every single
 *       distance between Ca atoms. We return a position-related array, where
 *       each position is the minimum distance between each pair of Ca.
 *@param string &strFileName, vector<double> 
 *@return  bool
 */
bool Spacer::NMRGetMinimumCADistanceVector(const string &strFileName, vector<double> *pNMR) {
    if (!pNMR)
        return false;

    ifstream nmrFile(strFileName.c_str());
    if (!nmrFile) {
        cerr << "NMR file exception: does not exist" << endl;
        return false;
    }

    PdbLoader pdb(nmrFile);
    pdb.setNoHAtoms();
    vector<Spacer*> arr;
    // Remember that models go from [1..models] and not [0..models-1].
    // We load in N spacers all the models stored in the NMR file.
    // In this way we can make statistics over Ca distances.
    unsigned int models = pdb.getMaxModelsFast();
    for (unsigned int modelNumber = 1; modelNumber <= models; ++modelNumber) {
        pdb.setModel(modelNumber);
        Spacer *psp = new Spacer;
        if (!psp)
            return false;
        psp->load(pdb);
        arr.push_back(psp);
    }
    // First of all I read all the information about Ca distances,
    // but only if we were able to read at least one model from the NMR file.
    if (arr.size() == 0) {
        cerr << "No NMR information loaded" << endl;
        return false;
    }

    vector< vector< vgVector3<double> > > coords;
    unsigned int length = arr[0]->sizeAmino();
    for (unsigned int pos = 0; pos < length; ++pos) {
        vector< vgVector3<double> > tmp;
        for (vector<Spacer*>::iterator it = arr.begin(); it != arr.end(); ++it)
            tmp.push_back((*it)->getAmino(pos)[CA].getCoords());
        coords.push_back(tmp);
    }


    // We delete the memory used by the spacers.
    for (vector<Spacer*>::iterator it = arr.begin(); it != arr.end(); ++it)
        delete *it;
    // Now I calculate the minimum distance, for each position, between every
    // Ca pair.
    for (vector< vector< vgVector3<double> > >::iterator it = coords.begin(); it != coords.end(); ++it) {
        double currMinDistance = DBL_MAX;
        double currMaxDistance = DBL_MIN;
        double avgDistance = 0;
        for (vector< vgVector3<double> >::iterator prev = it->begin(); prev != it->end() - 1; ++prev)
            for (vector< vgVector3<double> >::iterator curr = prev + 1; curr != it->end(); ++curr) {
                vgVector3<double> diff = *prev - *curr;
                double currDistance = diff.length();
                if (currDistance < currMinDistance)
                    currMinDistance = currDistance;
                if (currDistance > currMaxDistance)
                    currMaxDistance = currDistance;
                avgDistance += currDistance;
            }
        cout << "Minimum distance=" << currMinDistance << endl;
        cout << "Maximum distance=" << currMaxDistance << endl;
        double nn = coords[0].size();
        cout << "Average distance=" << avgDistance / ((nn * (nn - 1)) / 2) << endl;
        // We save the minimum distance calculated for each position "it".
        pNMR->push_back(currMinDistance);
    }

    return true;
}

/**
 *@Description Set the state from H bonds.
 *@param void
 *@return  void
 * 
 * H bond definition: Kabsch, Sander, Biopolymers. 1983 Dec;22(12):2577-637.
 
 */

void Spacer::getBackboneHbonds() {
    IntCoordConverter icc;

    // Variables for H bonds
    vgVector3<double> bondVector; // H-N----O=C
    double distance;
    vgVector3<double> DCoords; // Donor  (N)
    vgVector3<double> ACoords; // Acceptor (O)
    vgVector3<double> HCoords; // Hydrogen
    vgVector3<double> HVector;
    double Hangle; // HNO angle (D-H-A)

    // Varibles for Bends
    vgVector3<double> CA1Coords;
    vgVector3<double> CA2Coords;
    vgVector3<double> preVector; //   CA(i) -->  CA(i-2) 
    vgVector3<double> postVector; // CA(i+2) -->  CA(i) 
    double Sangle;

    // Initialize H bonds matrix
    backboneHbonds = new bool*[sizeAmino()];
    for (unsigned int i = 0; i < sizeAmino(); i++) {
        backboneHbonds[i] = new bool[sizeAmino()];
        for (unsigned int j = 0; j < sizeAmino(); j++) {
            backboneHbonds[i][j] = false;
        }
    }

    // Initialize ss vector
    for (unsigned int i = 0; i < sizeAmino(); i++) {
        set<char> tmp;
        ss.push_back(tmp);
    }

    // Get H bonds for the backbone
    for (unsigned int i = 0; i < (sizeAmino()); i++) {
        for (unsigned int j = 0; j < sizeAmino(); j++) {
            if ((i != j) &&
                    (getAmino(j).isMember(N)) &&
                    (getAmino(j).isMember(H)) &&
                    (getAmino(i).isMember(O))) {

                ACoords = getAmino(i)[O].getCoords();
                HCoords = getAmino(j)[H].getCoords();
                DCoords = getAmino(j)[N].getCoords();

                HVector = getAmino(j)[H].getTrans();
                bondVector = icc.calculateTrans(DCoords, ACoords);
                Hangle = RAD2DEG * icc.getAngle(bondVector, HVector);

                distance = getAmino(j)[N].distance(getAmino(i)[O]);

                if ((Hangle <= 63.0) && (distance <= 5.2)) {
                    // It is important to preserve the direction
                    backboneHbonds[i][j] = true; // Rigth,   O ---> H-N
                }
            }
        }

        // Calculate bends
        if ((i > 1) && (i < (sizeAmino() - 2))) {
            CA1Coords = getAmino(i)[CA].getCoords();
            CA2Coords = getAmino(i - 2)[CA].getCoords();
            preVector = icc.calculateTrans(CA1Coords, CA2Coords);
            CA2Coords = getAmino(i + 2)[CA].getCoords();
            postVector = icc.calculateTrans(CA2Coords, CA1Coords);
            Sangle = RAD2DEG * icc.getAngle(preVector, postVector);
            if (Sangle > 70.0) {
                ss[i].insert('S');
            }
        }
    }
}


/**
 *@Description Set the state from H bonds.
 *@param void
 *@return  void
 * 
 * SS definition: Kabsch, Sander, Biopolymers. 1983 Dec;22(12):2577-637.
 * 
 * n-Turn: T = (i,i+n) n=3,4,5
 *         3 = (n = 3)
 *         4 = (n = 4)
 *         5 = (n = 5)
 *             
 * Bridge: B = (i-1,i,i+1) (j-1,j,j+1) nonoverlapping
 *         P = parallel (i-1,j) && (j,i+1) || (j-1,i) && (i,j+1)
 *         A =     anti (i,j) && (j,i) || (i-1,j+1) && (j-1,i+1)
 * 
 * Helices at least 2 consecutive n-Turns (longer helices by overlapping)
 *         G = 3-helix     : (i-1,i+2) && (i,i+3)
 *         H = 4-helix at i: (i-1,i+3) && (i,i+4)
 *         I = 5-helix     : (i-1,i+4) && (i,i+5)
 *  
 * Ladder: B = Consecutive bridges of identical type
 *  Sheet: E = One or more ladders connected by shared residues
 *   Bend: S = angle { CA(i)-CA(i-2), CA(i+2)-CA(i) } > 70
 * 
 * Structure Priority (switched E and B):
 * H E B G I T S
 * 
 * 
 * TODO: check for chain breaks.
 * TODO: remove small helices.
 * TODO: check for bulges.
 * 
 */
// calculate SS from H bonds

void Spacer::setDSSP(bool verbose) {

    // It also calculate bends
    getBackboneHbonds();

    // calculate n-Turns and bridges
    for (unsigned int i = 0; i < sizeAmino(); i++) {
        for (unsigned int j = 0; j < sizeAmino(); j++) {
            // Helices
            if (backboneHbonds[i][j]) {
                if (j == (i + 3)) {
                    for (unsigned int l = i + 1; l < j; l++) {
                        ss[l].insert('3');
                        ss[l].insert('T');
                    }
                }
                if (j == (i + 4)) {
                    for (unsigned int l = i + 1; l < j; l++) {
                        ss[l].insert('4');
                        ss[l].insert('T');
                    }
                }
                if (j == (i + 5)) {
                    for (unsigned int l = i + 1; l < j; l++) {
                        ss[l].insert('5');
                        ss[l].insert('T');
                    }
                }
            }
            // Bridges
            if (abs(i - j) > 2) { // Avoid overlapping stretches of 3 residues
                // Parallel bridges
                if ((i > 0) && (i < (sizeAmino() - 1))) {
                    if ((backboneHbonds[i - 1][j]) && (backboneHbonds[j][i + 1])) {
                        ss[i].insert('B');
                        ss[j].insert('B');
                    }
                }
                if ((j > 0) && (j < (sizeAmino() - 1))) {
                    if ((backboneHbonds[j - 1][i]) && (backboneHbonds[i][j + 1])) {
                        ss[i].insert('B');
                        ss[j].insert('B');
                    }
                }
                // Antiparallel bridges
                if ((backboneHbonds[i][j]) && (backboneHbonds[j][i])) {
                    ss[i].insert('B');
                    ss[j].insert('B');
                }
                if ((j > 0) && (j < (sizeAmino() - 1)) && (i > 0) && (i < (sizeAmino() - 1))) {
                    if ((backboneHbonds[i - 1][j + 1]) && (backboneHbonds[j - 1][i + 1])) {
                        ss[i].insert('B');
                        ss[j].insert('B');
                    }
                }
            }
        }
    }

    // Assign helices and sheets
    string turns = "435";
    string helices = "HGI";
    for (unsigned int i = 1; i < sizeAmino(); i++) {
        if (!ss[i].empty()) {
            // Set helices
            for (unsigned int l = 0; l < turns.size(); l++) {
                set<char> ::iterator it = ss[i].find(turns[l]);
                if (it != ss[i].end()) { // found a n-turn
                    char pos_s = *it;
                    int pos = atoi(&(pos_s));
                    set<char> ::iterator it1 = ss[i + pos - 1].find(turns[l]);
                    if (it1 != ss[i + pos - 1].end()) { // found the same n-turn after n-1 positions
                        bool helixBreak = false;
                        for (unsigned int m = i; m < (i + pos); m++) {
                            set<char> ::iterator it2 = ss[m].find(turns[l]);
                            if (it2 == ss[m].end()) {
                                helixBreak = true;
                            }
                        }
                        if (!helixBreak) {
                            for (unsigned int m = i; m < (i + pos); m++) {
                                ss[m].insert(helices[l]);
                            }
                        }
                        break;
                    }
                }
            }
            // Set sheets (E)
            set<char> ::iterator it = ss[i].find('B');
            if (it != ss[i].end()) {
                set<char> ::iterator it1 = ss[i - 1].find('B');
                if (it1 != ss[i - 1].end()) { // found a B
                    ss[i].insert('E');
                    ss[i - 1].insert('E');
                }
                if (i > 1) {
                    set<char> ::iterator it2 = ss[i - 2].find('B');
                    if (it2 != ss[i - 2].end()) { // found a B (breaks)
                        ss[i - 2].insert('E');
                        ss[i].insert('E');
                        ss[i - 1].insert('E');
                    }
                }
            }
        }
    }

    // Filter SS labelling by the "states" scheme
    string states = "HEBGITS";
    for (unsigned int i = 0; i < sizeAmino(); i++) {
        if (verbose)
            cout << getPdbNumberFromIndex(i) + 1 << "\t" << getAmino(i).getType1L() << "  ";
        if (!ss[i].empty()) {
            string states_tmp = " ";
            for (set<char> ::iterator it1 = ss[i].begin(); it1 != ss[i].end(); ++it1) {
                states_tmp += (*it1);
                states_tmp += ' ';
            }
            for (unsigned int j = 0; j < states.size(); j++) {
                set<char> ::iterator it;
                it = ss[i].find(states[j]);
                if (it != ss[i].end()) {
                    if (verbose)
                        cout << (*it) << "\t";
                    ss[i].clear();
                    ss[i].insert(states[j]);
                    break;
                }
            }
            if (verbose)
                cout << states_tmp;
        }
        if (verbose)
            cout << "\n";
    }
}


