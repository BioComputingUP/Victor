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
#include <IntLoader.h>
#include <IoTools.h>
#include <AminoAcid.h>
#include <IntCoordConverter.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Victor; using namespace Victor::Biopool;

/**
 *   Private helper function to set bond structure after loading the spacer.
 *@param  
 *@return  
 */
void
IntLoader::setBonds(Spacer& sp) {
    sp.getAmino(0).setBondsFromPdbCode(true);
    for (unsigned int i = 1; i < sp.size(); i++) {
        sp.getAmino(i).setBondsFromPdbCode(true, &(sp.getAmino(i - 1)));
    }
}

/**
 *   Converts from Tinker format into cartesian. Place the Atom "at" in the 3D space given a distance, an angle, the chirality and a torsion angle.
 * The input Atoms need to be connected in this order: atbLP -> atbAP -> attAP -> at
 * chiral can be 0,1,-1. If chiral == 0 the torsionAngle is a Dihedral angle. Otherwise is a bondLength angle.
 * All angles for this function are in *DEGREES* not radiants!
 * 
 * @param 
 *        atblP (Atom&) = atombondLengthPartner
 *        bondLength (const double) 
 *        atbAP (Atom&) = atombondAnglePartner
 *        bondAngle (const double)
 *        attAP (Atom&) = atomtorsionAnglePartner
 *        tosionAngle (const double)
 *        chiral (const int)
 *        at (Atom&) = modified atom
 * 
 * 
 * 
 */

void IntLoader::zAtomToCartesian(Atom& atbLP, double bondLength, Atom& atbAP,
        double bondAngle, Atom& attAP, double torsionAngle, int chiral, Atom& at) {
    double sine, sine2, a, b, c;
    double cosine;
    vgVector3<double> u1, u2, u3, u4;

    // mathematical and geometrical constants
    double eps = 1e-8;
    double sin_1 = sin(bondAngle / 57.29577951308232088);
    double cos_1 = cos(bondAngle / 57.29577951308232088);
    double sin_2 = sin(torsionAngle / 57.29577951308232088);
    double cos_2 = cos(torsionAngle / 57.29577951308232088);

    vgVector3<double> aCatbAP = atbAP.getCoords();
    vgVector3<double> aCattAP = attAP.getCoords();
    vgVector3<double> aCatbLP = atbLP.getCoords();

    // first, the case where second angle is a dihedral angle
    if (chiral == 0) {
        u1 = (aCatbAP - aCattAP).normalize();
        u2 = (aCatbLP - aCatbAP).normalize();

        u3[0] = u1[1] * u2[2] - u1[2] * u2[1];
        u3[1] = u1[2] * u2[0] - u1[0] * u2[2];
        u3[2] = u1[0] * u2[1] - u1[1] * u2[0];
        cosine = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];

        if (fabs(cosine) < 1.0)
            sine = sqrt(1.0 - (cosine * cosine));
        else {
            ERROR("IntLoader::zAtomToCartesian : Undefined Dihedral Angle",
                    exception);
            sine = sqrt((cosine * cosine) - 1.0);
        }

        u3 /= sine;
        u4[0] = u3[1] * u2[2] - u3[2] * u2[1];
        u4[1] = u3[2] * u2[0] - u3[0] * u2[2];
        u4[2] = u3[0] * u2[1] - u3[1] * u2[0];

        at.setCoords(aCatbLP[0] + bondLength
                * (-u2[0] * cos_1 + u4[0] * sin_1 * cos_2 + u3[0] * sin_1 * sin_2),
                aCatbLP[1] + bondLength
                * (-u2[1] * cos_1 + u4[1] * sin_1 * cos_2 + u3[1] * sin_1 * sin_2),
                aCatbLP[2] + bondLength
                * (-u2[2] * cos_1 + u4[2] * sin_1 * cos_2 + u3[2] * sin_1 * sin_2));
    }// now, the case where second angle is a bondLength angle
    else if ((int)fabs(chiral) == 1) {
        u1 = (aCatbLP - aCattAP).normalize();
        u2 = (aCatbAP - aCatbLP).normalize();

        u3[0] = u1[1] * u2[2] - u1[2] * u2[1];
        u3[1] = u1[2] * u2[0] - u1[0] * u2[2];
        u3[2] = u1[0] * u2[1] - u1[1] * u2[0];
        cosine = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];

        if (fabs(cosine) < 1.0)
            sine2 = 1.0 - (cosine * cosine);
        else
            ERROR("IntLoader::zAtomToCartesian : Defining Atoms Colinear at Atom.",
                exception);
        a = (-cos_2 - cosine * cos_1) / sine2;
        b = (cos_1 + cosine * cos_2) / sine2;
        c = (a * cos_2 + 1. - b * cos_1) / sine2;

        if (c > eps)
            c = chiral * sqrt(c);
        else if (c < -eps) {
            double temp1 = a * u1[0] + b * u2[0];
            double temp2 = a * u1[1] + b * u2[1];
            double temp3 = a * u1[2] + b * u2[2];
            c = sqrt(temp1 * temp1 + temp2 * temp2 + temp3 * temp3);
            a /= c;
            b /= c;
            c = 0;
#ifdef DEBUG
            ERROR("IntLoader::zAtomToCartesian : Sum of bondLength Angles Too Large.",
                    exception);
#endif
        }
        else
            c = 0;
        at.setCoords(aCatbLP[0] + bondLength
                * (a * u1[0] + b * u2[0] + c * u3[0]),
                aCatbLP[1] + bondLength
                * (a * u1[1] + b * u2[1] + c * u3[1]),
                aCatbLP[2] + bondLength
                * (a * u1[2] + b * u2[2] + c * u3[2]));
    }
}

/**
 *   Loads a Group in INT format
 *@param  group reference
 *@return  void
 */

void IntLoader::loadGroup(Group& group) {
    vgVector3<double> tmpV;
    vgVector3<double> tmpV1;

    unsigned long nAtoms;
    input >> nAtoms;
    if (nAtoms < 1)
        ERROR("IntLoader::read(istream&) : number of atoms less than one.",
            exception);

    group.setType(readLine(input)); // read title line  
    for (unsigned long i = 0; i < nAtoms; i++) {
        unsigned long n;
        input >> n;
        n--;
        if (n != i)
            ERROR("IntLoader::loadGroup(Group&) : file format error.", exception);

        string typeName;
        int typeNo;
        Atom at;
        input >> typeName >> typeNo;
        at.setType(typeName);
        at.setNumber(n + 1);

        if (i == 0) {
            group.addAtom(at);
            group[0].setCoords(0, 0, 0); // first atom is placed at the origin
            continue;
        }

        int bondLengthPartner;
        double bondLength;
        input >> bondLengthPartner >> bondLength;
        bondLengthPartner--;

        if (i == 1) {
            group.addAtom(at);
            group[bondLengthPartner].bindOut(group[group.size() - 1]);
            // bind the bond angle partner to the actual atom 
            group[1].setCoords(0, 0, bondLength);
            // second atom is placed along the z-axis
            continue;
        }

        int bondAnglePartner = 0;
        double bondAngle = 0.0;
        double old_bondLength = 0.0;
        double old_bondAngle = 0.0;
        input >> bondAnglePartner >> bondAngle;
        bondAnglePartner--;

        if (i == 2) {
            // third atom is placed in the x,z-plane
            tmpV = group[i - 1].getCoords();
            group.addAtom(at);
            group[bondLengthPartner].bindOut(group[group.size() - 1]);
            group[2].setCoords(bondLength
                    * sin(bondAngle / 57.29577951308232088), 0.0,
                    tmpV[2] - (bondLength * cos(bondAngle / 57.29577951308232088)));
            old_bondLength = bondLength;
            // keep it for the special case i == 3 && bondLengthPartner == 1
            old_bondAngle = bondAngle;
            continue;
        }

        int torsionAnglePartner;
        double torsionAngle;
        int chiral;
        input >> torsionAnglePartner >> torsionAngle;
        input >> chiral;
        torsionAnglePartner--;

        if ((i == 3) && (bondLengthPartner == 1)) {
            // Special case for the 4th atom. It's connected with the first.
            tmpV = group[2].getCoords();
            group[2].setCoords(tmpV[0], tmpV[1], old_bondLength
                    * cos(old_bondAngle / 57.29577951308232088));
        }
        if (i > 2) {
            group.addAtom(at);
            group[bondLengthPartner].bindOut(group[group.size() - 1]);
            tmpV = group[i - 1].getCoords();
            if (tmpV[0] == 0.0) { // Special case : as long as atoms remain linear with the 
                // first two atoms, keep placing them along the z-axis
                tmpV = group[bondLengthPartner].getCoords();
                tmpV1 = group[bondAnglePartner].getCoords();
                double sign;
                if (tmpV[2] > tmpV1[2])
                    sign = 1.0;
                else
                    sign = -1.0;
                group[n].setCoords(bondLength
                        * sin(bondAngle / 57.29577951308232088), 0.0,
                        tmpV[2] - sign * bondLength
                        * cos(bondAngle / 57.29577951308232088));
            } else
                zAtomToCartesian(group[bondLengthPartner], bondLength,
                    group[bondAnglePartner], bondAngle, group[torsionAnglePartner],
                    torsionAngle, chiral, group[n]);
        }
    }
}

/**
 *   Loads a sidechain in INT format. 
 *@param  SideChain reference, AminoAcid pointer
 *@return  void
 */
void IntLoader::loadSideChain(SideChain& sc, AminoAcid* aaRef) {
    loadGroup(sc);
    if (aaRef != NULL)
        sc.setBackboneRef(aaRef);
}

/**
 *    Private helper function to determine if atom is backbone or sidechain. 
 *@param  AminoAcid reference , Atom reference
 *@return  bool
 */
bool IntLoader::inSideChain(const AminoAcid& aa, const Atom& at) {
    if (isBackboneAtom(at.getCode()))
        return false;
    if ((at.getType() == "H") || (at.getType() == "HN")
            || ((at.getType() == "HA") && (!aa.isMember(HA)))
            || (at.getType() == "1HA") || (at.getType() == "1H")
            || (at.getType() == "2H") || (at.getType() == "3H"))
        return false; // special case for GLY H (code HA)

    return true; // rest of aminoacid is its sidechain
}

/**
 *   Loads an AminoAcid in INT format. 
 *@param  AminoAcid reference
 *@return  void
 */

void
IntLoader::loadAminoAcid(AminoAcid& aa) {
    loadGroup(aa);

    if (checkForKeyword(input, "sidechain")) {
        loadSideChain(aa.getSideChain(), &aa);
        if (connect)
            if (aa.getSideChain().size())
                for (unsigned int i = 0; i < aa.size(); i++)
                    if (aa[i].getType().c_str()[0] == 'C') // connect backbone
                    {
                        if (!aa[i].isOutBond(aa.getSideChain()[0]))
                            aa[i].bindOut(aa.getSideChain()[0]);
                        break;
                    }
    } else
        DEBUG_MSG("IntLoader::loadAminoAcid: No sidechain found.");
}

/**
 *   Add atom in aminoacid and returns the pointer from the atom
 *@param   Atom reference
 *@return   AminoAcid pointer,  Atom pointer
 */

Atom&
IntLoader::addToAminoAcid(AminoAcid* aa, Atom* at) {
    if (inSideChain(*aa, *at)) {
        aa->getSideChain().addAtom(*at);
        return ((*aa)[aa->size() - 1]);
    } else {
        aa->addAtom(*at);
        return ((*aa)[aa->sizeBackbone() - 1]);
    }
}

/**
 *   Special case : as long as atoms remain linear with the first
 *                 two atoms, keep placing them along the z-axis
 *@param  int  , double , int  , double , Atom reference
 *@return  void
 */

void
IntLoader::keepPlacingAlongZAxis(int bLP, double bondLength, int bAP,
        double bondAngle, Atom* at) {
    // as long as atoms remain linear with the first
    //  two atoms, keep placing them along the z-axis
    double sign;
    vgVector3<double> tmpV = atomIndex[bLP].getCoords();
    vgVector3<double> tmpV1 = atomIndex[bAP].getCoords();
    if (tmpV[2] > tmpV1[2])
        sign = 1.0;
    else
        sign = -1.0;

    tmpV = atomIndex[bLP].getCoords();
    at->setCoords(bondLength * sin(bondAngle / 57.29577951308232088), 0.0,
            tmpV[2] - sign * bondLength * cos(bondAngle / 57.29577951308232088));
}

/**
 *    Load a spacer from a file in Tinker format. 
 *@param  Spacer reference
 *@return  void
 */

void IntLoader::loadSpacer(Spacer& sp) {
    vgVector3<double> tmpV;
    vgVector3<double> prevAtomCoords;
    vgVector3<double> vSave_i3_bLP1_forSpecialCase;
    bool specialCase_i3_bLP1 = false;
    double bondAngle_forSpecialCase = 0.0;
    double bondLength_forSpecialCase = 0.0;
    unsigned long nAtoms;
    input >> nAtoms;
    if (nAtoms < 1)
        ERROR("IntLoader::loadSpacer() : number of atoms less than one.",
            exception);

    sp.setType(readLine(input)); // read title line from spacer  
    unsigned long i = 0;
    AminoAcid* aa = NULL;
    Atom* at = NULL;

    do {
        unsigned long n;
        input >> n;
        n--;
        if (n != i)
            ERROR("IntLoader::loadSpacer() : file format error.", exception);

        int typeNo;
        string typeName;
        input >> typeName >> typeNo;

        if ((typeName == "N") || (aa == NULL)) // Is there a new aminoacid ?
        {
            aa = new AminoAcid;
            aa->getSideChain().setBackboneRef(aa);
            sp.insertComponent(aa);
        }

        at = new Atom;
        at->setType(typeName);

        if (i == 0) {
            at->setCoords(0.0, 0.0, 0.0); // first atom is placed at the origin
            atomIndex.push_back(addToAminoAcid(aa, at));
            i++;
            continue;
        }

        int bondLengthPartner;
        double bondLength;
        input >> bondLengthPartner >> bondLength;
        bondLengthPartner--;

        if (i == 1) { // second atom is placed along the z-axis
            at->setCoords(0.0, 0.0, bondLength);
            tmpV = at->getCoords();
            // keep the value for the next atom (i == 2)
            atomIndex.push_back(addToAminoAcid(aa, at));
            i++;
            continue;
        }

        int bondAnglePartner;
        double bondAngle;
        input >> bondAnglePartner >> bondAngle;
        bondAnglePartner--;

        if (i == 2) { // third atom is placed in the x,z-plane
            at->setCoords(bondLength * sin(bondAngle / 57.29577951308232088),
                    0.0, tmpV[2] - (bondLength * cos(bondAngle
                    / 57.29577951308232088)));
            atomIndex.push_back(addToAminoAcid(aa, at));
            bondLength_forSpecialCase = bondLength;
            bondAngle_forSpecialCase = bondAngle;
            vSave_i3_bLP1_forSpecialCase = at->getCoords();
            prevAtomCoords = vSave_i3_bLP1_forSpecialCase;
            i++;
            continue;
        }

        int torsionAnglePartner;
        double torsionAngle;
        int chiral;
        input >> torsionAnglePartner >> torsionAngle;
        input >> chiral;
        torsionAnglePartner--;

        if ((i == 3) && (bondLengthPartner == 1))
            specialCase_i3_bLP1 = true;

        if (i > 2) {
            if ((prevAtomCoords[0] == 0.0)) // special case
                keepPlacingAlongZAxis(bondLengthPartner, bondLength,
                    bondAnglePartner, bondAngle, at);
            else
                zAtomToCartesian(atomIndex[bondLengthPartner], bondLength,
                    atomIndex[bondAnglePartner], bondAngle,
                    atomIndex[torsionAnglePartner], torsionAngle, chiral, *at);
            atomIndex.push_back(addToAminoAcid(aa, at));
        }
        i++;
    } while (i < nAtoms);
    setBonds(sp);

    if (specialCase_i3_bLP1)
        sp.getAmino(0)[2].setCoords(vSave_i3_bLP1_forSpecialCase[0], vSave_i3_bLP1_forSpecialCase[1], bondLength_forSpecialCase * cos(bondAngle_forSpecialCase / 57.29577951308232088));

    if (sp.sizeAmino())
        sp.getAmino(0).adjustLeadingN();
}

/**
 *   Load a Ligand from a file in Tinker format. 
 *@param  Ligand reference
 *@return  void
 */
void
IntLoader::loadLigand(Ligand& l) {
    ERROR("Not implemented yet", exception);
}
