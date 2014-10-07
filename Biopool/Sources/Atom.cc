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
#include <Atom.h>
#include <Group.h>
#include <matrix3.h>
#include <math.h>
#include <IoTools.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:

Atom::Atom(unsigned int mI, unsigned int mO) : SimpleBond(mI, mO),
superior(NULL), type(X), coords(0, 0, 0), Bfac(0.0), trans(0, 0, 0), rot(1),
modified(false) {
    PRINT_NAME;
}

Atom::Atom(const Atom& orig) {
    PRINT_NAME;
    this->copy(orig);
}

Atom::~Atom() {
    PRINT_NAME;
}

// PREDICATES:

double
Atom::distance(Atom& other) {
    sync();
    other.sync();

    double tmp = sqrt(
            (coords[0] - other.coords[0]) * (coords[0] - other.coords[0])
            + (coords[1] - other.coords[1]) * (coords[1] - other.coords[1])
            + (coords[2] - other.coords[2]) * (coords[2] - other.coords[2])
            );
    return tmp;
}


// MODIFIERS:

void
Atom::setCoords(double _x, double _y, double _z) {
    vgVector3<double> tmp(_x, _y, _z);
    setCoords(tmp);
}

void
Atom::setCoords(vgVector3<double> c) {
    sync();
    vgVector3<double> oldTrans = trans;
    coords = c;

    if (isNotFirstAtomInStructure())
        trans = coords - getInBond(0).coords;
    else
        if (hasSuperior())
        trans = coords - getSuperior().getTrans();
    else
        trans = coords;

    // adjust coords of following atoms:
    for (unsigned int i = 0; i < sizeOutBonds(); i++) {
        getOutBond(i).setTrans(getOutBond(i).getTrans() + (oldTrans - trans));
    }

    setModified();
}

void
Atom::copy(const Atom& orig) {
    PRINT_NAME;
    SimpleBond::copy(orig);

    setNumber(orig.getNumber());

    superior = orig.superior;
    type = orig.type;
    coords = orig.coords;
    Bfac = orig.Bfac;

    trans = orig.trans;
    rot = orig.rot;
    modified = orig.modified;
}

void
Atom::sync() // synchronize coords with structure
{
    if (inSync())
        return;

    vgMatrix3<double> tmpMatrix(1);
    vgVector3<double> supTrans(0, 0, 0);

    if (isNotFirstAtomInStructure())
        getInBond(0).sync();
    else
        if (hasSuperior()) { // use superior's trans & rot
        supTrans = getSuperior().getTrans();
        rot = getSuperior().getRot() * rot;
        const_cast<Group&> (getSuperior()).setRot(tmpMatrix);
    }

    trans = rot * trans;
    coords = trans + supTrans; // set the relative position

    if (isNotFirstAtomInStructure())
        coords += getInBond(0).getCoords(); // make absolute position

    propagateRotation();
    modified = false;
}

void
Atom::setModified() {
    if (modified)
        return;
    modified = true;

    for (unsigned int i = 0; i < sizeOutBonds(); i++)
        getOutBond(i).setModified();

}

void
Atom::setUnModified() {
    modified = false;
}


// OPERATORS:

Atom&
        Atom::operator=(const Atom& orig) {
    PRINT_NAME;

    if (&orig != this)
        copy(orig);
    return *this;
}


// HELPERS:

void
Atom::propagateRotation() {
    vgMatrix3<double> tmpMatrix(1);
    if (rot == tmpMatrix)
        return;

    for (unsigned int i = 0; i < sizeOutBonds(); i++)
        if (!hasSuperior() ||
                (!((getSuperior().getType() == "PRO") && (getCode() == CD)
                && (getOutBond(i).getCode() == N))
                && !((getSuperior().getType() == "PHE") && (getCode() == CE2)
                && (getOutBond(i).getCode() == CZ))
                && !((getSuperior().getType() == "TYR") && (getCode() == CE2)
                && (getOutBond(i).getCode() == CZ))
                && !((getSuperior().getType() == "TRP") && (getCode() == NE1)
                && (getOutBond(i).getCode() == CE2))
                && !((getSuperior().getType() == "TRP") && (getCode() == CZ3)
                && (getOutBond(i).getCode() == CH2))
                && !((getSuperior().getType() == "HIS" && (getCode() == CE1)
                && (getOutBond(i).getCode() == NE2)))
                )
                ) {
            getOutBond(i).addRot(rot);
        }
    // Cyclic structures (PRO, PHE, TYR, TRP & HIS) are special cases 
    // when propagating the rotation, because the the re-conjunction of 
    // two branches the rotation would be added twice.
    rot = tmpMatrix;
}

inline bool
Atom::isNotFirstAtomInStructure() {
    if (sizeInBonds() && !((getSuperior().getType() == "PRO")
            && (getCode() == N) && (getInBond(0).getCode() == CD)))
        // NB: Proline's N is a special case because unbound PRO 
        // (ie. not part of a Spacer) has CD as the 1st inBond of N,
        // which would crash the program
        return true;
    else
        return false;
}
