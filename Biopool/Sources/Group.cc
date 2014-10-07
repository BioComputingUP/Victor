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
#include <Group.h>
#include <IoTools.h>
#include <limits.h>
#include <float.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:

Group::Group(unsigned int mI, unsigned int mO) : Monomer(mI, mO),
atoms(), trans(0, 0, 0), rot(1) {
}

Group::Group(const Group& orig) {
    this->copy(orig);
}

Group::~Group() {
    PRINT_NAME;
}

// PREDICATES:

// MODIFIERS:

void
Group::resetBoundaries() {
    if (size() < 1)
        return;

    vgVector3<double> tmpV(DBL_MAX - 1, DBL_MAX - 1, DBL_MAX - 1);
    lowerBound = tmpV;
    upperBound = -tmpV;

    for (unsigned int i = 0; i < atoms.size(); i++) {
        for (unsigned int j = 0; j < 3; j++)
            if (atoms[i].getCoords()[j] < lowerBound[j])
                lowerBound[j] = atoms[i].getCoords()[j];
        for (unsigned int j = 0; j < 3; j++)
            if (atoms[i].getCoords()[j] > upperBound[j])
                upperBound[j] = atoms[i].getCoords()[j];
    }

}

void
Group::setModified() {
    if (modified)
        return;
    modified = true;
    if (hasSuperior())
        getSuperior().setModified();

    if (atoms.size())
        atoms[0].setModified();
}

void
Group::removeAtom(Atom& a) {
    PRINT_NAME;
    PRECOND(atoms.size() > 0, exception);
    if (atoms[atoms.size() - 1] == a) {
        atoms.pop_back();
        a.Atom::setSuperior(NULL);
        Group::setModified();
    } else {
        vector<Atom>::iterator iter = find(atoms.begin(),
                atoms.end(), a);
        if (iter == atoms.end())
            DEBUG_MSG("Group::removeAtom(): WARNING: Atom not found.");
        else
            atoms.erase(iter);
    }
}

void
Group::copy(const Group& orig) {
    PRINT_NAME;
    Monomer::copy(orig);

    while (atoms.size() > 0)
        atoms.pop_back();
    // this solution is not perfectly clean, 
    // but avoids orig getting stripped of its bonds
    // unfortunately there's no (obvious) better solution.
    // this has no bonds after this operation (they will be added below)
    for (unsigned int i = 0; i < orig.atoms.size(); i++) {
        atoms.push_back(orig.atoms[i]);
        const_cast<Atom&> (orig.atoms[i]) = atoms[atoms.size() - 1];
    }

    // copy all bonds
    for (unsigned int i = 0; i < size(); i++)
        for (unsigned int j = 0; j < orig[i].sizeOutBonds(); j++) {
            unsigned int tmp = orig.pGetIndex(orig[i].getOutBond(j).getNumber());
            // in order to copy all bonds from orig the index has to be computed 
            if (tmp < atoms.size())
                (*this)[i].bindOut((*this)[tmp]);
        }

    trans = orig.trans;
    rot = orig.rot;
    modified = orig.modified;

    id = orig.id;
    // not initializing Bond in order to avoid copying the inBond 
    // and outBond references, which would be invalid for the copy

    // let superior reference new group:
    for (unsigned int i = 0; i < atoms.size(); i++)
        atoms[i].Atom::setSuperior(this);
}

void
Group::sync() // synchronize coords with structure
{
    if (!modified)
        return;

    for (unsigned int i = 0; i < atoms.size(); i++)
        atoms[i].sync();

    modified = false;
    vgMatrix3<double> tmp(1);
    rot = tmp;

    resetBoundaries();
}

void
Group::addAtom(Atom& a) {
    PRINT_NAME;
    atoms.push_back(a);
    atoms[atoms.size() - 1].Atom::setSuperior(this);

    for (unsigned int j = 0; j < 3; j++)
        if (a.getCoords()[j] < lowerBound[j])
            lowerBound[j] = a.getCoords()[j];
    for (unsigned int j = 0; j < 3; j++)
        if (a.getCoords()[j] > upperBound[j])
            upperBound[j] = a.getCoords()[j];
}


// OPERATORS:

Group&
        Group::operator=(const Group& orig) {
    PRINT_NAME;
    if (&orig != this)
        copy(orig);
    return *this;
}

// HELPERS:

Atom*
Group::pGetAtom(const AtomCode& ac) const {
    for (unsigned int i = 0; i < atoms.size(); i++)
        if (getAtom(i).getCode() == ac)
            return &(const_cast<Atom&> (getAtom(i)));
    return NULL;
}

unsigned int
Group::pGetIndex(const unsigned long cmp) const {
    for (unsigned int i = 0; i < atoms.size(); i++) {
        if ((*this)[i].getNumber() == cmp)
            return i;
    }
    return size() + 1;
}


