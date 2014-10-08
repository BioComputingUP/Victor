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


#include <algorithm>
#include <SimpleBond.h>

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:

/**
 * @Descripton Basic constructor
 * @param mI, maximum in bonds
 * @param mO, maximum out bonds
 */
SimpleBond::SimpleBond(unsigned int mI, unsigned int mO) : inBonds(),
outBonds(), maxIn(mI), maxOut(mO), id("X") {
    PRINT_NAME;
    inBonds.reserve(mI);
    outBonds.reserve(mO);
}

SimpleBond::SimpleBond(const SimpleBond& orig) {
    PRINT_NAME;
    copy(orig);
}

SimpleBond::~SimpleBond() {
    PRINT_NAME;
    while (!inBonds.empty())
        unbindIn(*inBonds.front());
    while (!outBonds.empty())
        unbindOut(*outBonds.front());
}


// PREDICATES:

/**
 *@Description
 *    Checks if this is indirectly bonded to c.
 *    (ie. A to C if A bond B and B bond C) 
 */
bool SimpleBond::isIndirectBond(const SimpleBond& c) const {
    for (unsigned int i = 0; i < sizeInBonds(); i++)
        if (getInBond(i).isBond(c))
            return true;
    for (unsigned int i = 0; i < sizeOutBonds(); i++)
        if (getOutBond(i).isBond(c))
            return true;

    return false;
}

/**
 *@Description
 *    Checks if this is indirectly bonded to c.
 *    (ie. A to C if A bond B and B bond C) */
bool SimpleBond::isIndirectInBond(const SimpleBond& c) const {
    for (unsigned int i = 0; i < sizeInBonds(); i++)
        if (getInBond(i).isBond(c))
            return true;

    return false;
}

/**
 *@Description
 *    Checks if this is indirectly bonded to c.
 *    (ie. A to C if A bond B and B bond C) */
bool SimpleBond::isIndirectOutBond(const SimpleBond& c) const {
    for (unsigned int i = 0; i < sizeOutBonds(); i++)
        if (getOutBond(i).isBond(c))
            return true;

    return false;
}

/**
 *@Description
 *    Checks if this is torsion bond to c.
 *    (ie. A to D if A indirect bond C and C bond D) */
bool SimpleBond::isTorsionBond(const SimpleBond& c) const {
    for (unsigned int i = 0; i < sizeOutBonds(); i++)
        if (getOutBond(i).isIndirectBond(c))
            return true;

    for (unsigned int i = 0; i < sizeInBonds(); i++)
        if (getInBond(i).isIndirectBond(c))
            return true;

    return false;
}

// MODIFIERS:

/**
 * @Description
 *    Sets and in-bond from this to c.
 */
void SimpleBond::bindIn(SimpleBond& c) {
    PRINT_NAME;
    if (isInBond(c))
        return;
    PRECOND((sizeInBonds() < maxIn), exception);
    PRECOND((c.sizeOutBonds() < c.maxOut), exception);
    inBonds.push_back(&c);
    c.outBonds.push_back(this);
}

/**
 * @Description Sets and out-bond from this to c
 */
void SimpleBond::bindOut(SimpleBond& c) {
    PRINT_NAME;
    if (isOutBond(c))
        return;
    PRECOND((sizeOutBonds() < maxOut), exception);
    PRECOND((c.sizeInBonds() < c.maxIn), exception);
    outBonds.push_back(&c);
    c.inBonds.push_back(this);
}

/**
 *@Description
 *    Removes an in-bond from this to c.*/
void SimpleBond::unbindIn(SimpleBond& c) {
    PRINT_NAME;
    pUnbindIn(c);
    c.pUnbindOut(*this);
}

/**
 *@Description
 *    Removes an out-bond from this to c.*/
void SimpleBond::unbindOut(SimpleBond& c) {
    PRINT_NAME;
    pUnbindOut(c);
    c.pUnbindIn(*this);
}

/**
 *@Description
 *    Private method to find the matching in-bond from c to this to remove.*/
void SimpleBond::pUnbindIn(SimpleBond& c) {
    PRINT_NAME;
    vector<SimpleBond*>::iterator iter = find(inBonds.begin(),
            inBonds.end(), &c);
    if (iter == inBonds.end())
        DEBUG_MSG("SimpleBond::pUnbindIn(): WARNING: Bond not found.");
    else
        inBonds.erase(iter);
}

/**
 *@Description  Private method to find the matching out-bond from c to this to remove.
 */
void SimpleBond::pUnbindOut(SimpleBond& c) {
    PRINT_NAME;
    vector<SimpleBond*>::iterator iter = find(outBonds.begin(),
            outBonds.end(), &c);
    if (iter == outBonds.end())
        DEBUG_MSG("SimpleBond::pUnbindOut(): WARNING: Bond not found.");
    else
        outBonds.erase(iter);
}

/** 
 *@Description
 *    Copy operator.
 *    Attention: copy() strips orig from its bonds and attaches them to 
 *               the new bond. 
 */
void SimpleBond::copy(const SimpleBond& orig) {
    PRINT_NAME;
    maxIn = orig.maxIn;
    maxOut = orig.maxOut;
    id = orig.id;

    // Remove old bonds.
    vector<SimpleBond*>::iterator iter = inBonds.begin();
    while (iter != inBonds.end())
        unbindIn(**iter);

    iter = outBonds.begin();
    while (iter != outBonds.end())
        unbindOut(**iter);

    // Create new bonds.
    for (int i = orig.sizeInBonds() - 1; i >= 0; i--) {
        SimpleBond& tmp = const_cast<SimpleBond&> (orig.getInBond(i));
        const_cast<SimpleBond&> (orig).unbindIn(tmp);
        bindIn(tmp);
    }
    for (int i = orig.sizeOutBonds() - 1; i >= 0; i--) {
        SimpleBond& tmp = const_cast<SimpleBond&> (orig.getOutBond(i));
        const_cast<SimpleBond&> (orig).unbindOut(tmp);
        bindOut(tmp);
    }
}

