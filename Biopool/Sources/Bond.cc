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
 
 *@Description:     Defines chemical and abstract bonds between objects 
 *                  which are compositions of atoms. Eg. 'bond' between 
 *                  two amino acids.
 *                  Attention: copy() strips orig from its bonds and 
 *                  attaches them to the new bond. 
 * 
 */

#include <algorithm>
#include <Bond.h>

using namespace Biopool;

Bond::Bond(unsigned int mI, unsigned int mO) : SimpleBond(mI, mO),
inRef(), outRef() {
    PRINT_NAME;
    inRef.reserve(mI);
    outRef.reserve(mO);
}

Bond::Bond(const Bond& orig) {
    PRINT_NAME;
    copy(orig);
}

Bond::~Bond() {
    PRINT_NAME;
    while (!inRef.empty())
        inRef.pop_back();
    while (!outRef.empty())
        outRef.pop_back();
}


// PREDICATES

/**
 *@Description  Returns reference to the i-th open in-bond atom.
 *    NB: This is undefined (ie. returns error) for Bond, but is needed for
 *    later structures (AminoAcid, etc.).
 *@param unsigned int 
 *@return  Atom reference
 */

Atom&
Bond::getOpenInBondRef(unsigned int n) {
    ERROR("getOpenInBondRef() undefined for this class.", exception);
    return *inRef[n];
}

const Atom&
Bond::getOpenInBondRef(unsigned int n) const {
    ERROR("getOpenInBondRef() undefined for this class.", exception);
    return *inRef[n];
}

/**
 *@Description  Returns reference to the i-th open out-bond atom.
 *    NB: This is undefined (ie. returns error) for Bond, but is needed for
 *    later structures (AminoAcid, etc.).
 *@param unsigned int
 *@return   atom reference
 */
Atom& Bond::getOpenOutBondRef(unsigned int n) {
    ERROR("getOpenOutBondRef() undefined for this class.", exception);
    return *outRef[n];
}

const Atom&
Bond::getOpenOutBondRef(unsigned int n) const {
    ERROR("getOpenOutBondRef() undefined for this class.", exception);
    return *outRef[n];
}

// MODIFIERS

/**
 *@Description   In-going bond from this to c, from this' atom _this to c's _other atom. 
 *@param  atom reference, bond reference, atom reference
 *@return    void
 */

void Bond::bindIn(Atom& _this, Bond& c, Atom& _other) {
    if (isInBond(c)) // already bound
    {
        _this.bindIn(_other); // check if atoms are already bound as well
        return;
    };

    PRECOND((sizeInBonds() < maxIn), exception);
    PRECOND((c.sizeOutBonds() < c.maxOut), exception);

    inBonds.push_back(&c);
    c.outBonds.push_back(this);
    inRef.push_back(&_this);
    c.outRef.push_back(&_other);

    _this.bindIn(_other); // bind the atoms
}

/**
 *@Description   Out-going bond from this to c, from this' atom _this to c's _other atom. 
 *@param   atom reference, bond reference, atom reference
 *@return  void  
 */
void
Bond::bindOut(Atom& _this, Bond& c, Atom& _other) {
    if (isOutBond(c)) // already bound
    {
        _this.bindOut(_other); // check if atoms are already bound as well
        return;
    };

    PRECOND((sizeOutBonds() < maxOut), exception);
    PRECOND((c.sizeInBonds() < c.maxIn), exception);

    outBonds.push_back(&c);
    c.inBonds.push_back(this);
    outRef.push_back(&_this);
    c.inRef.push_back(&_other);

    _this.bindOut(_other); // bind the atoms
}

/**
 *@Description   Remove in-going bond from this to c. 
 *@param     bond reference 
 *@return    void
 */
void Bond::unbindIn(Bond& c) {
    PRINT_NAME;
    pUnbindIn(c, true);
    c.pUnbindOut(*this);
}

/**
 *@Description   Remove out-going bond from this to c.
 *@param   bond reference 
 *@return    void
 */
void Bond::unbindOut(Bond& c) {
    PRINT_NAME;
    pUnbindOut(c, true);
    c.pUnbindIn(*this);
}


/**
 *@Description    Private helper method for removing in-bonds. 
 *@param   bond reference , bool
 *@return    void
 */
//
// ----------------------------------------------------------------------------

void Bond::pUnbindIn(Bond& c, bool unbind) {
    PRINT_NAME;
    int index = 0;

    for (unsigned int i = 0; i < inBonds.size(); i++)
        if (*inBonds[i] == dynamic_cast<SimpleBond&> (c)) {
            index = i;
            inBonds.erase(inBonds.begin() + i);
            break;
        }

    if (index == 0)
        DEBUG_MSG("Bond::pUnbindIn(Bond& c): WARNING: Bond not found.");
    else {
        if (unbind) {
            // find bound atom in c
            int index2 = 0;
            for (unsigned int i = 0; i < c.outRef.size(); i++)
                if (c.outRef[i] == inRef[index]) {
                    index2 = i;
                    break;
                }
            inRef[index]->unbindIn(*c.outRef[index2]);
        }
        inRef.erase(inRef.begin() + index); // delete refIn  
    }
}

/**
 *@Description   Private helper method for removing out-bonds. 
 *@param   bond reference , bool
 *@return    void
 */
void Bond::pUnbindOut(Bond& c, bool unbind) {
    PRINT_NAME;
    int index = 0;

    for (unsigned int i = 0; i < outBonds.size(); i++)
        if (*outBonds[i] == dynamic_cast<SimpleBond&> (c)) {
            index = i;
            outBonds.erase(outBonds.begin() + i);
            break;
        }

    if (index == 0)
        DEBUG_MSG("Bond::pUnbindOut(Bond& c): WARNING: Bond not found.");
    else {
        if (unbind) {
            // find bound atom in c
            int index2 = 0;
            for (unsigned int i = 0; i < c.inRef.size(); i++)
                if (c.inRef[i] == outRef[index]) {
                    index2 = i;
                    break;
                }
            outRef[index]->unbindOut(*c.inRef[index2]);
        }
        outRef.erase(outRef.begin() + index); // delete refOut  
    }
}

/**
 *@Description   Copy operator. Attention: copy() strips orig from its bonds and  attaches them to the new bond. 
 *@param   bond reference 
 *@return    void
 */
void
Bond::copy(const Bond& orig) {
    PRINT_NAME;
    SimpleBond::copy(orig);

    inRef = orig.inRef;
    outRef = orig.outRef;
}

