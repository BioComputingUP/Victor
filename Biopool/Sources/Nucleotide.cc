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
#include <Nucleotide.h>
#include <NucleotideCode.h>
#include <Debug.h>
#include <limits.h>
#include <float.h>

// Global constants, typedefs, etc. (to avoid):
using namespace Biopool;

/**
 *@Description Constructor
 *@param none
 */
Nucleotide::Nucleotide() : Group(1, 1), type(XX), icc() {
}

/**
 *@Description Constructor
 *@param Nucleotide
 */
Nucleotide::Nucleotide(const Nucleotide& orig) {
    PRINT_NAME;
    this->copy(orig);
}

/**
 *@Description DESTRUCTOR
 *@param none
 */
Nucleotide::~Nucleotide() {
    PRINT_NAME;
}
/**
 *@Description Returns the atom corresponding to N, 
 * aminoacids have only a single possible, hard coded, open in-bond
 *@param unsigned int 
 *@return const Atom& 
 */


// MODIFIERS:

/**
 *@Description  Copies an aa
 *@param const Nucleotide& (copy from the orig)
 *@return  void
 */
void Nucleotide::copy(const Nucleotide& orig) {
    PRINT_NAME;
    Group::copy(orig);

    type = orig.type;

    // set absolute position to orig's:
    if (orig[0].sizeInBonds()) {
        setTrans(const_cast<Nucleotide&> (orig)[0].getInBond(0).getCoords());
    }
}

/**
 *@Description  Clone the aa
 *@param none
 *@return  Component* 
 */
Component* Nucleotide::clone() {
    Nucleotide* tmp = new Nucleotide;
    tmp->copy(*this);
    return tmp;
}

// OPERATORS:

/**
 *@Description  Operator =, assign the aa
 *@param Nucleotide reference
 *@return  Nucleotide
 */
Nucleotide& Nucleotide::operator=(const Nucleotide& orig) {
    PRINT_NAME;
    if (&orig != this)
        copy(orig);
    return *this;
}


