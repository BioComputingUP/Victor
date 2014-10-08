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
#include <Protein.h>
#include <iostream>
using namespace std;
using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:

Protein::Protein() : Polymer(1, 1) {
    PRINT_NAME;
}

Protein::Protein(const Protein& orig) {
    PRINT_NAME;
    this->copy(orig);
}

Protein::~Protein() {
    PRINT_NAME;
}

// PREDICATES:

/**
 * @Description Return a pointer to the Spacer with the requested chainID
 * @param c (char), the chainID
 * @return A pointer to the Spacer
 */
Spacer* 
Protein::getSpacer(char c) {
    unsigned int n = getChainNum(c);
    return getSpacer(n);
}

/**
 * @Description Get the Polymer (Spacer + LigandSet) by an index.
 * @param n (unsigned int), the index of the Polymer in the Components vector 
 * @return The reference to the Polymer
 */
Polymer&
Protein::getPolymer(unsigned int n) {
    if (n > components.size() - 1)
        ERROR("Index out of range", exception);
    return *(dynamic_cast<Polymer*> (components[n]));
}

/**
 * @Description Return a pointer to the Spacer by an index in the Components vector.
 * @param n (unsigned int), the index of the Polymer in the Components vector
 * @return A pointer to the Spacer
 */
Spacer*
Protein::getSpacer(unsigned int n) {
    if (n > sizeProtein() - 1)
        ERROR("Index out of range", exception);
    Polymer& p = getPolymer(n);
    return &(dynamic_cast<Spacer&> (p[0]));
}

/**
 * @Description Return a pointer to the LigandSet with the requested chainID
 * @param c (char), the chainID
 * @return A pointer to the LiganSet
 */
LigandSet* //return a pointer to the ligandSet (if exists) with the correct chainID
Protein::getLigandSet(char c) {
    unsigned int n = getChainNum(c);
    return getLigandSet(n);
}
/**
 * @Description Return a pointer to the LiganSet by an index in the Components vector.
 * @param n (unsigned int), the index of the Polymer in the Components vector
 * @return A pointer to the LiganSet
 */
LigandSet*
Protein::getLigandSet(unsigned int n) {
    if (n > sizeProtein() - 1)
        ERROR("Index out of range", exception);
    Polymer& p = getPolymer(n);
    if (p.size() < 2)
        return NULL; //No Ligands for this chain
    return &(dynamic_cast<LigandSet&> (p[1]));
}
/**
 * @Description Returns the chain index by chainID
 * @param c (char) chainID
 * @return The chainID index in the chains vector
 */
unsigned int
Protein::getChainNum(char c) {
    for (unsigned int i = 0; i < chains.size(); i++)
        if (chains[i] == c)
            return i;
    ERROR("Chain not found", exception);
}
/**
 * @Description Returns the chainID by index
 * @param i (unsigned int) index
 * @return The chainID (char)
 */
char
Protein::getChainLetter(unsigned int i) {
    if (i > chains.size())
        ERROR("Index out of range", exception);
    return chains[i];
}

// MODIFIERS:

void
Protein::removeComponent(Component* c) {
    Polymer::removeComponent(c);
    setModified();
}

void
Protein::deleteComponent(Component* c) {
    Polymer::deleteComponent(c);
    setModified();
}

void //not finished because Component::copy is not finished too
Protein::copy(const Protein& orig) {
    PRINT_NAME;
    chains = orig.chains;
    Polymer::copy(orig);
}

Protein*
Protein::clone() {
    Protein* tmp = new Protein;
    tmp->copy(*this);
    return tmp;
}

 /**
 *@Description Insert a Polymer in the Protein.
 *               It is related to a single chain, and must
 *               contain a Spacer and possibly a LigandSet
 *@param p (Component*)
 */ void
Protein::insertComponent(Component* p) {
    if (p->hasSuperior())
        ERROR("Component does have a superior", exception);
    if (p->getClassName() != "Polymer")
        ERROR("Protein can insert directly only Polymer", exception);
    Polymer::insertComponent(p);
}

// OPERATORS:

Protein&
        Protein::operator=(const Protein& orig) {
    PRINT_NAME;
    if (&orig != this)
        copy(orig);
    return *this;
}
// HELPERS:

void
Protein::printComponents() {
    ERROR("Not finished yet", exception);
}
