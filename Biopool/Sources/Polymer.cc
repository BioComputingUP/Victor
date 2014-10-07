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
#include <Polymer.h>
#include <Debug.h>
#include <AminoAcid.h> // for testing only

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:

Polymer::Polymer(unsigned int mI, unsigned int mO) : Component(mI, mO) {
    PRINT_NAME;
}

Polymer::Polymer(const Polymer& orig) {
    PRINT_NAME;
    this->copy(orig);
}

Polymer::~Polymer() {
    PRINT_NAME;
    if (hasSuperior())
        getSuperior().removeComponent(this);
    setSuperior(NULL);
    while (size() > 0)
        deleteComponent(components[0]);
}

// PREDICATES:

// MODIFIERS:

/**
 *@Description 
 *@param 
 */
void Polymer::copy(const Polymer& orig) {
    PRINT_NAME;
    Component::copy(orig);
}

/**
 *@Description 
 *@param 
 */
Component* Polymer::clone() {
    return new Polymer;
}

/**
 *@Description 
 *@param 
 */
void Polymer::insertComponent(Component* c) {
    PRECOND(c != NULL, exception);
    c->setSuperior(this);
    components.push_back(c);
}

/**
 *@Description 
 *@param 
 */
void Polymer::removeComponent(Component* c) {
    if (c == NULL)
        return;
    for (unsigned int i = 0; i < size(); i++)
        if (components[i] == c) {
            components[i]->setSuperior(NULL);
            components.erase(components.begin() + i);
            return;
        }

    DEBUG_MSG("Polymer::removeComponent: component not found.");
}

/**
 *@Description 
 *@param 
 */
void Polymer::removeComponentFromIndex(unsigned int i) {
    if (i > size())
        DEBUG_MSG("Index out of bound");
    components[i]->setSuperior(NULL);
    components.erase(components.begin() + i);
}

/**
 *@Description 
 *@param 
 */
void Polymer::deleteComponent(Component* c) {
    if (c == NULL)
        return;
    for (unsigned int i = 0; i < size(); i++)
        if (components[i] == c) {

            components[i]->setSuperior(NULL);
            components.erase(components.begin() + i);
            delete c;
            return;
        }
    DEBUG_MSG("Polymer::removeComponent: component not found.");
}


// OPERATORS:

// HELPERS:



