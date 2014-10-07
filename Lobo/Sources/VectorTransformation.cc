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
 *@Class:              VectorTransformation
 *@Description:This class allows to store transformation steps for transforming a 
 *     vector v into v' according to the series of steps performed earlier.
 *      
 */

// Includes:
#include <VectorTransformation.h>
#include <IoTools.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;

// CONSTRUCTORS/DESTRUCTOR:

/**
 *@Description basic constructor
 */
VectorTransformation::VectorTransformation() : rot(0), trans(0) {
    addNewElem();
}

/**
 *@Description Constructor that copies a original vector transformation 
 */
VectorTransformation::VectorTransformation(const VectorTransformation& orig) {
    this->copy(orig);
}

/**
 *@Description basic destructor
 */
VectorTransformation::~VectorTransformation() {
    PRINT_NAME;
    rot.clear();
    trans.clear();
}


// PREDICATES:

/**
 *@Description transform the original vector
 *@param  original vector(vgVector3<float> )
 *@return  corresponding value after the transformation(vgVector3<float> )
 */
vgVector3<float> VectorTransformation::transform(vgVector3<float> orig) {
    PRECOND(rot.size() == trans.size(), exception);

    for (unsigned int i = rot.size(); i > 0; i--) {
        orig = rot[i - 1] * orig;
        orig += trans[i - 1];
    }

    return orig;
}



// MODIFIERS:

/**
 *@Description Adds rotation 
 *@param   matrix to rotate(vgMatrix3<float>)
 *@return  changes are made internally(void)
 */
void VectorTransformation::addRot(vgMatrix3<float> rm) {
    PRECOND(rot.size() == trans.size(), exception);

    if (trans[trans.size() - 1].length() != 0)
        addNewElem();

    rot[rot.size() - 1] = rot[rot.size() - 1] * rm;

}

/**
 *@Description Adds translation 
 *@param   matrix to rotate(vgVector3<float>)
 *@return  changes are made internally(void)
 */
void VectorTransformation::addTrans(vgVector3<float> t) {
    PRECOND(rot.size() == trans.size(), exception);

    vgMatrix3<float> idM(1);
    if (rot[rot.size() - 1] != idM)
        addNewElem();

    trans[trans.size() - 1] += t;

}

/**
 *@Description Clears the object data
 *@param  none
 *@return  changes are made internally(void)
 */
void VectorTransformation::clear() {
    rot.clear();
    trans.clear();
    addNewElem();
}

/**
 *@Description Copies the original vector transformation
 *@param  original vector transformation
 *@return   changes are made internally(void)
 */
void VectorTransformation::copy(const VectorTransformation& orig) {
    PRINT_NAME;
    rot.clear();
    trans.clear();

    for (unsigned int i = 0; i < orig.rot.size(); i++)
        rot.push_back(orig.rot[i]);
    for (unsigned int i = 0; i < orig.trans.size(); i++)
        trans.push_back(orig.trans[i]);
}


// OPERATORS:

/**
 *@Description Assigns a transformation vector into another one.
 *@param  reference to the original transformation vector(VectorTransformation& )
 *@return  reference to the new transformation vector(VectorTransformation& )
 */
VectorTransformation& VectorTransformation::operator=(const VectorTransformation& orig) {
    if (&orig != this)
        copy(orig);
    return *this;
}

