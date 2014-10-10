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
 *@Class:              ranking_helper2
 * 
 *@Description:Just contains an int and double value which gets sorted (by a STL set)
 *    in order to determine the ranking of low rms solutions with special 
 *    filters. This container differs from the ranking_helper by the fact,
 *    that it has the < operator overloaded differently 
 *      
 */

//Includes
#include<ranking_helper2.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Victor;

using namespace Victor::Lobo;
// CONSTRUCTORS/DESTRUCTOR:

/**
 *@Description constructor  sets the values of the index that contains the rms ranking of the solution , 
 * and the value from a filter (like the propensity or the collision)
 *@param index(int) , value(double)
 */
ranking_helper2::ranking_helper2(int ind, double val) {
    index = ind;
    value = val;
}

/**
 *@Description Constructor that copies the info from another object
 *@param  reference to original object(const ranking_helper2& )
 */
ranking_helper2::ranking_helper2(const ranking_helper2& c) {
    this->copy(c);
}

// PREDICATES:

// MODIFIERS:

// OPERATORS:

/**
 *@Description Operator that allows to verify if its lower than other
 *@param  reference to the object(const ranking_helper2 &)
 *@return  result of the verification
 */
bool ranking_helper2::operator<(const ranking_helper2 &name) const {
    return value < name.get_value();
}

/**
 *@Description Operator that allows to assign one object to other
 *@param  reference to the object(const ranking_helper2 &)
 *@return  result of the verification
 */
ranking_helper2& ranking_helper2::operator=(const ranking_helper2& orig) {
    if (&orig != this) {
        copy(orig);
    }
    return *this;
}

/**
 *@Description copies the info from another object
 *@param  reference to original object(const ranking_helper2& )
 *@return changes are made internally(void)
 */
void ranking_helper2::copy(const ranking_helper2 & c) {
    index = c.index;
    value = c.value;
}

// HELPERS:
