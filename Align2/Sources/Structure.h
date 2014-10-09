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


#ifndef __Structure_H__
#define __Structure_H__

#include <SubMatrix.h>
#include <math.h>
#include <string>

namespace Victor { namespace Align2{

    /** @brief   Base class for structural scores.
     * 
     * @Description  
     * @This 
     **/
    class Structure {
    public:

        // CONSTRUCTORS:

        /// Default constructor.

        Structure(SubMatrix *subStr) : subStr(subStr) {
        }

        /// Copy constructor.

        Structure(const Structure &orig) {
            copy(orig);
        }

        /// Destructor.

        virtual ~Structure() {
        }


        // OPERATORS:

        /// Assignment operator.
        Structure& operator =(const Structure &orig);


        // PREDICATES:

        /// Calculate scores to create matrix values.
        virtual double scoringStr(int i, int j) = 0;


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const Structure &orig);

        /// Construct a new "deep copy" of this object.
        virtual Structure* newCopy() = 0;

        /// Reverse template structural components.

        virtual void reverse() {
        }


        // ATTRIBUTES:

        SubMatrix *subStr; ///< Structural substitution matrix.


    protected:


    private:

    };

    // -----------------------------------------------------------------------------
    //                                  Structure
    // -----------------------------------------------------------------------------

    // OPERATORS:

    inline Structure&
            Structure::operator =(const Structure &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // MODIFIERS:

    inline void
    Structure::copy(const Structure &orig) {
        subStr = orig.subStr->newCopy();
    }

}} // namespace

#endif
