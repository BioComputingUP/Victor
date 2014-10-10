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


#ifndef __GapFunction_H__
#define __GapFunction_H__

#include <Debug.h>
#include <iostream>

namespace Victor { namespace Align2{

    /** @brief Base class for gap functions.
     * 
     *   
     * @This 
     **/
    class GapFunction {
    public:

        // CONSTRUCTORS:

        /// Default constructor.

        GapFunction() {
        }

        /// Copy constructor.

        GapFunction(const GapFunction &orig) {
            copy(orig);
        }

        /// Destructor.

        virtual ~GapFunction() {
        }


        // OPERATORS:

        /// Assignment operator.
        GapFunction& operator =(const GapFunction &orig);


        // PREDICATES:

        /// Return open gap penalty for template position p.
        virtual double getOpenPenalty(int p) = 0;

        /// Return extension gap penalty for template position p.
        virtual double getExtensionPenalty(int p) = 0;


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const GapFunction &orig);

        /// Construct a new "deep copy" of this object.
        virtual GapFunction* newCopy() = 0;

        /// Set open gap penalty.
        virtual void setOpenPenalty(double pen) = 0;

        /// Set extension gap penalty.
        virtual void setExtensionPenalty(double pen) = 0;


    protected:


    private:

    };

    // -----------------------------------------------------------------------------
    //                                 GapFunction
    // -----------------------------------------------------------------------------

    // OPERATORS:

    inline GapFunction&
            GapFunction::operator =(const GapFunction &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }


    // MODIFIERS:

    inline void
    GapFunction::copy(const GapFunction &orig) {
    }

}} // namespace

#endif
