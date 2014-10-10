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

#ifndef __Traceback_H__
#define __Traceback_H__

#include <Debug.h>
#include <iostream>

namespace Victor { namespace Align2{

    /** @brief    Reconstruct the path in the alignment matrix.
     * 
     *   

     **/
    class Traceback {
    public:

        // CONSTRUCTORS:

        /// Default constructor.

        Traceback() : i(-1), j(-1) {
        }

        /// Constructor assigning i and j.

        Traceback(int i, int j) : i(i), j(j) {
        }

        /// Copy constructor.

        Traceback(const Traceback &orig) {
            copy(orig);
        }

        /// Destructor.

        virtual ~Traceback() {
        }


        // OPERATORS:

        /// Assignment operator.
        Traceback& operator =(const Traceback &orig);

        /// Comparison operator.
        friend bool operator ==(const Traceback &left, const Traceback &right);


        // PREDICATES:

        /// Check if object is invalid.
        static bool isInvalidTraceback(const Traceback &tb);

        /// Return invalid position.
        static Traceback getInvalidTraceback();

        /// Return true if other and this object are equivalent.
        virtual bool compare(const Traceback &other) const;


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const Traceback &orig);

        /// Construct a new "deep copy" of this object.
        virtual Traceback* newCopy();


        // ATTRIBUTES:

        int i; ///< Position (row).
        int j; ///< Position (column).


    protected:


    private:

    };

    // -----------------------------------------------------------------------------
    //                                  Traceback
    // -----------------------------------------------------------------------------

    // OPERATORS:

    inline Traceback&
            Traceback::operator =(const Traceback &orig) {
        if (&orig != this)
            copy(orig);
        POSTCOND((orig == *this), exception);
        return *this;
    }

    inline bool
    operator ==(const Traceback &left, const Traceback &right) {
        return left.compare(right);
    }


    // PREDICATES:

    inline bool
    Traceback::isInvalidTraceback(const Traceback &tb) {
        return ((tb.i < 0) || (tb.j < 0));
    }

    inline Traceback
    Traceback::getInvalidTraceback() {
        Traceback tb;
        tb.i = -1;
        tb.j = -1;
        return tb;
    }

    inline bool
    Traceback::compare(const Traceback &other) const {
        return ((i == other.i) && (j == other.j));
    }


    // MODIFIERS:

    inline void
    Traceback::copy(const Traceback &orig) {
        i = orig.i;
        j = orig.j;
    }

    inline Traceback*
    Traceback::newCopy() {
        Traceback *tmp = new Traceback(*this);
        return tmp;
    }

}} // namespace

#endif
