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

#ifndef __Sec_H__
#define __Sec_H__

#include <AlignmentData.h>
#include <Structure.h>

namespace Victor { namespace Align2{

    /** @brief   Calculate structural scores with info derived from secondary
     *                  structure.
     * 
     *   

     **/
    class Sec : public Structure {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        Sec(SubMatrix *subStr, AlignmentData *ad, double cSec);

        /// Copy constructor.
        Sec(const Sec &orig);

        /// Destructor.
        virtual ~Sec();


        // OPERATORS:

        /// Assignment operator.
        Sec& operator =(const Sec &orig);


        // PREDICATES:

        /// Calculate structural scores to create matrix values.
        virtual double scoringStr(int i, int j);


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const Sec &orig);

        /// Construct a new "deep copy" of this object.
        virtual Sec* newCopy();

        /// Reverse template secondary structure.
        virtual void reverse();


    protected:


    private:

        // ATTRIBUTES:

        string sec1; ///< Target secondary structure.
        string sec2; ///< Template secondary structure.
        double cSec; ///< Coefficient for secondary structure alignment.

    };

}} // namespace

#endif
