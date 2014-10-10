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


#ifndef __Ss2_H__
#define __Ss2_H__

#include <AlignmentData.h>
#include <Structure.h>
#include <Ss2Input.h>

namespace Victor { namespace Align2{

    /** @brief Calculate structural scores with info derived from PSI-PRED.  
     * 
     *   

     **/
    class Ss2 : public Structure {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        Ss2(SubMatrix *subStr, AlignmentData *ad, Ss2Input *psipred1,
                Ss2Input *psipred2, double cSs2);

        /// Copy constructor.
        Ss2(const Ss2 &orig);

        /// Destructor.
        virtual ~Ss2();


        // OPERATORS:

        /// Assignment operator.
        Ss2& operator =(const Ss2 &orig);


        // PREDICATES:

        /// Calculate structural scores to create matrix values.
        virtual double scoringStr(int i, int j);


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const Ss2 &orig);

        /// Construct a new "deep copy" of this object.
        virtual Ss2* newCopy();

        /// Reverse template secondary structure.
        virtual void reverse();


    protected:


    private:

        // ATTRIBUTES:

        string sec1; ///< Target secondary structure.
        string sec2; ///< Template secondary structure.
        Ss2Input *psipred1; ///< Target PSI-PRED input file.
        Ss2Input *psipred2; ///< Template PSI-PRED input file.
        double cSs2; ///< Coefficient for PSI-PRED prediction.

    };

}} // namespace

#endif
