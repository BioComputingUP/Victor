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


#ifndef __Threading_H__
#define __Threading_H__

#include <AlignmentData.h>
#include <Structure.h>
#include <ThreadingInput.h>

namespace Victor { namespace Align2{

    /** @brief    Calculate structural scores with info derived from
     *                  threading.
     * 
     *   
     * @This 
     **/
    class Threading : public Structure {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        Threading(AlignmentData *ad, ThreadingInput *thread, double cThr);

        /// Copy constructor.
        Threading(const Threading &orig);

        /// Destructor.
        virtual ~Threading();


        // OPERATORS:

        /// Assignment operator.
        Threading& operator =(const Threading &orig);


        // PREDICATES:

        /// Calculate structural scores to create matrix values.
        virtual double scoringStr(int i, int j);


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const Threading &orig);

        /// Construct a new "deep copy" of this object.
        virtual Threading* newCopy();


    protected:


    private:

        // ATTRIBUTES:

        string seq1; ///< Target sequence.
        ThreadingInput *thread; ///< Template threading input file.
        double cThr; ///< Coefficient for threading.

    };

}} // namespace

#endif
