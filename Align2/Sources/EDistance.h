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


#ifndef __EDistance_H__
#define __EDistance_H__

#include <Profile.h>
#include <ScoringFunction.h>

namespace Victor { namespace Align2{

    /** @brief Calculate scores for profile to profile alignment using
     *                  euclidean distance.
     * 
     *   
     * @This 
     **/
    class EDistance : public ScoringFunction {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        EDistance(Profile *pro1, Profile *pro2);

        /// Constructor assigning offset.
        EDistance(Profile *pro1, Profile *pro2, double offset);

        /// Copy constructor.
        EDistance(const EDistance &orig);

        /// Destructor.
        virtual ~EDistance();


        // OPERATORS:

        /// Assignment operator.
        EDistance& operator =(const EDistance &orig);


        // PREDICATES:

        /// Calculate scores to create matrix values.
        virtual double scoringSeq(int i, int j);

        /// Return offset.
        virtual double getOffset();


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const EDistance &orig);

        /// Construct a new "deep copy" of this object.
        virtual EDistance* newCopy();

        /// Set offset.
        virtual void setOffset(double off);


    protected:


    private:

        // ATTRIBUTES:

        Profile *pro1; ///< Target profile.
        Profile *pro2; ///< Template profile.
        double offset; ///< Offset.

    };

    // -----------------------------------------------------------------------------
    //                                  EDistance
    // -----------------------------------------------------------------------------

    // PREDICATES:

    inline double
    EDistance::getOffset() {
        return offset;
    }


    // MODIFIERS:

    inline void
    EDistance::setOffset(double off) {
        offset = off;
    }

}} // namespace

#endif
