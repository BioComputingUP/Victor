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


#ifndef __Panchenko_H__
#define __Panchenko_H__

#include <Profile.h>
#include <PssmInput.h>
#include <ScoringFunction.h>

namespace Victor { namespace Align2{

    /** @brief  Calculate scores for profile to profile alignment using
     *                  Panchenko method. 
     * 
     *   Some explanations can be found in:
     *
     *                  Panchenko AR.
     *                  Finding weak similarities between proteins by sequence
     *                  profile comparison.
     *                  Nucleic Acids Res. 2003 Jan 15;31(2):683-9.
     * @This 
     **/
    class Panchenko : public ScoringFunction {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        Panchenko(Profile *pro1, Profile *pro2, PssmInput *pssm1, PssmInput *pssm2);

        /// Copy constructor.
        Panchenko(const Panchenko &orig);

        /// Destructor.
        virtual ~Panchenko();


        // OPERATORS:

        /// Assignment operator.
        Panchenko& operator =(const Panchenko &orig);


        // PREDICATES:

        /// Calculate scores to create matrix values.
        virtual double scoringSeq(int i, int j);

        /// Return the number of different aminoacids in column i.
        int returnAaColumnTarget(int i);

        /// Return the number of different aminoacids in column i.
        int returnAaColumnTemplate(int i);


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const Panchenko &orig);

        /// Construct a new "deep copy" of this object.
        virtual Panchenko* newCopy();


    protected:


    private:

        // ATTRIBUTES:

        Profile *pro1; ///< Target profile.
        Profile *pro2; ///< Template profile.
        PssmInput *pssm1; ///< Target PSSM.
        PssmInput *pssm2; ///< Template PSSM.

    };

}} // namespace

#endif
