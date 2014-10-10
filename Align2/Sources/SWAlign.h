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


#ifndef __SWAlign_H__
#define __SWAlign_H__

#include <Align.h>

namespace Victor { namespace Align2{

    /** @brief  Implement Smith-Waterman local alignment.
     * 
     *   

     **/
    class SWAlign : public Align {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        SWAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss);

        /// Constructor with weighted alignment positions.
        SWAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss,
                const vector<unsigned int> &v1, const vector<unsigned int> &v2);

        /// Copy constructor.
        SWAlign(const SWAlign &orig);

        /// Destructor.
        virtual ~SWAlign();


        // OPERATORS:

        /// Assignment operator.
        SWAlign& operator =(const SWAlign &orig);


        // PREDICATES:

        /// Return two-element array containing an alignment with maximal score.
        virtual void getMultiMatch();


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const SWAlign &orig);

        /// Construct a new "deep copy" of this object.
        virtual SWAlign* newCopy();


        // HELPERS:

        /// Update/create matrix values.
        virtual void pCalculateMatrix(bool update = true);

        /// Update/create weighted matrix values.
        virtual void pCalculateMatrix(const vector<unsigned int> &v1,
                const vector<unsigned int> &v2, bool update = true);


    protected:


    private:

    };

}} // namespace

#endif
