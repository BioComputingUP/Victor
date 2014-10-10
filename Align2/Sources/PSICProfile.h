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


#ifndef __PSICProfile_H__
#define __PSICProfile_H__

#include <Profile.h>

namespace Victor { namespace Align2{

    /** @brief   Calculate a frequency profile or PSSM using PSIC weighting
     *                  scheme.
     * 
     * @Description  Some explanations can be found in:
     *
     *                  Guoli Wang, Roland L. Dunbrack jr.
     *                  Scoring profile-to-profile sequence alignments.
     *                  Institute for Cancer Research, Fox Chase Cancer Center,
     *                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
     * @This 
     **/
    class PSICProfile : public Profile {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        PSICProfile();

        /// Copy constructor.
        PSICProfile(const PSICProfile &orig);

        /// Destructor.
        virtual ~PSICProfile();


        // OPERATORS:

        /// Assignment operator.
        PSICProfile& operator =(const PSICProfile &orig);


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const PSICProfile &orig);

        /// Construct a new "deep copy" of this object.
        virtual PSICProfile* newCopy();


        // HELPERS:

        /// Calculate alignment weights.
        virtual void pCalculateWeight(Alignment &ali);

        /// Calculate the raw (ie. unnormalized) aminoacids frequencies for position i.
        virtual void pCalculateRawFrequency(vector<double> &freq, double &gapFreq,
                Alignment &ali, unsigned int i);

        /// Construct data from alignment.
        virtual void pConstructData(Alignment &ali);


    protected:


    private:

        vector< vector<double> > aliWeight; ///< Alignment weights.

    };

}} // namespace

#endif
