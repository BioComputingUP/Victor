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
 
#ifndef __SeqDivergenceProfile_H__
#define __SeqDivergenceProfile_H__

#include <AGPFunction.h>
#include <NWAlign.h>
#include <Profile.h>
#include <ScoringS2S.h>
#include <SequenceData.h>
#include <SubMatrix.h>
#include <math.h>

namespace Biopool
{/** @brief    Calculate a frequency profile or PSSM using SeqDivergence
*                  weighting scheme.
 * 
* @Description  Some explanations can be found in:
*
*                  Guoli Wang, Roland L. Dunbrack jr.
*                  Scoring profile-to-profile sequence alignments.
*                  Institute for Cancer Research, Fox Chase Cancer Center,
*                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
*
* @This 
 **/

class SeqDivergenceProfile : public Profile
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	SeqDivergenceProfile();

	/// Copy constructor.
	SeqDivergenceProfile(const SeqDivergenceProfile &orig);

	/// Destructor.
	virtual ~SeqDivergenceProfile();


// OPERATORS:

	/// Assignment operator.
	SeqDivergenceProfile& operator = (const SeqDivergenceProfile &orig);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const SeqDivergenceProfile &orig);

	/// Construct a new "deep copy" of this object.
	virtual SeqDivergenceProfile* newCopy();


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

	vector<double> aliWeight;    ///< Alignment weights.

};

} // namespace

#endif
