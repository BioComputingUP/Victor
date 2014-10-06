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
 
#ifndef __HenikoffProfile_H__
#define __HenikoffProfile_H__

#include <Profile.h>

namespace Biopool
{
/** @brief Calculate a frequency profile or PSSM using Henikoff
*                  weighting scheme. 
 * 
* @Description  Some explanations can be found in:
*
*                  Guoli Wang, Roland L. Dunbrack jr.
*                  Scoring profile-to-profile sequence alignments.
*                  Institute for Cancer Research, Fox Chase Cancer Center,
*                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
* @This 
 **/
class HenikoffProfile : public Profile
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	HenikoffProfile();

	/// Copy constructor.
	HenikoffProfile(const HenikoffProfile &orig);

	/// Destructor.
	virtual ~HenikoffProfile();


// OPERATORS:

	/// Assignment operator.
	HenikoffProfile& operator = (const HenikoffProfile &orig);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const HenikoffProfile &orig);

	/// Construct a new "deep copy" of this object.
	virtual HenikoffProfile* newCopy();


// HELPERS:

	/// Calculate alignment weights.
	//to save computational time we suggest cLen=25. Francesco Lovo 2012
	virtual void pCalculateWeight(Alignment &ali,unsigned int cLen=50);

	/// Calculate the raw (ie. unnormalized) aminoacids frequencies for position i.
	virtual void pCalculateRawFrequency(vector<double> &freq, double &gapFreq,
		Alignment &ali, unsigned int i);

	/// Construct data from alignment.
	virtual void pConstructData(Alignment &ali);


protected:


private:

	vector< vector<double> > aliWeight;    ///< Alignment weights.

};

} // namespace

#endif
