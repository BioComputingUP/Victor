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
 

#ifndef __ScoringS2S_H__
#define __ScoringS2S_H__

#include <ScoringScheme.h>

namespace Biopool
{
/** @brief  Calculate scores for sequence to sequence alignment.
 * 
* @Description  
* @This 
 **/
class ScoringS2S : public ScoringScheme
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ScoringS2S(SubMatrix *sub, AlignmentData *ad, Structure *str, double cSeq);

	/// Copy constructor.
	ScoringS2S(const ScoringS2S &orig);

	/// Destructor.
	virtual ~ScoringS2S();


// OPERATORS:

	/// Assignment operator.
	ScoringS2S& operator = (const ScoringS2S &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoring(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ScoringS2S &orig);

	/// Construct a new "deep copy" of this object.
	virtual ScoringS2S* newCopy();

	/// Reverse template sequence.
	virtual void reverse();


protected:


private:

// ATTRIBUTES:

	string seq1;    ///< Target sequence.
	string seq2;    ///< Template sequence.
	double cSeq;    ///< Coefficient for sequence alignment.

};

} // namespace

#endif
