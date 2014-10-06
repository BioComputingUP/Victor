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
 

#ifndef __ScoringFunction_H__
#define __ScoringFunction_H__

#include <math.h>
#include <string>

namespace Biopool
{
/** @brief    Base class for scoring functions.
 * 
* @Description  
* @This 
 **/
class ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ScoringFunction()
	{ }

	/// Copy constructor.
	ScoringFunction(const ScoringFunction &orig)
	{
		copy(orig);
	}

	/// Destructor.
	virtual ~ScoringFunction()
	{ }


// OPERATORS:

	/// Assignment operator.
	ScoringFunction& operator = (const ScoringFunction &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j) = 0;


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ScoringFunction &orig);

	/// Construct a new "deep copy" of this object.
	virtual ScoringFunction* newCopy() = 0;


protected:


private:

};

// -----------------------------------------------------------------------------
//                               ScoringFunction
// -----------------------------------------------------------------------------

// OPERATORS:

inline ScoringFunction&
ScoringFunction::operator = (const ScoringFunction &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// MODIFIERS:

inline void
ScoringFunction::copy(const ScoringFunction &orig)
{ }

} // namespace

#endif
