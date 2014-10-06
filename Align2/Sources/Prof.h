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
 

#ifndef __Prof_H__
#define __Prof_H__

#include <ProfInput.h>
#include <Structure.h>

namespace Biopool
{
/** @brief   Calculate structural scores with info derived from PHD.
 * 
* @Description  
* @This 
 **/
class Prof : public Structure
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	Prof(SubMatrix *subStr, ProfInput *phd1, ProfInput *phd2, double cPrf);

	/// Copy constructor.
	Prof(const Prof &orig);

	/// Destructor.
	virtual ~Prof();


// OPERATORS:

	/// Assignment operator.
	Prof& operator = (const Prof &orig);


// PREDICATES:

	/// Calculate structural scores to create matrix values.
	virtual double scoringStr(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const Prof &orig);

	/// Construct a new "deep copy" of this object.
	virtual Prof* newCopy();


protected:


private:

// ATTRIBUTES:

	ProfInput *phd1;    ///< Target PHD input file.
	ProfInput *phd2;    ///< Template PHD input file.
	double cPrf;        ///< Coefficient for PHD prediction.

};

} // namespace

#endif
