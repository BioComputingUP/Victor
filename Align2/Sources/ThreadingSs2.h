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
 

#ifndef __ThreadingSs2_H__
#define __ThreadingSs2_H__

#include <AlignmentData.h>
#include <Ss2Input.h>
#include <Structure.h>
#include <ThreadingInput.h>

namespace Biopool
{
/** @brief   Calculate structural scores with info derived from
*                  threading and PSI-PRED.
 * 
* @Description  
* @This 
 **/
class ThreadingSs2 : public Structure
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ThreadingSs2(SubMatrix *subStr, AlignmentData *ad, ThreadingInput *thread,
		Ss2Input *psipred1, Ss2Input *psipred2, double cThr, double cSs2);

	/// Copy constructor.
	ThreadingSs2(const ThreadingSs2 &orig);

	/// Destructor.
	virtual ~ThreadingSs2();


// OPERATORS:

	/// Assignment operator.
	ThreadingSs2& operator = (const ThreadingSs2 &orig);


// PREDICATES:

	/// Calculate structural scores to create matrix values.
	virtual double scoringStr(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ThreadingSs2 &orig);

	/// Construct a new "deep copy" of this object.
	virtual ThreadingSs2* newCopy();


protected:


private:

// ATTRIBUTES:

	string seq1;               ///< Target sequence.
	string sec1;               ///< Target secondary structure.
	string sec2;               ///< Template secondary structure.
	ThreadingInput *thread;    ///< Template threading input file.
	Ss2Input *psipred1;        ///< Target PSI-PRED input file.
	Ss2Input *psipred2;        ///< Template PSI-PRED input file.
	double cThr;               ///< Coefficient for threading.
	double cSs2;               ///< Coefficient for PSI-PRED prediction.

};

} // namespace

#endif
