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
 

#ifndef __ProfInput_H__
#define __ProfInput_H__

#include <IoTools.h>
#include <stringtools.h>
#include <Substitution.h>
#include <iostream>
#include <string>

namespace Biopool
{
/** @brief  Implement I/O objects for handling PHD files.
 * 
* @Description  
* @This 
 **/
class ProfInput
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ProfInput();

	/// istream constructor.
	ProfInput(istream &is);

	/// Copy constructor.
	ProfInput(const ProfInput &orig);

	/// Destructor.
	virtual ~ProfInput();


// OPERATORS:

	/// Assignment operator.
	ProfInput& operator = (const ProfInput &orig);

	/// Output operator.
	friend ostream& operator << (ostream &os, const ProfInput &object);

	/// Input operator.
	friend istream& operator >> (istream &is, ProfInput &object);


// PREDICATES:

	/// Return the SSE prediction of Prof.
	char getProfSS(int i);

	/// Return the ASA prediction of Prof.
	char getProfBE(int i);

	/// Return the charachter for SSE/ASA combination.
	char getProfMixSSBE(char ssChar, char beChar);

	/// Return the size of the object referred as the number of aminoacids.
	virtual unsigned int size() const;


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ProfInput &orig);

	/// Construct a new "deep copy" of this object.
	virtual ProfInput* newCopy();


// HELPERS:

	/// Helper function used to write a string construct.
	static void pWriteString(ostream &os, string data1, string data2,
		string data3);

	/// Helper function used to read a string construct.
	static void pReadString(istream &is, string &data1, string &data2,
		string &data3);


protected:


private:

// ATTRIBUTES:

	string seq;           ///< Sequence.
	string profSSPred;    ///< SSE prediction.
	string profBEPred;    ///< ASA prediction.

};

// -----------------------------------------------------------------------------
//                                  ProfInput
// -----------------------------------------------------------------------------

// PREDICATES:

inline char
ProfInput::getProfSS(int i)
{
	return profSSPred[i];
}


inline char
ProfInput::getProfBE(int i)
{
	return profBEPred[i];
}


inline unsigned int
ProfInput::size() const
{
	return seq.size();
}

} // namespace

#endif
