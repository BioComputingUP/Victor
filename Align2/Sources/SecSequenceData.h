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
 

#ifndef __SecSequenceData_H__
#define __SecSequenceData_H__

#include <AlignmentData.h>

namespace Biopool
{
/** @brief   Print alignment of two sequences considering also secondary
*                  structure.
 * 
* @Description  
* @This 
 **/
class SecSequenceData : public AlignmentData
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	SecSequenceData(int n, const string &seq1, const string &seq2,
		const string &sec1, const string &sec2, const string &name1 = "Seq1",
			const string &name2 = "Seq2");

	/// Copy constructor.
	SecSequenceData(const SecSequenceData &orig);

	/// Destructor.
	virtual ~SecSequenceData();


// OPERATORS:

	/// Assignment operator.
	SecSequenceData& operator = (const SecSequenceData &orig);


// PREDICATES:

	/// Return the sequence at position n of the vector.
	virtual string getSequence(int n);

	/// Calculate single match positions.
	virtual void calculateMatch(int i, int tbi, int j, int tbj);

	/// Reverse the strings of the vector.
	virtual void getMatch();

	/// Control if the strings of the vector are similar and print them.
	virtual void outputMatch(ostream &os, bool fasta = false);

	/// Generate and return an ensemble of suboptimal alignments.
	virtual Alignment& generateMatch(double score = 0.00);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const SecSequenceData &orig);

	/// Construct a new "deep copy" of this object.
	virtual SecSequenceData* newCopy();

	/// Set the sequence at position n of the vector.
	virtual void setSequence(string s, int n);


protected:


private:

// ATTRIBUTES:

	string seq1;    ///< Target sequence.
	string seq2;    ///< Template sequence.
	string sec1;    ///< Target secondary structure.
	string sec2;    ///< Template secondary structure.

};

} // namespace

#endif
