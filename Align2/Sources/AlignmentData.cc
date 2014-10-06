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
// --*- C++ -*------x-----------------------------------------------------------
//
// Description:     Base class for printing alignments.
//
// -----------------x-----------------------------------------------------------

#include <AlignmentData.h>

namespace Biopool
{

// CONSTRUCTORS:

AlignmentData::AlignmentData(int n, const string &_n1, const string &_n2)
	: match(n), n(n), name1(_n1), name2(_n2)
{ }


AlignmentData::AlignmentData(const AlignmentData &orig)
{
	copy(orig);
}


AlignmentData::~AlignmentData()
{ }


// OPERATORS:

AlignmentData&
AlignmentData::operator = (const AlignmentData &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

bool
AlignmentData::similar(char a, char b)
{
	// polar
	if ((a == 'S') && (b == 'T'))
		return true;
	if ((a == 'T') && (b == 'S'))
		return true;

	// polar & positive
	if ((a == 'K') && (b == 'R'))
		return true;
	if ((a == 'R') && (b == 'K'))
		return true;

	// polar & negative
	if ((a == 'D') && ((b == 'E') || (b == 'N')))
		return true;
	if ((a == 'E') && ((b == 'D') || (b == 'Q')))
		return true;
	if ((a == 'Q') && ((b == 'E') || (b == 'N')))
		return true;
	if ((a == 'N') && ((b == 'D') || (b == 'Q')))
		return true;

	// non-polar
	if ((a == 'L') && ((b == 'I') || (b == 'V') || (b == 'A')))
		return true;
	if ((a == 'I') && ((b == 'L') || (b == 'V') || (b == 'A')))
		return true;
	if ((a == 'V') && ((b == 'L') || (b == 'I') || (b == 'A')))
		return true;
	if ((a == 'A') && ((b == 'L') || (b == 'I') || (b == 'V')))
		return true;

	// other
	if ((a == 'F') && (b == 'Y'))
		return true;
	if ((a == 'Y') && (b == 'F'))
		return true;

	return false;
}


// MODIFIERS:

void
AlignmentData::copy(const AlignmentData &orig)
{
	match.clear();
	for (unsigned int i = 0; i < orig.match.size(); i++)
		match.push_back(orig.match[i]);

	n = orig.n;
	name1 = orig.name1;
	name2 = orig.name2;
}

} // namespace
