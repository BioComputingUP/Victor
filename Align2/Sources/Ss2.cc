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
// Description:     Calculate structural scores with info derived from PSI-PRED.
//
// -----------------x-----------------------------------------------------------

#include <Ss2.h>

namespace Biopool
{

// CONSTRUCTORS:

Ss2::Ss2(SubMatrix *subStr, AlignmentData *ad, Ss2Input *psipred1,
	Ss2Input *psipred2, double cSs2) : Structure(subStr),
		sec1(ad->getSequence(3)), sec2(ad->getSequence(4)), psipred1(psipred1),
			psipred2(psipred2), cSs2(cSs2)
{ }


Ss2::Ss2(const Ss2 &orig) : Structure(orig)
{
	copy(orig);
}


Ss2::~Ss2()
{ }


// OPERATORS:

Ss2&
Ss2::operator = (const Ss2 &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
Ss2::scoringStr(int i, int j)
{
	const string ss2_indices = "HEC";


	// Target PSI-PRED input file

	double weigthH = psipred1->score((i-1), 1);
	double weigthE = psipred1->score((i-1), 2);
	double weigthC = psipred1->score((i-1), 0);
	int substiH = subStr->score[sec2[j-1]][ss2_indices[0]];
	int substiE = subStr->score[sec2[j-1]][ss2_indices[1]];
	int substiC = subStr->score[sec2[j-1]][ss2_indices[2]];
	double tmp1 = weigthH * substiH + weigthE * substiE + weigthC * substiC;


	// Template PSI-PRED input file

	weigthH = psipred2->score((j-1), 1);
	weigthE = psipred2->score((j-1), 2);
	weigthC = psipred2->score((j-1), 0);
	substiH = subStr->score[sec1[i-1]][ss2_indices[0]];
	substiE = subStr->score[sec1[i-1]][ss2_indices[1]];
	substiC = subStr->score[sec1[i-1]][ss2_indices[2]];
	double tmp2 = weigthH * substiH + weigthE * substiE + weigthC * substiC;


	return cSs2 * ((tmp1 + tmp2) / 2);
}


// MODIFIERS:

void
Ss2::copy(const Ss2 &orig)
{
	Structure::copy(orig);
	sec1 = orig.sec1;
	sec2 = orig.sec2;
	psipred1 = orig.psipred1->newCopy();
	psipred2 = orig.psipred2->newCopy();
	cSs2 = orig.cSs2;
}


Ss2*
Ss2::newCopy()
{
	Ss2 *tmp = new Ss2(*this);
	return tmp;
}


void
Ss2::reverse()
{
	string tmp = "";
	for (unsigned int i = sec2.length(); i > 0; i--)
		tmp.push_back(sec2[i-1]);
	sec2 = tmp;
}

} // namespace
