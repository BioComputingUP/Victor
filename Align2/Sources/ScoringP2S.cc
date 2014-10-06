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
// Description:     Calculate scores for profile to sequence alignment.
//
// -----------------x-----------------------------------------------------------

#include <ScoringP2S.h>

namespace Biopool
{

// CONSTRUCTORS:

ScoringP2S::ScoringP2S(SubMatrix *sub, AlignmentData *ad, Structure *str,
	Profile *pro, double cSeq) : ScoringScheme(sub, ad, str),
		seq1(ad->getSequence(1)), seq2(ad->getSequence(2)), pro(pro), cSeq(cSeq)
{ }


ScoringP2S::ScoringP2S(const ScoringP2S &orig) : ScoringScheme(orig)
{
	copy(orig);
}


ScoringP2S::~ScoringP2S()
{ }


// OPERATORS:

ScoringP2S&
ScoringP2S::operator = (const ScoringP2S &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
ScoringP2S::scoring(int i, int j)
{
	double s = 0.00;

	for (AminoAcidCode amino = ALA; amino < TYR; amino++)
	{
		char aminoacid = aminoAcidOneLetterTranslator(amino);
		s += (sub->score[seq2[j-1]][aminoacid]) *
			(pro->getAminoFrequencyFromCode(amino, (i-1)));
	}
	s *= cSeq;

	if (str != 0)
		s += str->scoringStr(i, j);

	return s;
}


// MODIFIERS:

void
ScoringP2S::copy(const ScoringP2S &orig)
{
	ScoringScheme::copy(orig);
	seq1 = orig.seq1;
	seq2 = orig.seq2;
	pro = orig.pro->newCopy();
	cSeq = orig.cSeq;
}


ScoringP2S*
ScoringP2S::newCopy()
{
	ScoringP2S *tmp = new ScoringP2S(*this);
	return tmp;
}


void
ScoringP2S::reverse()
{
	ScoringScheme::reverse();

	string tmp = "";
	for (unsigned int i = seq2.length(); i > 0; i--)
		tmp.push_back(seq2[i-1]);
	seq2 = tmp;
}

} // namespace
