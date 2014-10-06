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
 
//
// Description:     Print alignment of two sequences.
//
// -----------------x-----------------------------------------------------------

#include <SequenceData.h>

namespace Biopool
{

// CONSTRUCTORS:

SequenceData::SequenceData(int n, const string &_seq1, const string &_seq2,
	const string &_n1, const string &_n2) : AlignmentData(n, _n1, _n2),
		seq1(_seq1), seq2(_seq2)
{ }


SequenceData::SequenceData(const SequenceData &orig) : AlignmentData(orig)
{
	copy(orig);
}


SequenceData::~SequenceData()
{ }


// OPERATORS:

SequenceData&
SequenceData::operator = (const SequenceData &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

string
SequenceData::getSequence(int n)
{
	if (n == 1)
		return seq1;
	else
		return seq2;
}


void
SequenceData::calculateMatch(int i, int tbi, int j, int tbj)
{
	string res1 = "";
	string res2 = "";

	if (i == tbi)
		res1 = res1 + "-";
	else
		res1 = res1 + seq1[i-1];

	if (j == tbj)
		res2 = res2 + "-";
	else
		res2 = res2 + seq2[j-1];

	add(res1, 0);
	add(res2, 1);
}


void
SequenceData::getMatch()
{
	for (int i = 0; i < n; i++)
	{
		string tmp = match[i];
		reverse(tmp.begin(), tmp.end());
		match[i] = tmp;
	}
}


void
SequenceData::outputMatch(ostream &os, bool fasta)
{
	string temp = "";
	unsigned int cons = 0;
	unsigned int sim = 0;

	for (unsigned int j = 0; j < match[0].length(); j++)
		if (match[0][j] == match[1][j])
		{
			temp += match[0][j];
			cons++;
		}
		else
			if (similar(match[0][j], match[1][j]))
			{
				temp += '+';
				sim++;
			}
			else
				temp += " ";

	if (!fasta)
	{
		os << "Conservation = " << cons << " / " << cons + sim << " / "
		   << match[0].length() << "\t\t" << setw(5) << setprecision(3)
		   << (static_cast<double>(cons) / match[0].length() * 100) << "%  /  "
		   << setw(5) << setprecision(3)
		   << (static_cast<double>(cons + sim) / match[0].length() * 100) << "%\n\n"
		   << match[0] << "\n"
		   << temp << "\n"
		   << match[1] << "\n"
		   << "--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--xx--" << endl;
	}
	else
	{
		AlignmentBase::saveFasta(match[0], name1, os);
		AlignmentBase::saveFasta(match[1], name2, os);
	}

	clear();
}


Alignment&
SequenceData::generateMatch(double score)
{
	Alignment *ali = new Alignment;
	ali->setTarget(match[0], name1);
	ali->setTemplate(match[1], name2, score);
	clear();
	return *ali;
}


// MODIFIERS:

void
SequenceData::copy(const SequenceData &orig)
{
	AlignmentData::copy(orig);
	seq1 = orig.seq1;
	seq2 = orig.seq2;
}


SequenceData*
SequenceData::newCopy()
{
	SequenceData *tmp = new SequenceData(*this);
	return tmp;
}


void
SequenceData::setSequence(string s, int n)
{
	if (n == 1)
		seq1 = s;
	else
		seq2 = s;
}

} // namespace
