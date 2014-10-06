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
// Description:     Implement I/O objects for handling BLAST PSSM (Position
//                  Specific Score Matrix).
//
// -----------------x-----------------------------------------------------------

#include <PssmInput.h>

namespace Biopool
{

// CONSTRUCTORS:

PssmInput::PssmInput()
{ }


PssmInput::PssmInput(istream &is)
{
	is >> *this;
}


PssmInput::PssmInput(const PssmInput &orig)
{
	PssmInput::copy(orig);
}


PssmInput::~PssmInput()
{ }


// OPERATORS:

PssmInput&
PssmInput::operator = (const PssmInput &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


ostream&
operator << (ostream &os, const PssmInput &object)
{
	PssmInput::pWriteDoubleVector(os, object.residuescores, object.allPosition,
		object.allAa);
	return os;
}


istream&
operator >> (istream &is, PssmInput &object)
{
	PssmInput::pReadDoubleVector(is, object.residuescores, object.allPosition,
		object.allAa);
	return is;
}


// MODIFIERS:

void
PssmInput::copy(const PssmInput &orig)
{
	residuescores.clear();
	residuescores.reserve(orig.residuescores.size());
	allPosition.clear();
	allPosition.reserve(orig.allPosition.size());
	allAa.clear();
	allAa.reserve(orig.allAa.size());
	vector<string> tmp1, tmp2;
	tmp1.reserve(orig.allPosition.size());
	tmp2.reserve(orig.allAa.size());

	for (unsigned int i = 0; i < orig.residuescores.size(); i++)
	{
		vector<double> tmp;
		tmp.reserve(orig.residuescores[i].size());
		tmp1.push_back(orig.allPosition[i]);
		tmp2.push_back(orig.allAa[i]);
		for (unsigned int j = 0; j < orig.residuescores[i].size(); j++)
			tmp.push_back(orig.residuescores[i][j]);
		residuescores.push_back(tmp);
	}
}


PssmInput*
PssmInput::newCopy()
{
	PssmInput *tmp = new PssmInput(*this);
	return tmp;
}


// HELPERS:

template<class T> void
PssmInput::pWriteDoubleVector(ostream &os, vector< vector<T> > data,
	vector<string> data1, vector<string> data2)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		os << data1[i] << " "
		   << data2[i] << " ";
		for (unsigned int j = 0; j < 20; j++)
			os << data[i][j] << " ";
		os << endl;
	}
	os << endl;
}


template<class T> void
PssmInput::pReadDoubleVector(istream &is, vector< vector<T> > &data,
	vector<string> &data1, vector<string> &data2)
{
	for (unsigned int h = 0; h < 25; h++) // 5 words + 20 aa char
	{
		string Tmp;
		is >> Tmp;
	}

	while (is)
	{
		string position;
		is >> position;
		if (position == "K")
			break;

		data1.push_back(position);
		string aa;
		is >> aa;

		data2.push_back(aa);
		vector<T> row;
		row.reserve(20);

		for (unsigned int j = 0; j < 20; j++)
		{
			T tmp;
			is >> tmp;
			row.push_back(tmp);
		}

		data.push_back(row);
	}
}


template void
PssmInput::pReadDoubleVector(istream &is, vector< vector<int> > &data,
	vector<string> &data1, vector<string> &data2);


template void
PssmInput::pReadDoubleVector(istream &is, vector< vector<double> > &data,
	vector<string> &data1, vector<string> &data2);

} // namespace
