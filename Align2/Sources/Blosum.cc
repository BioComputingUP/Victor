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
// --*- C++ -*------x---------------------------------------------------------
//
// Description:     This class implemenents a standard substitution matrix.
// 
// -----------------x-------------------x-------------------x-----------------

#include <Blosum.h>

// CONSTRUCTORS:

Blosum::Blosum() : Substitution()
{ }


Blosum::Blosum(const Blosum& orig) : Substitution(orig),
    residuescores(orig.residuescores), residues(orig.residues)
{
  copy(orig);
}


Blosum::Blosum(istream& is) : Substitution()
{
  is >> *this;
}


Blosum::~Blosum()
{ }


// OPERATORS:

Blosum& Blosum::operator = (const Blosum& orig)
{
  if (&orig != this)
    copy(orig);
  return *this;
}


ostream& 
operator << (ostream& os, const Blosum& object)
{
  os << object.residues << endl;
  Substitution::pWriteDoubleVector(os, object.residuescores);
  return os;
}


istream& 
operator >> (istream& is, Blosum& object)
{
  is >> object.residues; 
  Substitution::pReadDoubleVector(is, object.residuescores);
  object.buildscore(object.residues, object.residuescores); 
  return is;
}


// MODIFIERS:

void 
Blosum::copy(const Blosum& orig)
{
  Substitution::copy(orig); // copy score matrix
  residuescores = orig.residuescores;
  residues = orig.residues;
}




