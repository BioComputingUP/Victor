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
// Description:     Calculate score for sequence to sequence alignment 
//                  considering also the secondary structure.
// 
// -----------------x-------------------x-------------------x-----------------

#include <ScoringS2Ssec.h>
namespace Biopool
{

// CONSTRUCTORS:

ScoringS2Ssec::ScoringS2Ssec(Blosum* sub, Blosum* sec, 
AlignmentData *a): ScoringScheme(sub,a), secstr(sec), seq1(a->getSequence(1)), 
seq2(a->getSequence(2)),sec1(a->getSequence(3)), sec2(a->getSequence(4))
{	
  //check secondary sctructure
  if ( (sec1.size() != seq1.size()) || (sec2.size() != seq2.size()) 
       || (!checkSequence(sec1)) || (!checkSequence(sec2)) ) 
    {
      cout << "Illegal sequence: " << sec1 
	   << " " << sec2 << " " << sub->getResidues();
      ERROR("Error checking sequence.", exception);
    }
  
}


ScoringS2Ssec::ScoringS2Ssec(const ScoringS2Ssec& orig) : ScoringScheme(orig)
{
  copy(orig);
}


ScoringS2Ssec::~ScoringS2Ssec()
{ }


// OPERATORS:

ScoringS2Ssec& ScoringS2Ssec::operator = (const ScoringS2Ssec& orig)
{
  if (&orig != this)
  {
    copy(orig);
  }
//  POSTCOND( (orig == *this), exception);
  return *this;
}


// PREDICATES:

double
ScoringS2Ssec::scoring(int i, int j)
{
  double s = sub->score[seq1[i-1]][seq2[j-1]] 
    + secstr->score[sec1[i-1]][sec2[j-1]];
  return s;
}


// MODIFIERS:

void 
ScoringS2Ssec::copy(const ScoringS2Ssec& orig)
{
  ScoringScheme::copy(orig);
  secstr = orig.secstr;
  seq1 = orig.seq1;
  seq2 = orig.seq2;
  sec1 = orig.sec1;
  sec2 = orig.sec2;
}


ScoringS2Ssec* 
ScoringS2Ssec::newCopy()
{
  ScoringS2Ssec* tmp = new ScoringS2Ssec(*this);
  return tmp;
}


void 
ScoringS2Ssec::reverse()
{
  string tmp = "";
  for (unsigned int i = seq2.length(); i > 0; i--)
    tmp.push_back(seq2[i-1]);
  seq2 = tmp;

  string tmpS = "";
  for (unsigned int i = sec2.length(); i > 0; i--)
    tmpS.push_back(sec2[i-1]);
  sec2 = tmpS;
}


} // namespace



