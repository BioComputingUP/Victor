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
// Description:     Calculate score for profile to profile alignment with
//                  secondary structure. 
// 
// -----------------x-------------------x-------------------x-----------------

#include <ScoringP2Psec.h>
#include <math.h>

namespace Biopool
{

// CONSTRUCTORS:

ScoringP2Psec::ScoringP2Psec(Blosum* sub, Blosum* sec, 
 AlignmentData *a, Profile* p1, Profile* p2): ScoringScheme(sub,a),
 secstr(sec), seq1(a->getSequence(1)), seq2(a->getSequence(2)), 
 sec1(a->getSequence(3)), sec2(a->getSequence(4)), p1(p1), p2(p2)
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
  

ScoringP2Psec::ScoringP2Psec(const ScoringP2Psec& orig): ScoringScheme(orig), 
 p1(orig.p1), p2(orig.p2)
{
  copy(orig);
}


ScoringP2Psec::~ScoringP2Psec()
{ }


// OPERATORS:

ScoringP2Psec& ScoringP2Psec::operator = (const ScoringP2Psec& orig)
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
ScoringP2Psec::scoring(int i, int j)
{
  double s = 0.00;
  double tmp1 = 0.00; 
  double tmp2 = 0.00; 

  for (AminoAcidCode Amino = ALA; Amino <= TYR; Amino++)
    { 
      tmp1 = 0.0;
      for (AminoAcidCode Amino1 = ALA; Amino1 <= TYR; Amino1++)
	{  
	  tmp1 += (sub->score[aminoAcidOneLetterTranslator(
			Amino)][aminoAcidOneLetterTranslator(Amino1)])
	    * (p2->getAminoFrequencyFromCode(Amino1,(j-1))); 
	}
      tmp2 += (p1->getAminoFrequencyFromCode(Amino,(i-1)) * tmp1); 
    }
  //  s = log(tmp2+1);
  s = tmp2;
  s += secstr->score[sec1[i-1]][sec2[j-1]];
  
  return s;
}

// MODIFIERS:

void 
ScoringP2Psec::copy(const ScoringP2Psec& orig)
{	
  ScoringScheme::copy(orig);
  p1 = orig.p1->newCopy();
  p2 = orig.p2->newCopy();   
  secstr = orig.secstr;
  seq1 = orig.seq1;
  seq2 = orig.seq2;
  sec1 = orig.sec1;
  sec2 = orig.sec2;  
}


ScoringP2Psec* 
ScoringP2Psec::newCopy()
{
  ScoringP2Psec* tmp = new ScoringP2Psec(*this);
  return tmp;
}


void 
ScoringP2Psec::reverse()
{
  p2->reverse();

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



