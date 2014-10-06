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
 

#ifndef __ScoringP2Psec_H__
#define __ScoringP2Psec_H__

#include <ScoringScheme.h>
#include <Profile.h>
#include <string>


namespace Biopool
{/** @brief  Calculate score for profile to profile alignment with
*                  secondary structure. 
 * 
* @Description  
* @This 
 **/

class ScoringP2Psec : public ScoringScheme
{
public:

// CONSTRUCTORS:
  /// Default constructor.
  ScoringP2Psec(Blosum* sub, Blosum* sec, AlignmentData *a, 
		Profile* p1, Profile* p2);
  /// Copy constructor.
  ScoringP2Psec(const ScoringP2Psec& orig);
  /// Destructor.
  virtual ~ScoringP2Psec();
  
// OPERATORS:
  /// Assignment operator.
  ScoringP2Psec& operator = (const ScoringP2Psec& orig);
  
// PREDICATES:
  /// Calculate scores to create Matrix values.
  double scoring(int i, int j);
  
// MODIFIERS:
  /// Copies orig object to this object ("deep copy").
  virtual void copy(const ScoringP2Psec& orig);
  /// Constructs a new ("deep") copy of this object.
  virtual ScoringP2Psec* newCopy();
  /// Reverse second component (sequence or profile)
  virtual void reverse();

protected:
  
private:
  
// ATTRIBUTES:
  Blosum* secstr;        ///< Secondary structure substitution matrix. 
  string seq1;          ///< First sequence.
  string seq2;          ///< Second sequence.
  string sec1;	        ///< First secondary structure.
  string sec2;	        ///< Second secondary structure.
  Profile* p1;           ///< First profile.
  Profile* p2;		///< Second profile.
};

}

#endif
