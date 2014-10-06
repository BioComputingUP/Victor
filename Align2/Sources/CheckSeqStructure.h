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
// Description:     This class is used to control the correspondence between
//                  a Fasta sequence file and a Pdb structure file. 
//                  It assign a progressive number to the aa of the sequence
//                  according to the aa of the structures.
// 
// -----------------x-------------------x-------------------x-----------------

#ifndef __CheckSeqStructure_H__
#define __CheckSeqStructure_H__

// Includes

#include <Spacer.h>
#include <iostream>
#include <string>

namespace Biopool
{
 
    /** @brief  This class is used to control the correspondence between a Fasta 
* sequence file and a Pdb structure file. 
 * 
* @Description  It assign a progressive number 
* to the aa of the sequence according to the aa of the structures.
* @This 
 **/
    class CheckSeqStructure 
{
public:
  
// CONSTRUCTORS:
  /// Default constructor.
  CheckSeqStructure(const Spacer& sp, const string& seq);
  /// Copy constructor.
  CheckSeqStructure(const CheckSeqStructure& orig);
  /// Destructor.
  virtual ~CheckSeqStructure();
  
// OPERATORS: 
  /// Assigment operator.
  CheckSeqStructure& operator = (const CheckSeqStructure& orig);
  
// PREDICATES:
  // Returns sequence.
  string getSequence();
  // Returns PDB number corresponding to sequence number.
  int getPdbSeqNumber(unsigned int i);
  
protected:
  
// MODIFIERS:
  /// Copies orig object to this object ("deep copy").
  virtual void copy(const CheckSeqStructure& orig);
  /// Defines start offset between seq and sp.
  void SetStartOffset( string& seq, Spacer& sp );
  /// Defines end offset between seq and sp.
  void SetEndOffset( string& seq, Spacer& sp );
  void setPdbSeqNumber(int index, int pbdnumber);
  /// Sets confidence threshold in AA for "identity".
  void setConfidence(int n);
  /// Generates equivalence array.
  void pConstructArray();
  /// Special case in equivalence array generation.
  void pSpecialCase();
  
// ATTRIBUTES:
  
  Spacer *sp;                ///< PDB Structure.
  string seq;               ///< Sequence to check.
  int endoffset[2];         ///< Array of end offsets for both sequences.
  int startoffset[2];       ///< Array of start offsets for both sequences.
  vector<int> pdbseqnumber; 
  ///< Correspondence PDB numbers from the structure assigned to the sequence.
  unsigned int conf;        ///< Confidence threshold in AA for "identity".
  
private:
};

// -----------------x-------------------x-------------------x-----------------
//                               CheckSeqStructure
// -----------------x-------------------x-------------------x-----------------

// PREDICATES:

inline string
CheckSeqStructure::getSequence()
{
  return seq;
}


inline int
CheckSeqStructure::getPdbSeqNumber(unsigned int i)
{
  return pdbseqnumber[i];
}


inline void 
CheckSeqStructure::setPdbSeqNumber(int index, int pdbnumber)
{
  pdbseqnumber[index] = pdbnumber;
}


inline void
CheckSeqStructure::setConfidence(int n)
{
  conf = n;
}
  


}
#endif /* __CheckSeqStructure_H__ */
  




