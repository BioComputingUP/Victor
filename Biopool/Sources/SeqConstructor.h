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
/**
 *@Description:
*    This class builds a spacer by concatenating the same aminoacid type for 
*    n times.
*/
#ifndef _SIDECHAIN_CONSTRUCTOR_H_
#define _SIDECHAIN_CONSTRUCTOR_H_

// Includes:
#include <Spacer.h>
#include <AminoAcid.h>
#include <Debug.h>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/**@brief This class builds a spacer by concatenating the same aminoacid type for 
*    n times.
 * 
*@Description  
*@This 
 * */
class SeqConstructor{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  SeqConstructor(istream& _refInput = cin) :
    refInput(_refInput), loaded(false), refAmino() { }
  virtual ~SeqConstructor() { PRINT_NAME; }  

// MODIFIERS:
  virtual Spacer& makeSpacer(string type, unsigned int n); 
  // make spacer by concatenating code n-times 

protected:
  // HELPERS:
  virtual void loadReference();
  void searchReference(AminoAcid& aa, string type);
  void buildAminoAcid(AminoAcid& aa, AminoAcid* prev);

private:
  istream& refInput;             // reference stream
  bool loaded;                   // is reference loaded yet?
  vector<AminoAcid> refAmino;    // reference types
};
 
} // namespace
#endif //_SEQ_CONSTRUCTOR_H_
