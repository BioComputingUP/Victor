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
#ifndef _SIDECHAIN_CONSTRUCTOR_H_
#define _SIDECHAIN_CONSTRUCTOR_H_

// Includes:
#include <Spacer.h>
#include <SideChain.h>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/**@brief base reference stream
*@This 
 * */
class SideChainConstructor {
public: 

// CONSTRUCTORS/DESTRUCTOR:
  SideChainConstructor(istream& _refInput = cin) :
    refInput(_refInput), loaded(false), refSpacer() { }
  virtual ~SideChainConstructor() { PRINT_NAME; }  

// MODIFIERS:
  virtual SideChain& makeSideChain(string code); 

protected:
  // HELPERS:
  virtual void loadReference();

private:
  istream& refInput;             // reference stream
  int loaded;                    
  Spacer refSpacer;  // reference types
};

} // namespace
#endif //_SIDECHAIN_CONSTRUCTOR_H_
