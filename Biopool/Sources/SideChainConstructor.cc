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
 * */

// Includes:
#include <SideChainConstructor.h>
#include <PdbLoader.h>
#include <vector3.h>

// Global constants, typedefs, etc. (to avoid):

using namespace Biopool;


// CONSTRUCTORS/DESTRUCTOR:


// MODIFIERS:

void SideChainConstructor::loadReference(){
  PdbLoader pl(refInput);
  refSpacer.load(pl);
  loaded = true;
}

SideChain& SideChainConstructor::makeSideChain(string type){
  if (!loaded)
    loadReference();
  for (unsigned int i = 0; i < refSpacer.sizeAmino(); i++)
    if (refSpacer.getAmino(i).getSideChain().getType() == type)      {
	SideChain& sc = refSpacer.getAmino(i).getSideChain();
	return sc;
      }
  ERROR("No reference structure found.", exception);
  SideChain* sc = new SideChain;
  return *sc;
}

