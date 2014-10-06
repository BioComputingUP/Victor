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

#ifndef _XYZ_SAVER_H_
#define _XYZ_SAVER_H_

// Includes:
#include <Group.h>
#include <SideChain.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <Ligand.h>
#include <Saver.h>
#include <Debug.h>
#include <string>
#include <iostream>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/**@brief Saves components (Atoms, Groups, etc.) in carthesian format 
 * 
*@Description 
*    Internal format is made of type, coords & bonds of each atom,
*    one per line. Keywords "aminoacid" and "sidechain" delimit these
*    structures.
 * */
class XyzSaver : public Saver{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  XyzSaver(ostream& _output = cout, int _offset = 1) 
    : output(_output), delimit(true), offset(_offset) { }
  // this class uses the implicit copy operator.
  virtual ~XyzSaver() { PRINT_NAME; }  

// MODIFIERS:
  virtual void saveGroup(Group& gr);
  virtual void saveSideChain(SideChain& sc);
  virtual void saveAminoAcid(AminoAcid& aa);
  virtual void saveSpacer(Spacer& sp);
  virtual void saveLigand(Ligand& l);
  void setDelimit(bool _d) { delimit = _d; }

protected:
  virtual void pSaveAtomVector(vector<Atom>& va);

private:
  ostream& output;   // output stream
  bool delimit;      // write delimiters ("aminoacid", "sidechain", etc.)  
  int offset;        // ID offset for saving (optional) 
};
}
#endif //_XYZ_SAVER_H_
