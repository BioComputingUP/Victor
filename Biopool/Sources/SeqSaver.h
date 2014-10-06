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
*    Loads components (Atoms, Groups, etc.) in SEQ format.
*    SEQ format lists the aminoacids, one per line, followed by the 
*    torsion angle settings. 
*    Note: saveGroup() is not implemented, as it has no valid use.
*/

#ifndef _SEQ_SAVER_H_
#define _SEQ_SAVER_H_

// Includes:
#include <SideChain.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <Saver.h>
#include <Debug.h>
#include <string>
#include <iostream>
#include <Ligand.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/**@brief Loads components (Atoms, Groups, etc.) in SEQ format.
 * 
*@Description SEQ format lists the aminoacids, one per line, followed by the 
*    torsion angle settings. 
*@This 
 * */
class SeqSaver : public Saver{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  SeqSaver(ostream& _output = cout) 
    : output(_output), writeChi(true) { }
  // this class uses the implicit copy operator.
  virtual ~SeqSaver() { PRINT_NAME; }  

// MODIFIERS:
  void setWriteChi(bool _w) { writeChi = _w; }
  virtual void saveSideChain(SideChain& sc, bool header = 1);
  virtual void saveAminoAcid(AminoAcid& aa);
  virtual void saveSpacer(Spacer& sp);
  virtual void saveLigand(Ligand& l);

protected:

private:
  ostream& output;   // output stream
  bool writeChi;     // switch: write chi angles?
};

} // namespace
#endif //_SEQ_SAVER_H_
