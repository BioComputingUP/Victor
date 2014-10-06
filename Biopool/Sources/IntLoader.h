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

#ifndef _INT_LOADER_H_
#define _INT_LOADER_H_

// Includes:
#include <Spacer.h>
#include <LigandSet.h>
#include <Loader.h>
#include <Debug.h>
#include <string>
#include <iostream>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/**@brief Loads components (Atoms, Groups, etc.) in internal format.
 * 
*@Description Internal format is defined by listing type, bond length partner &
 *    bond length, bond angle partner & bond angle, torsion angle partner
 *    & torsion angle plus a chirality (0 if it is a 'true' torsion angle,
 *    +1 or -1 if the 'torsion angle' is a second bond angle), for each
 *    atom, one per line.
 *    NB: Only chirality 0 is currently supported.
*@This 
 * */
class IntLoader : public Loader
{
public: 

// CONSTRUCTORS/DESTRUCTOR:
  IntLoader(istream& _input = cin) : input(_input), atomIndex() { }
  // this class uses the implicit copy operator.
  virtual ~IntLoader() { PRINT_NAME; }  

// MODIFIERS:
  virtual void loadGroup(Group& group);
  virtual void loadSideChain(SideChain& sc,AminoAcid* aaRef);
  virtual void loadAminoAcid(AminoAcid& aa);
  virtual void loadSpacer(Spacer& sp);
  virtual void loadLigand(Ligand& l);

protected:
  void setBonds(Spacer& sp);
 
private:
  void zAtomToCartesian(Atom& atbLP, double bondLength, Atom& atbAP, 
			double bondAngle, Atom& attAP, double torsionAngle, 
			int chiral,Atom& at);
  bool inSideChain(const AminoAcid& aa, const Atom& at);
  Atom& addToAminoAcid(AminoAcid* aa,Atom* at);
  void keepPlacingAlongZAxis(int bLP, double bondLength, int bAP,
				 double bondAngle, Atom* at);
  istream& input;   // input stream
  vector<Atom> atomIndex;
  bool connect; // are segments to be connected to each other?
};

} // namespace
#endif //_INT_LOADER_H_











