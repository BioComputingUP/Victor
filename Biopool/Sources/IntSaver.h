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


#ifndef _INT_SAVER_H_
#define _INT_SAVER_H_

// Includes:
#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <AminoAcid.h>
#include <Saver.h>
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
     * 
     **/
    class IntSaver : public Saver {
    public:

        // CONSTRUCTORS/DESTRUCTOR:

        IntSaver(ostream& _output = cout) : output(_output) {
        };
        // this class uses the implicit copy operator.

        virtual ~IntSaver() {
            PRINT_NAME;
        };

        // MODIFIERS:
        virtual void saveGroup(Group& group);
        virtual void saveSideChain(SideChain& node);
        virtual void saveAminoAcid(AminoAcid& node);
        virtual void saveSpacer(Spacer& sp);
        virtual void saveLigand(Ligand& l);

    protected:
        void pSaveAtomVector(vector<Atom>& va, unsigned int offset = 0);

    private:
        ostream& output; // output stream
        // unsigned int count;
    };

} // namespace
#endif //_INT_SAVER_H_
