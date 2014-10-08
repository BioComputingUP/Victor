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
#ifndef _REL_SAVER_H_
#define _REL_SAVER_H_

// Includes:
#include <vector>
#include <Group.h>
#include <SideChain.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <Saver.h>
#include <Debug.h>
#include <string>
#include <iostream>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

    /**@brief Save components (Atoms, Groups, etc.) in relative format 
     * 
     *@Description The coordinates here are relative to the previous connected Component (Atom, Group, etc..:) 
     * rather than to the absolute Cartesian coordinate system.
     * 
     * */
    class RelSaver : public Saver {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        /**
         * @Description Basic constructor.
         * @param _output The output stream object
         * @param _offset Not implemented
         */
        RelSaver(ostream& _output = cout, int _offset = 1)
        : output(_output), offset(_offset) {
        }
        // this class uses the implicit copy operator.

        virtual ~RelSaver() {
            PRINT_NAME;
        }

        // MODIFIERS:
        virtual void saveGroup(Group& gr);
        virtual void saveSideChain(SideChain& sc);
        virtual void saveAminoAcid(AminoAcid& aa);
        virtual void saveSpacer(Spacer& sp);

    protected:
        virtual void pSaveAtomVector(vector<Atom>& va);

    private:
        ostream& output; // output stream
        int offset; // ID offset for saving (optional) -- currently disabled
    };

} // namespace
#endif //_REL_SAVER_H_
