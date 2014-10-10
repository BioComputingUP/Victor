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


#ifndef _XYZ_LOADER_H_
#define _XYZ_LOADER_H_

// Includes:
#include <Group.h>
#include <Ligand.h>
#include <SideChain.h>
#include <Loader.h>
#include <AminoAcid.h>
#include <Spacer.h>

// Global constants, typedefs, etc. (to avoid):

namespace Victor { namespace Biopool { 

    /**
     * @brief Loads components (Atoms, Groups, etc.) in XYZ (carthesian) format.
     * 
     *  
     *    Internal format is made of type, coords & bonds of each atom,
     *    one per line. Keywords "aminoacid" and "sidechain" delimit these
     *    structures.
     * */
    class XyzLoader : public Loader {
    public:

        // CONSTRUCTORS/DESTRUCTOR:

        XyzLoader(istream& _input = cin) : input(_input), connect(true) {
        }
        // this class uses the implicit copy operator.

        virtual ~XyzLoader() {
            PRINT_NAME;
        }

        // MODIFIERS:

        void connectSegments(bool c) {
            connect = c;
        }
        // determines if segments (aminoacids, sidechains) have to be connected
        virtual void loadGroup(Group& gr);
        virtual void loadSideChain(SideChain& sc, AminoAcid* aaRef = NULL);
        virtual void loadAminoAcid(AminoAcid& aa);
        virtual void loadSpacer(Spacer& sp);
        virtual void loadLigand(Ligand& l);

    protected:

    private:
        istream& input; // input stream
        bool connect; // are segments to be connected to each other?
    };
}} //namespace
#endif //_XYZ_LOADER_H_
