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

#ifndef _SAVER_H_
#define _SAVER_H_

// Includes:
#include <Debug.h>

// Global constants, typedefs, etc. (to avoid):

namespace Victor { namespace Biopool { 

    class Group;
    class SideChain;
    class AminoAcid;
    class Spacer;
    class Ligand;
    class LigandSet;
    class Protein;
    class Nucleotide;

    /**@brief Base class for saving components (Atoms, Groups, etc.).
     * 
     * */
    class Saver {
    public:

        // CONSTRUCTORS/DESTRUCTOR:

        Saver() {
        };
        // this class uses the implicit copy operator.

        virtual ~Saver() {
        };

        // MODIFIERS:

        virtual void saveGroup(Group& group) {
        };

        virtual void saveSideChain(SideChain& node) {
        };

        virtual void saveAminoAcid(AminoAcid& node) {
        };

        virtual void saveSpacer(Spacer& node) {
        };

        virtual void saveLigand(Ligand& node) {
        };

        virtual void saveLigandSet(LigandSet& node) {
        };

        virtual void saveProtein(Protein& node) {
        };

        virtual void saveNucleotide(Nucleotide& node) {
        };
    protected:

    private:

    };

}} //namespace
#endif //_SAVER_H_










