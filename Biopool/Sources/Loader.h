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

#ifndef _LOADER_H_
#define _LOADER_H_

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

    /**@brief Base class for loading components (Atoms, Groups, etc.).
     * 
     * */
    class Loader {
    public:

        // CONSTRUCTORS/DESTRUCTOR:

        Loader() {
            PRINT_NAME;
        }
        // this class uses the implicit copy operator.

        virtual ~Loader() {
            PRINT_NAME;
        };

        // MODIFIERS:

        virtual void loadGroup(Group& group) {
            PRINT_NAME;
        };

        virtual void loadSideChain(SideChain& node) {
            PRINT_NAME;
        };

        virtual void loadAminoAcid(AminoAcid& node) {
            PRINT_NAME;
        };

        virtual void loadSpacer(Spacer& node) {
            PRINT_NAME;
        };

        virtual void loadLigand(Ligand& node) {
            PRINT_NAME;
        };

        virtual void loadLigandSet(LigandSet& node) {
            PRINT_NAME;
        };

        virtual void loadProtein(Protein& node) {
            PRINT_NAME;
        };

        virtual void loadNucleotide(Nucleotide& node) {
            PRINT_NAME;
        };
    protected:

    private:

    };


}} //namespace
#endif //_LOADER_H_









