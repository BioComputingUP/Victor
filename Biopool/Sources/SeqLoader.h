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


#ifndef _SEQ_LOADER_H_
#define _SEQ_LOADER_H_

// Includes:
#include <SideChain.h>
#include <AminoAcid.h>
#include <Loader.h>
#include <Debug.h>
#include <string>
#include <vector>

// Global constants, typedefs, etc. (to avoid):

namespace Victor { namespace Biopool { 

    /**
     * @brief Loads components (Atoms, Groups, etc.) in SEQ format. 
     * 
     *   SEQ format lists the aminoacids, one per line, followed by the 
     *    torsion angle settings. In order to construct the protein, a 
     *    reference file with sample aminoacids has to be loaded first. 
     * */
    class SeqLoader : public Loader {
    public:

        // CONSTRUCTORS/DESTRUCTOR:

        SeqLoader(istream& _input = cin, istream& _refInput = cin) :
        input(_input), refInput(_refInput), loaded(false), refAmino() {
        }
        // this class uses the implicit copy operator.

        ~SeqLoader() {
            PRINT_NAME;
        }

        // MODIFIERS:
        virtual void loadAminoAcid(AminoAcid& node, AminoAcid* prev = NULL);
        virtual void loadSpacer(Spacer& node);
        virtual void loadLigand(Ligand& node);

    protected:
        // HELPERS:
        virtual void loadReference();
        virtual void setStructure(AminoAcid& aa, string type);

        // ATTRIBUTES:
        istream& input; // input stream
        istream& refInput; // reference stream
        int loaded; // is reference loaded yet?
        vector<AminoAcid> refAmino; // reference amino acids

    private:

    };

}} //namespace
#endif //_SEQ_LOADER_H_
