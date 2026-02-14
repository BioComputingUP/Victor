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

#ifndef CIFLOADER_H
#define	CIFLOADER_H

// Includes:
#include <string.h>
#include <utility>
#include <Loader.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <LigandSet.h>
#include <Protein.h>

#include "CifStructure.h"


// Global constants, typedefs, etc. (to avoid):

namespace Victor {
    namespace Biopool {
        
        /**
         * @brief Loads components (Atoms, Groups, Spacer, etc.) in standard CIF format.
         */
        class CifLoader : public Loader {
        public:

            // CONSTRUCTORS/DESTRUCTOR:
            /**
             * Constructor.
             * @param _input = CIF file object
	     * @param _output = log file
             * @param _permissive = if true, allows loading residues with missing atoms
             * @param _noHAtoms = if true, doesn't load Hydrogens
             * @param _noHetAtoms = if true, doesn't load het atoms
             * @param _noSecondary = if true, doesn't load secondary structure (neither the one calculated from torsional angles nor the DSSP)
             * @param _noConnection = if true, doesn't connect residues
             * @param _noWater = if true, doesn't load water atoms
             * @param _verb = if true, verbose mode
             * @param _allChains = if true, loads all chains
             * @param _NULL = the name of the chain to be loaded, if not provided only loads the first chain
             * @param _onlyMetal = if true, load only metals as ligands
             * @param _noNucleotideChains = if true, doesn't load DNA/RNA chains
             */
            CifLoader(istream& _input = cin, ostream& output = cout, bool _permissive = false,
                    bool _noHAtoms = false, bool _noHetAtoms = false,
                    bool _noSecondary = false, bool _noConnection = false,
                    bool _noWater = true, bool _verb = true, bool _allChains = false,
                    string _NULL = "", bool _onlyMetal = false,
                    bool _noNucleotideChains = true);

            // this class uses the implicit copy operator.

            virtual ~CifLoader();

            // PREDICATES:

            bool isValid();
            void checkModel(); //to check input values ​​
            void checkAndSetChain(); //chosen by the user
            unsigned int getMaxModels();
            vector<char> getAllChains();

            // MODIFIERS:

            void setPermissive();
            void setNonPermissive();
            void setVerbose();
            void setNoVerbose();
            void setChain(char _ch);
            void setModel(unsigned int _mod);
            void setAltAtom(char _a);
            void setNoHAtoms();
            void setNoHetAtoms();
            void setOnlyMetalHetAtoms();
            void setNoSecondary();
            void setWithSecondary();
            void setNoConnection();
            void setWithConnection();
            void setWater();
            void setAllChains();

            //virtual void loadSpacer(Spacer& sp);
            //virtual void loadLigandSet(LigandSet& l);

            virtual void loadProtein(Protein& prot);

            //virtual void loadNucleotideChainSet(NucleotideChainSet& ns); //new class, new code by Damiano

        protected:
            // HELPERS:
            bool setBonds(Spacer& sp);
            bool inSideChain(const AminoAcid& aa, const Atom& at);
            void assignSecondary(Spacer& sp);
            int parseCifline(string atomLine, string tag, Ligand* lig, AminoAcid* aa);

            // ATTRIBUTES 
        private:
            istream& input; // input stream
	    ostream& output;
            bool permissive; //
            bool valid; //
            bool noHAtoms; //
            bool noHetAtoms; // hetatms contain water, simpleMetalIons and cofactors
            bool onlyMetalHetAtoms; // with this flag we select only 2nd cathegory
            bool noSecondary;
            bool noConnection; // skip connecting aminoacids
            bool noWater; // 
            bool verbose;
            bool allChains; //
            char chain; // chain ID to be loaded
            unsigned int model; // model number to be loaded
            char altAtom; // ID of alternate atoms to be loaded
            bool noNucleotideChains; // does not load nucleotide atoms

            string helixCode; // parallel vector of helix data, chain name for each helixData element
            string sheetCode;

            vector<pair<int, int> > helixData; // begin and end of the helix
            vector<pair<int, int> > sheetData;

            CifStructure* cif;
        };

        inline bool CifLoader::isValid() {
            return valid;
        }

        inline void CifLoader::setPermissive() {
            permissive = true;
        }

        inline void CifLoader::setNonPermissive() {
            permissive = false;
        }

        inline void CifLoader::setVerbose() {
            verbose = true;
        }

        inline void CifLoader::setNoVerbose() {
            verbose = false;
        }

        inline void CifLoader::setChain(char _ch) {
            chain = _ch;
        }

        inline void CifLoader::setModel(unsigned int _mod) {
            model = _mod;
        }

        inline void CifLoader::setAltAtom(char _a) {
            altAtom = _a;
        }

        inline void CifLoader::setNoHAtoms() {
            noHAtoms = true;
        }

        inline void CifLoader::setNoHetAtoms() {
            noHetAtoms = true;
        }

        inline void CifLoader::setNoSecondary() {
            noSecondary = true;
        }

        inline void CifLoader::setWithSecondary() {
            noSecondary = false;
        }

        inline void CifLoader::setNoConnection() {
            noConnection = true;
        }

        inline void CifLoader::setWithConnection() {
            noConnection = false;
        }

        inline void CifLoader::setAllChains() {
            allChains = true;
        }

    } // namespace Biopool
} // namespace Victor

#endif	/* CIFLOADER_H */

