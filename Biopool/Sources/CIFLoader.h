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

/*
 * File:   CifLoader.h
 * Author: Francesco Menniti
 *
 * Created on 30 novembre 2015, 15.40
 */

#ifndef CIFLOADER_H
#define CIFLOADER_H

// Includes:
#include <AminoAcid.h>
#include <AminoAcidHydrogen.h>
#include <AtomCode.h>
#include <CIFDescriptor.h>
#include <CIFFileManager.h>
#include <IoTools.h>
#include <Ligand.h>
#include <LigandSet.h>
#include <Loader.h>
#include <Nucleotide.h>
#include <Protein.h>
#include <Spacer.h>
#include <String2Number.h>
#include <string.h>
#include <utility>
#include <vector3.h>

namespace Victor {
    namespace Biopool {

        /**@brief Loads components (Atoms, Groups, Spacer, etc.) in standard CIF format.
         */
        class CIFLoader : public Loader {
        public:

            // CONSTRUCTORS/DESTRUCTOR:
            /**
             * Constructor.
             * @param _input = the CIF file object
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
            CIFLoader(istream& _input = cin, bool _permissive = false,
                    bool _noHAtoms = false, bool _noHetAtoms = false, bool _noSecondary = false,
                    bool _noConnection = false, bool _noWater = true,
                    bool _verb = true, bool _allChains = false, string _NULL = "", bool _onlyMetal = false,
                    bool _noNucleotideChains = true);

            // this class uses the implicit copy operator.

            virtual ~CIFLoader() {
                PRINT_NAME;
            }

            // PREDICATES:

            bool isValid() {
                return valid;
            }

            /**
             * If user selected a Model, it check validity of this choice,
             * otherwise it select first available chain. To check input values.
             */
            void checkModel();

            /**
             * If user selected a chain, it check validity of this choice,
             * otherwise it select first available chain. Chosen by the user.
             */
            void checkAndSetChain();

            /**
             * Reads in the maximum allowed number of NMR models, zero otherwise
             */
            unsigned int getMaxModels();

            /**
             * Returns all available chain IDs for a CIF file.
             *
             */
            vector<char> getAllChains();

            // MODIFIERS:

            void setPermissive() {
                permissive = true;
            }

            void setNonPermissive() {
                permissive = false;
            }

            void setVerbose() {
                verbose = true;
            }

            void setNoVerbose() {
                verbose = false;
            }

            void setChain(char _ch) {
                chain = _ch;
            }

            void setModel(unsigned int _mod) {
                model = _mod;
            }

            void setAltAtom(char _a) {
                altAtom = _a;
            }

            void setNoHAtoms() {
                noHAtoms = true;
            }

            void setNoHetAtoms() {
                noHetAtoms = true;
            }
            void setOnlyMetalHetAtoms();

            void setNoSecondary() {
                noSecondary = true;
            }

            void setWithSecondary() {
                noSecondary = false;
            }

            void setNoConnection() {
                noConnection = true;
            }

            void setWithConnection() {
                noConnection = false;
            }

            void setWater();

            void setAllChains() {
                allChains = true;
            }

            /**
             * Core function for CIF file parsing.
             * @param prot (Protein&)
             */
            virtual void loadProtein(Protein& prot);

        protected:
            // HELPERS:
            /**
             * Private helper function to set bond structure after loading the spacer.
             * @param   Spacer reference
             * @return  bool
             */
            bool setBonds(Spacer& sp);

            /**
             * Private helper function to determine if atom is backbone or sidechain.
             * @param   Spacer reference
             * @return  bool
             */
            bool inSideChain(const AminoAcid& aa, const Atom& at);
            void loadSecondary();

            /**
             * Try to assigns the secondary structure from the CIF header. If not present
             * uses Spacer's setStateFromTorsionAngles().
             * @param   Spacer reference
             */
            void assignSecondary(Spacer& sp);

            /**
             * Parse a single line of a CIF file.
             * @param atomLine (string) the whole CIF line as it is
             * @param tag (string) = the first field (keyword) in a CIF line
             * @param lig (Ligand) pointer
             * @param aa (AminoAcid) pointer
             * @return Residue number read from the CIF line (int)
             */
            int parseCIFline(int atomLineIndex, string tag, Ligand* lig, AminoAcid* aa);

            // ATTRIBUTES
        private:
            CIFFileManager* cifFM;
            istream& input; //input stream
            bool permissive; //
            bool valid; //
            bool noHAtoms; //
            bool noHetAtoms; //hetatms contain water, simpleMetalIons and cofactors
            bool onlyMetalHetAtoms; //with this flag we select only 2nd cathegory
            bool noSecondary;
            bool noConnection; //skip connecting aminoacids
            bool noWater; //
            bool verbose;
            bool allChains; //
            char chain; //chain ID to be loaded
            unsigned int model; //model number to be loaded
            char altAtom; //ID of alternate atoms to be loaded

            bool noNucleotideChains; //does not load nucleotide atoms

            string helixCode; // parallel vector of helix data, chain name for each helixData element
            string sheetCode;

            vector<pair<int, int> > helixData; //inizio e fine dell'elica
            vector<pair<int, int> > sheetData;

            /**
             * A loadProtein subroutine. Analyze and parse the helix information.
             */
            void parseHelixData();
            
            /**
             * A loadProtein subroutine. Analyze and parse the sheet information.
             */
            void parseSheetData();
            
            /**
             * A loadProtein subroutine. Analyze and parse the information for atoms and hetero atoms.
             * @param currentChain The current chain.
             * @param aaNum The current aminiacid serial number
             * @param oldAaNum The last aminiacid serial number
             * @param aa The aminoacid to create
             * @param sp The spacer to create
             * @param lig The ligand to create
             * @param ls The ligand set to create
             */
            void parseAtomHetatmData(const char currentChain, int &aaNum, int &oldAaNum, AminoAcid* aa, Spacer* sp, Ligand* lig, LigandSet* ls);

        };
    }
}

#endif /* CIFLOADER_H */

