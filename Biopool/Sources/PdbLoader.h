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

#ifndef _PDB_LOADER_H_
#define _PDB_LOADER_H_

// Includes:
#include <string.h>
#include <utility>
#include <Loader.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <LigandSet.h>
#include <Protein.h>

// Global constants, typedefs, etc. (to avoid):

namespace Victor { namespace Biopool { 

    /**@brief Loads components (Atoms, Groups, Spacer, etc.) in standard PDB format.
     *   
     * */
    class PdbLoader : public Loader {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        /**
         *   Constructor.
         * @param _input = the PDB file object
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
        PdbLoader(istream& _input = cin, bool _permissive = false,
                bool _noHAtoms = false, bool _noHetAtoms = false, bool _noSecondary = false,
                bool _noConnection = false, bool _noWater = true,
                bool _verb = true, bool _allChains = false, string _NULL = "", bool _onlyMetal = false,
                bool _noNucleotideChains = true)
        : input(_input), permissive(_permissive), valid(true),
        noHAtoms(_noHAtoms), noHetAtoms(_noHetAtoms), noSecondary(_noSecondary),
        noConnection(_noConnection), noWater(_noWater), verbose(_verb),
        allChains(_allChains), chain(' '), model(999), altAtom('A'), helixCode(_NULL),
        //sheetCode(_NULL), helixData(), sheetData(), onlyMetalHetAtoms(_onlyMetal), 
        sheetCode(_NULL), onlyMetalHetAtoms(_onlyMetal),
        noNucleotideChains(_noNucleotideChains) {
        }

        // this class uses the implicit copy operator.

        virtual ~PdbLoader() {
            PRINT_NAME;
        }

        // PREDICATES:

        bool isValid() {
            return valid;
        }
        void checkModel(); //to check input values ​​
        void checkAndSetChain(); //chosen by the user
        unsigned int getMaxModels();
        unsigned int getMaxModelsFast();
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



        //virtual void loadSpacer(Spacer& sp);
        //virtual void loadLigandSet(LigandSet& l);
        virtual void loadProtein(Protein& prot);



        //virtual void loadNucleotideChainSet(NucleotideChainSet& ns); //new class, new code by Damiano

    protected:
        // HELPERS:
        bool setBonds(Spacer& sp);
        bool inSideChain(const AminoAcid& aa, const Atom& at);
        void loadSecondary();
        void assignSecondary(Spacer& sp);
        int parsePDBline(string atomLine, string tag, Ligand* lig, AminoAcid* aa);



        // ATTRIBUTES 
    private:
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

    };

}} //namespace
#endif //_PDB_LOADER_H_

