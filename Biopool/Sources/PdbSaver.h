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

#ifndef _PDB_SAVER_H_
#define _PDB_SAVER_H_

// Includes:
#include <Group.h>
#include <SideChain.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <LigandSet.h>
#include <Ligand.h>
#include <Saver.h>
#include <Debug.h>
#include <string>
#include <iostream>
#include <Protein.h>
// Global constants, typedefs, etc. (to avoid):

namespace Victor { namespace Biopool { 

    /**@brief Saves components (Atoms, Groups, etc.) in standard PDB format
     * 
     * */
    class PdbSaver : public Saver {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        /**
         *   Basic constructor. By default it writes sequence, secondary structure and the term line. 
         * @param _output (ostream&) the output file object
         */
        PdbSaver(ostream& _output = cout)
        : output(_output), writeSeq(true), writeSecStr(true), writeTer(true),
        atomOffset(0), aminoOffset(0), ligandOffset(0), chain(' ') {
        }
        // this class uses the implicit copy operator.

        virtual ~PdbSaver() {
            PRINT_NAME;
        }

        // PREDICATES:

        void endFile() {
            output << "END\n";
        }

        // MODIFIERS:

        void setWriteSecondaryStructure() {
            writeSecStr = true;
        }

        void setDoNotWriteSecondaryStructure() {
            writeSecStr = false;
        }

        void setWriteSeqRes() {
            writeSeq = true;
        }

        void setDoNotWriteSeqRes() {
            writeSeq = false;
        }

        void setWriteAtomOnly() {
            writeSecStr = false;
            writeSeq = false;
            writeTer = false;
        }

        void setWriteAll() {
            writeSecStr = true;
            writeSeq = true;
            writeTer = true;
        }

        void setChain(char _ch) {
            chain = _ch;
        }
        virtual void saveGroup(Group& gr);
        virtual void saveSideChain(SideChain& sc);
        virtual void saveAminoAcid(AminoAcid& aa);
        virtual void saveSpacer(Spacer& sp);
        virtual void saveLigand(Ligand& l);
        virtual void saveLigandSet(LigandSet& l);
        virtual void saveProtein(Protein& prot);

    protected:

    private:

        // HELPERS:
        void writeSeqRes(Spacer& sp); // writes SEQRES entry
        void writeSecondary(Spacer& sp);
        // writes secondary entries (SHEET, HELIX, etc.)
        // ATTRIBUTES 
        ostream& output; // output stream
        bool writeSeq, writeSecStr, writeTer;
        unsigned int atomOffset, ligandOffset;
        int aminoOffset;
        char chain; // chain ID
        // offsets that determine at which atom, aminoacid and ligand number to start
    };

}} //namespace
#endif //_PDB_SAVER_H_


