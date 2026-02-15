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
 * File: CIFSaver.h
 * Author: Francesco Menniti
 *
 * Created on 22 dicembre 2015, 15.28
 */

#ifndef CIFSAVER_H
#define CIFSAVER_H

#include <Group.h>
#include <SideChain.h>
#include <AminoAcid.h>
#include <Spacer.h>
#include <LigandSet.h>
#include <Ligand.h>
#include <Saver.h>
#include <Debug.h>
#include <Protein.h>
#include <StringUtils.h>

namespace Victor {
    namespace Biopool {

        /**@brief Saves components (Atoms, Groups, etc.) in standard CIF format
         */
        class CIFSaver : public Saver {
        public:

            // CONSTRUCTORS/DESTRUCTOR:

            /**
             *Basic constructor. By default it writes sequence, secondary structure and the term line. 
             * @param _output (ostream&) the output file object
             */
            CIFSaver(ostream& _output = cout)
            : output(_output), writeSeq(true), writeSecStr(true), writeTer(true),
            atomOffset(0), aminoOffset(0), ligandOffset(0), chain(' '), atomSiteWrited(false) {
            }

            virtual ~CIFSaver() {
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

            /**
             *Saves a group in CIF format.
             *@param group reference 
             *@return void
             */
            virtual void saveGroup(Group& gr);

            /**
             *Saves a sidechain in CIF format. 
             *@param sideChain reference 
             *@return void
             */
            virtual void saveSideChain(SideChain& sc);

            /**
             *Saves an aminoacid in CIF format.
             *@param AminoAcid reference 
             *@return void
             */
            virtual void saveAminoAcid(AminoAcid& aa);

            /**
             *Saves a spacer in CIF format. 
             *@param Spacer reference 
             *@return void
             */
            virtual void saveSpacer(Spacer& sp);

            /**
             *Saves a Ligand in CIF format. 
             *@param Ligand reference 
             *@return void
             */
            virtual void saveLigand(Ligand& l);

            /**
             *Saves a LigandSet in CIF format. 
             *@param LigandSet reference 
             *@return void
             */
            virtual void saveLigandSet(LigandSet& l);

            /**
             *Saves a Protein in CIF format. 
             *@param Protein reference 
             *@return void
             */
            virtual void saveProtein(Protein& prot);

        protected:

        private:

            // HELPERS:
            /**
             *Writes the SEQRES entry (CIF format) for a spacer.
             *@param Spacer reference 
             *@return void
             */
            void writeSeqRes(Spacer& sp); // writes SEQRES entry

            /**
             *Writes the secondary information (CIF format) for a spacer, e.g. HELIX, SHEET, etc.
             *@param sideChain reference 
             *@return void
             */
            void writeSecondary(Spacer& sp);
            // ATTRIBUTES 
            ostream& output; // output stream
            bool writeSeq, writeSecStr, writeTer, atomSiteWrited;
            unsigned int atomOffset, ligandOffset;
            int aminoOffset;
            char chain; // chain ID
            // offsets that determine at which atom, aminoacid and ligand number to start
        };

    }
}
#endif /* CIFSAVER_H */

