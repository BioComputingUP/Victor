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
 * File:   CIFFileManager.h
 * Author: Francesco Menniti
 *
 * Created on 1 dicembre 2015, 9.52
 */

#ifndef CIFFILEMANAGER_H
#define CIFFILEMANAGER_H

#include <iostream>
#include <CIFLoop.h>
#include <StringUtils.h>
#include <IoTools.h>
#include <CIFDescriptor.h>
#include <map>

namespace Victor {
    namespace Biopool {
        
        /**@brief This class read the information of the CIF file. After the init process all the information are stored in the CIFFileManager object.
         */
        class CIFFileManager {
        public:
            
            /**
             * Constructor with the stream of the CIF file to parse
             * @param _input The CIF file to parse
             */
            CIFFileManager(istream& _input):input(_input){}
            
            /**
             * Standard copier
             */
            CIFFileManager(const CIFFileManager& orig);

            /**
             * This method return the value of one of non loop property in the CIF file
             * 
             * @param key The name of the property to get
             * @return the value of the input key
             */
            string GetInLineValue(string key);
            
            /**
             * This method start the parsing of the CIF file. After this operation the file information are stored in the CIFFileManager object.
             * Remark: call this method before any other operation.
             */
            void Init();
            
            /**
             * Store all the loops different from Atom, Helix and Sheet loop.
             */
            vector<CIFLoop> OtherLoops;
            
            /**
             * Store the information about the atom sequences
             */
            CIFLoop AtomLoop;
            
            /**
             * Store details about the backbone conformation of a segment of polymer.
             */
            CIFLoop HelixLoop;
            
            /**
             * Store details about the residue ranges that form a beta sheet.
             */
            CIFLoop SheetLoop;
            
        private:
            /**
             * Standard constructor
             */
            CIFFileManager();
            
            /**
             * Parse all the information abaout a single loop.
             * 
             * @return The loop parsed.
             */
            CIFLoop processLoop();
            
            /**
             * Parse a single line of the CIF file. This method process the non loop information of the CIF file.
             * Add the parsed information in the inLineValues HashMap.
             * 
             * @param line the line to parse
             */
            void processLine(string line);
            
            /**
             * Store the inLine information of the CIF file.
             */
            map<string,string> inLineValues;
            
            /**
             * The CIF file to parse
             */
            istream& input;
            
        };
    }
}
#endif /* CIFFILEMANAGER_H */

