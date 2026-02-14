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

#ifndef CIFSTRUCTURE_H
#define	CIFSTRUCTURE_H


// Includes:
#include <string>
#include <vector>
#include <iostream>

using std::string;
using std::istream;
using std::ostream;
using std::vector;

// Global constants, typedefs, etc. (to avoid):

namespace Victor {
    namespace Biopool {

        /**
         * @brief Helper class used to hold information from CIF file.
         */
        class CifStructure {
        public:
            /**
             * Constructor.
             * @param input input file stream
	     * @param output output file stream
             */
            CifStructure(istream& input, ostream& output = cout);
	    
	    /**
	     * Constructor.
             * @param output output file stream
	     * @param input input file stream
             */
	    CifStructure(ostream& output);
            
            /**
             * Destructor.
             */
            virtual ~CifStructure();
            
            /**
            * Returns the correct collection by group name.
            * @param name name of the CIF group
            * @return reference to the collection
            */
            vector<string>& getGroup(string name);
            
            /**
            * Returns the tag by name.
            * @param name name of tag
            * @return CIF tag
            */
            string getTag(string name);
            
            /**
            * Returns the column number of the field.
            * @param name name of the group
            * @param field name of the field
            * @return field column number
            */
            int getGroupColumnNumber(string name, string field);
            
            /**
            * Returns the field of the line at the columnNum column.
            * @param name name of the group
            * @param line line of the CIF
            * @param columnNum number of column
            * @return field at columnNum column
            */
            string getGroupField(string name, string& line, int columnNum);
            
	    /**
	     * Returns the field present in the CIF line.
             * @param line
             * @return 
             */
	    string getInlineField(string& line);
	    
            /**
            * Parses group of CIF fields and creates a vector with columns positions.
            * @param name name of the group 
            */
            void parseGroup(string group, string& line);
            
            /**
            * Sets the flag of the parsed group.
            * @param name name of the group
            */
            void setParsedFlag(string name);
            
            /**
             * Returns true if the group name is parsed, false otherwise.
             * @param name name of the group
             * @return true if group is parsed, false otherwise
             */
            bool isGroupParsed(string name);
	    
	    /**
	     * Prints group records names into output stream.
             * @param name name of the group
             */
	    void printGroup(string name);
	    
	    /**
	     * Returns the input stream.
             * @return input stream
             */
	    istream& getInput();
            
        private:
	    /**
	     * Initializes data members.
             */
	    void initData();
	    
            // CIF file
            istream& input;
	    ostream& output;

            // CIF tags
            string header;
	    
	    // atom group
            string atom;
	    string atomId;
	    string chain;
	    string asymId;
	    string entityId;
	    string residueIns;
	    string x;
            string y;
            string z;
	    string occupancy;
	    string tempFactor;
            string residueNum;
	    string residueName;
            string atomName;
            string model;
	    
	    // helix group
	    string helix;
            string helixStart;
            string helixEnd;
            string helixChainId;
	    
	    // sheet group
            string sheet;
            string sheetOrder;
            string sheetRange;
            string sheetHbond;
            string sheetStart;
            string sheetEnd;
            string sheetChainId;
            
            // collections of CIF group fields
            vector<string> atomGroup;
            vector<string> helixGroup;
            vector<string> sheetGroup;
            vector<string> sheetOrderGroup;
            vector<string> sheetRangeGroup;
            vector<string> sheetHbondGroup;

            // flags
            bool atomGroupParsed;
            bool helixGroupParsed;
            bool sheetGroupParsed;
            bool sheetOrderGroupParsed;
            bool sheetRangeGroupParsed;
            bool sheetHboundgroupParsed;
        };
	
	inline istream& CifStructure::getInput() {
	    return input;
	}
    } // namespace Biopool
} // namespace Victor 

#endif	/* CIFSTRUCTURE_H */

