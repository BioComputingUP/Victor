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
 * File:   CIFLoop.h
 * Author: Francesco Menniti
 *
 * Created on 01 december 2015, 13.34
 */

#ifndef CIFLOOP_H
#define CIFLOOP_H

#include <iostream>
#include <string>
#include <vector>
#include <StringUtils.h>
#include <IoTools.h>

using namespace std;

namespace Victor {
    namespace Biopool {

        /**
         * @brief This class store and manage the CIF file lopps.
         * Loops are like table and stores values. Loops have 2 section: the header and the rows.
         * The column are self-defined in the header of the loop.
         * The rows contains the value of the table.
         */
        class CIFLoop {
        public:

            /**
             * The standard constructors.
             */
            CIFLoop() {};
            
            /**
             * This costructor make a categorized loop.
             * @param _category The category of the loop read from CIF file
             */
            CIFLoop(string _category);

            /**
             * Extracts the category of a loop from the first line of the loop.
             * @param line The first line of the loop.
             * @return The category parsed.
             */
            static string DecodeCategory(string line);

            /**
             * The category of the loop.
             * @return The category of the loop.
             */
            string GetCategoryGroup() {
                return categoryGroup;
            }

            /**
             * Parse a header line of the loop and extract the name of the field-column
             * @param _line The line to precess.
             */
            void AddColumn(string _line);
            
            /**
             * This method return a value of the loop row.
             * @param _colName The value of the column.
             * @param _rowNumber The index of the loop row 0 based.
             * @return The value as string
             */
            string GetCellValue(string _colName, int _rowNumber);
            
            /**
             * Like GetCellValue. Remove the white space, if exists, before and after the cell value.
             * @param _colName The value of the column.
             * @param _rowNumber The index of the loop row 0 based.
             * @return The value as string
             */
            string GetCellValueTrim(string _colName, int _rowNumber);
            
            /**
             * Check if a column name belong to the loop.
             * @param _colName The column name to check.
             * @return True if the loop has the column name.
             */
            bool HasColumn(string _colName);

            /**
             * Add a record to the loop.
             * @param _line The line of the CIF file to store in the loop.
             */
            void AddRow(string _line) {
                rows.push_back(splitRow(_line));
            }

            /**
             * Get how many records the loop store.
             * @return The record nmber of the loop.
             */
            int RowsCount() {
                return rows.size();
            }

            /**
             * Get how many columns the loop has.
             * @return The columns nmber of the loop.
             */
            int ColumnsCount() {
                return columns.size();
            }

        private:
            /**
             * The columns list of the loop.
             */
            vector<string> columns;
            
            /**
             * The record list of the loop.
             * Each record is rappresented as a list of values.
             */
            vector<vector<string> > rows;
            
            /**
             * The category of the loop.
             */
            string categoryGroup;

            /**
             * This methoed get a line of the CIF file and split the values on each white space.
             * @param _line The line to process
             * @return The values list splitted.
             */
            vector<string> splitRow(string _line);
            
            /**
             * Get a header line of a loop and return the namen of the column defined.
             * @param _line The line of CIF file to parse
             * @return The column name store in the header loop line in input.
             */
            string decodeColName(string _line);
            
            /**
             * Get the index of the input column name.
             * @param _colName The name to find.
             * @return The index of the colunm. -1 if not present.
             */
            int getColNumber(string _colName);

        }; //class CIFLoop
    }//namespace Biopool
}//namespace Victor

#endif // CIFLOOP_H
