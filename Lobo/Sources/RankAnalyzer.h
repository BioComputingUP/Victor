/* 
 *  This file is part of Victor.
 *
 *    Victor is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    Victor is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.

 *    You should have received a copy of the GNU General Public License
 *    along with Victor.  If not, see <http://www.gnu.org/licenses/>.
 */
/** 
 *@Class:              RankAnalyzer
 */

#ifndef _RANKANALYZER_H_
#define _RANKANALYZER_H_

// Includes:
#include <string>
#include <vector>
using namespace std;

namespace Biopool {

    // Global constants, typedefs, etc. (to avoid):

/**@brief  This class implements functions to analyze the ranking output of lobo.
     *  
     * */
    class RankAnalyzer {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        RankAnalyzer();
        RankAnalyzer(const RankAnalyzer& orig);
        virtual ~RankAnalyzer();

        // PREDICATES:
        void printHeader();
        void printCorrelation(char* intro, unsigned int index);
        void printCorrelation(char* intro, vector<double> data);
        void printTopXResults(unsigned int top, string topFile,
                double maxScore = 1000.0);

        // MODIFIERS:
        void load(istream& inFile, unsigned int select);
        void calcStatistics();

        void copy(const RankAnalyzer& orig);

        // OPERATORS:

        //  protected:

        //  private:

        // HELPERS: 

        // ATTRIBUTES:
        static const unsigned int MAX_COL = 11;
        static const unsigned int MAX_LOOP = 1000;

        double avg[MAX_COL];
        double sd[MAX_COL];
        vector< vector<double> > zScore; // Z-scores of the data ("col")
        vector< vector<double> > col; // raw input data (i.e. "column x")
        vector< vector<double> > result; // RMSD result, used to get Top X 
        vector< vector<double> > resScore; // score of the RMSD results
    };

} // namespace
#endif //_RANKANALYZER_H_


