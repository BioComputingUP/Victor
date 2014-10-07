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
/**
 *@Class               EnergyAnalyzer
 *@Project    Victor
 */

#ifndef _ENERGYANALYZER_H_
#define _ENERGYANALYZER_H_

// Includes:
#include <string>
#include <vector>

using namespace std;

namespace Biopool {
    /**@brief This class implements functions to analyze the ranking output of Energy. 
     * 
     *@Description  
     * */
    // Global constants, typedefs, etc. (to avoid):

    class EnergyAnalyzer {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        EnergyAnalyzer();
        EnergyAnalyzer(const EnergyAnalyzer& orig);
        virtual ~EnergyAnalyzer();

        // PREDICATES:
        void printHeader();
        void printAllCorrelations();
        void printCorrelation(char* intro, unsigned int index);
        void printCorrelation(char* intro, vector<double> data);
        void printTopXResults(unsigned int top, string topFile,
                double maxScore = 1000.0);

        // MODIFIERS:
        void load(istream& inFile, unsigned int select, bool mult = false);
        void calcStatistics();
        vector<double> linearOptimization(double a1, double a2, double a3, double b1,
                double b2, double b3, double c1, double c2,
                double c3, double d1, double d2, double d3);
        vector<double> nonLinearOptimization(unsigned int max1, unsigned int max2);

        void copy(const EnergyAnalyzer& orig);

        // OPERATORS:

        //  protected:

        //  private:

        // HELPERS: 

        // ATTRIBUTES:
        static const unsigned int MAX_COL = 6;
        static const unsigned int MAX_DATA = 5000;

        double avg[MAX_COL];
        double sd[MAX_COL];
        vector<vector<double> > zScore; // Z-scores of the data ("col")
        vector<vector<double> > col; // raw input data (i.e. "column x")
        vector<vector<double> > result; // GDT_TS result, used to get Top X 
        vector<vector<double> > resScore; // score of the GDT_TS results
    };

} // namespace
#endif //_ENERGYANALYZER_H_


