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


#ifndef _globalStatistic_h_
#define _globalStatistic_h_

// Includes:
#include "LoopExtractor.h"
#include <set>
#include "ranking_helper2.h"

// Global constants, typedefs, etc. (to avoid):
namespace Victor { namespace Lobo {

    /**@brief  Methods to manages the global statistic data.
     * 
     *@Description  Used to work with the data for the global statistic in the 
     *   prop_calibration program.
     * */
    class globalStatistic {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        globalStatistic();
        ~globalStatistic();

        // PREDICATES:
        void outputArray(ofstream &outF);

        static const int LOOP_MIN = 3; // lower limit of examined Loops
        static const int LOOP_MAX = 13; // upper limit of examined loops

        // MODIFIERS:
        void updateVPArray(vector<double> &vdw_percent,
                vector<double> &prop_percent, int loopNr);
        void updateRMSArray(int loopNr, multiset<ranking_helper2> rmsh);
        void updateBadValuesArray(int loopNr,
                multiset<ranking_helper> prop_ranking,
                multiset<ranking_helper2> vdw_ranking,
                multiset<ranking_helper2> rmsh);
        void updateBadConsistencyArray(int loopNr, vector<int> consistency,
                multiset<ranking_helper2> rmsh);

        void updateEndRMSArray(vector<double> &endRMS_percent, int loopNr);
        void updateCompactnessArray(vector<double> &compactness_percent,
                int loopNr);
        void updateEnergyArray(vector<double> &energy_percent, int loopNr);

        void VdwCutoffGenerator(multiset<ranking_helper2> &rmsh,
                multiset<ranking_helper2> &rhs, int count,
                int loopNr);
        void PropCutoffGenerator(multiset<ranking_helper2> &rmsh,
                multiset<ranking_helper> &rhs, int count,
                int loopNr);
        void CompactnessCutoffGenerator(multiset<ranking_helper2> &rmsh,
                multiset<ranking_helper2> &rhs, int count,
                int loopNr);
        void EndRMSCutoffGenerator(multiset<ranking_helper2> &rmsh,
                multiset<ranking_helper2> &rhs, int count,
                int loopNr);
        void EnergyCutoffGenerator(multiset<ranking_helper2> &rmsh,
                multiset<ranking_helper2> &rhs, int count,
                int loopNr);
        void FilterCutoffGenerator(vector<double> values, int count, int loopNr);
        void updateRankedRmsArray(int loopNr, vector<double> rms);
        void updateSidechainArray(vector<int> collisionCount,
                multiset<ranking_helper2> rmsh, int loopNr);


        // OPERATORS:

    protected:

    private:

        // HELPERS:
        int calcRmsPercent(int index, multiset<ranking_helper2> rmsh);
        double calcDeviation(vector<double> *values);
        double calcMedian(vector<double> *values);
        double giveRms(int index, multiset<ranking_helper2> &rmsh);

        // ATTRIBUTES:

        static const int BAD_PROPENSITY; // value from which we have
        // a value to put into the bad prop statistic
        static const int BAD_VDW;
        // contains the value from which we have a value to put into the bad 
        // VDW statistic
        static const int BAD_CONSISTENCY;
        // contains the start value from which we put
        // a consistency value into the badConsistency array

        vector<double> * filterCutoff[LOOP_MAX - LOOP_MIN + 1][20];
        // used to calculate the filter-combination cutoff 

        double statArrayRMS[LOOP_MAX - LOOP_MIN + 1][6];
        // contains the statistic for the RMS @ 6!!
        int statArrayRMSCount[LOOP_MAX - LOOP_MIN + 1][6];
        // counts the entries in statArrayRMS @ 6!!

        double statArrayRankedRms[LOOP_MAX - LOOP_MIN + 1][6];
        // contains the rms values of the surviving solutions 
        // ranked by the filter results
        int statArrayRankedRmsCount[LOOP_MAX - LOOP_MIN + 1][6];
        // contains the number of entries in statArrayRankedRms

        long solutionSum; // contains the sum of all surviving solutions
        int loopCount; // contains the number of loops we have examined

        vector<double> *deviationRms[LOOP_MAX - LOOP_MIN + 1][6];
        // used to calculate the standard deviation for the RMS
        vector<double> *deviationRankedRms[LOOP_MAX - LOOP_MIN + 1][6];
        // used to calculate the standard deviation for the ranked RMS


        // counts the entries in statArrayEndRMS @ 6!!

    };


}}
#endif
