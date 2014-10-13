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


// Includes:
#include <globalStatistic.h>


// Global constants, typedefs, etc. (to avoid):
using namespace Victor;
using namespace Victor::Lobo;
using namespace Victor::Biopool;
const int globalStatistic::BAD_PROPENSITY = 0;
const int globalStatistic::BAD_VDW = 3000;
const int globalStatistic::BAD_CONSISTENCY = 1000;


// CONSTRUCTORS/DESTRUCTOR:

/**
 *@Description Basic contructor
 */
globalStatistic::globalStatistic() {

    solutionSum = 0;
    loopCount = 0;
    LoopExtractor le;

    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
        for (int j = 0; j < le.BEST_COUNT; j += 1) {
            statArrayRMS[i][j] = 0;
            statArrayRMSCount[i][j] = 0;
            statArrayRankedRms[i][j] = 0;
            statArrayRankedRmsCount[i][j] = 0;
            deviationRms[i][j] = new vector<double>;
            deviationRankedRms[i][j] = new vector<double>;
        }
    }
    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
        for (int j = 0; j < 20; j += 1) {

            filterCutoff[i][j] = new vector<double>;
        }
    }

}

/**
 *@Description Basic destructor
 */
globalStatistic::~globalStatistic() {
    LoopExtractor le;
    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
        for (int j = 0; j < le.BEST_COUNT; j += 1) {
            delete deviationRms[i][j];
            delete deviationRankedRms[i][j];
        }
    }

}

// PREDICATES:

/**
 *@Description The arrays in this class are written to the file the filehandle refers to.
 *@param  reference to the output file(ofstream &)
 *@return   changes are made internally(void)
 */
void globalStatistic::outputArray(ofstream &outF) {

    LoopExtractor le;
    outF << "\nRMS" << "\n---\n\n";
    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
        outF << "Loop der Laenge: " << i + LOOP_MIN << "\n";
        for (int j = 0; j < le.BEST_COUNT; j += 1) {
            outF << "Durchschnittliche RMS fuer Loesung Nr. " << j + 1 << ": ";
            if (statArrayRMSCount[i][j] != 0)
                outF << (statArrayRMS[i][j] / static_cast<double> (
                    statArrayRMSCount[i][j])) << "\n";
            else
                outF << "kein Eintrag\n";
        }
    }

    outF << "\n" << "RMS der Loesungen nach dem Filterranking" << "\n";
    outF << "----------------------------------------" << "\n" << "\n";
    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
        outF << "Loop der Laenge: " << i + LOOP_MIN << "\n";
        for (int j = 0; j < 6; j += 1) {
            outF << "Durchschnittliche RMS fuer Loesung Nr. " << j + 1 << ": ";
            if (statArrayRankedRmsCount[i][j] != 0)
                outF << (statArrayRankedRms[i][j]
                    / static_cast<double> (statArrayRankedRmsCount[i][j])) << "\n";
            else
                outF << "kein Eintrag\n";
        }
    }

    outF << "\nBest RMSof the n best Filter ranked values\n"
            << "------------------------------------\n\n";
    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
        outF << "Loop der Laenge: " << i + LOOP_MIN << "\n";
        for (int j = 0; j < 20; j += 1) {
            outF << "Median for the " << j + 1 << " best values"
                    << calcMedian(filterCutoff[i][j]) << "\n";
        }
        outF << "\n";
    }

    outF << "\nNumber of examined loops"
            << "\n------------------------\n";
    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
        outF << "Loop der Laenge: " << i + LOOP_MIN << " Anzahl: "
                << statArrayRMSCount[i][0] << "\n";
    }

    outF << "\n" << "As an average we have " <<
            static_cast<double> (solutionSum) / static_cast<double> (loopCount) <<
            " solutions remaining" << "\n";

    outF << "\n" << "Standard Deviation for the RMS" << "\n";
    outF << "------------------------------" << "\n";
    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
        outF << "Loop der Laenge: " << i + LOOP_MIN << "\n";
        for (int j = 0; j < le.BEST_COUNT; j += 1) {
            outF << "Standard Deviation for solution nr. " << j + 1 << ": ";
            outF << calcDeviation(deviationRms[i][j]) << "\n";
        }
    }

    outF << "\n" << "Standard Deviation for the RankedRMS" << "\n";
    outF << "------------------------------------" << "\n";
    for (int i = 0; i < LOOP_MAX - LOOP_MIN + 1; i += 1) {
        outF << "Loop der Laenge: " << i + LOOP_MIN << "\n";
        for (int j = 0; j < le.BEST_COUNT; j += 1) {
            outF << "Standard Deviation for solution nr. " << j + 1 << ": ";
            outF << calcDeviation(filterCutoff[i][j]) << "\n";
            //      outF << calcDeviation(deviationRankedRms[i][j]) << "\n";
        }
    }

    outF.close();
}

/**
 *@Description This function iterates over the rmsh array until BEST_COUNT is reached.
 *  In each loop the values of the RMS of the best BEST_COUNT solutions
 *  are stored. Also the array which contains the number of entries is updated.
 *  If we don't have le.BEST_COUNT solutions left, we fill the rest of the 
 *  arrays with the worst rms value.
 *@param  loop index(int), a set of ranking helpers(multiset<ranking_helper2> )
 *@return   changes are made internally(void)
 */
void globalStatistic::updateRMSArray(int loopNr, multiset<ranking_helper2> rmsh) {
    int zaehler = 0;
    set<ranking_helper2>::iterator pos; // used to iterate over the rmsh
    LoopExtractor le; // needed in order to obtain BEST_COUNT
    double worst_value = 0; // used to fill the rest of the arrays

    for (pos = rmsh.begin(); pos != rmsh.end(); ++pos) {
        if (zaehler >= le.BEST_COUNT)
            break; // we calculate these values only for the best solutions

        (deviationRms[loopNr][zaehler])->push_back(pos->get_value());
        statArrayRMS[loopNr][zaehler] += pos->get_value();
        statArrayRMSCount[loopNr][zaehler] += 1;
        worst_value = pos->get_value();
        zaehler += 1;
    }

    //  check whether we have updated all le.BEST_COUNT solutions
    // if not,   adds the worst of all solutions to it
    for (int i = zaehler; i < le.BEST_COUNT; i += 1) {
        statArrayRMS[loopNr][i] += worst_value;
        (deviationRms[loopNr][i])->push_back(worst_value);
        statArrayRMSCount[loopNr][i] += 1;
    }
}

/**
 *@Description This function iterates over the rmsh array until BEST_COUNT is reached.
 *  In each loop the values of the RMS of the best BEST_COUNT solutions
 *  are stored. Also the array which contains the number of entries is updated
 *  If we have less than 6 surviving solutions, we take the worst solution
 *  and add it to the remaining arrays. Besides this, we update the statistic for the average number of surviving
 *  solutions.
 *@param  index of the loop (int), pointer to vector that contains the rms values(vector<double>)
 *@return  changes are made internally(void)
 */
void globalStatistic::updateRankedRmsArray(int loopNr, vector<double> rms) {
    if (rms.size() == 0)
        ERROR("We don't have a solution!", exception);

    unsigned int count = rms.size();
    double worstValue = 0;

    if (rms.size() >= 7)
        count = 6;

    // now find a starting value for the worst value in the rms vector
    worstValue = rms[0];

    for (unsigned int i = 0; i < count; i += 1) {
        statArrayRankedRms[loopNr][i] += rms[i];
        (deviationRankedRms[loopNr][i])->push_back(rms[i]);
        statArrayRankedRmsCount[loopNr][i] += 1;

        if (worstValue < rms[i])
            worstValue = rms[i];
    }

    // fill rest of array in case less than 6 solutions left
    for (unsigned int i = count; i < 6; i += 1) {
        statArrayRankedRms[loopNr][i] += worstValue;
        (deviationRankedRms[loopNr][i])->push_back(worstValue);
        statArrayRankedRmsCount[loopNr][i] += 1;
    }

    // update statistics for average number of surviving solutions
    loopCount += 1;
    solutionSum += rms.size();
}

/**
 *@Description We check the first count values in the vector values for the minimum.
//  This minimum is put into the proper array.
 *@param   pointer to vector that contains the values(vector<double>) , elements in the vector(int),index for the filter cut off(int)
 *@return  changes are made internally(void)
 */
void globalStatistic::FilterCutoffGenerator(vector<double> values, int count,
        int loopNr) {
    if ((count > 20) || (values.size() == 0))
        return;

    double min = values[0];

    for (int i = 1; i < count; i += 1)
        if (values[i] < min)
            min = values[i];

    (filterCutoff[loopNr][count - 1])->push_back(min);
}

/**
 *@Description Calculates the standard deviation of the values in the given Vector
 *@param   pointer to vector that contains the values(vector<double>)
 *@return  Deviation value (double)
 */
double globalStatistic::calcDeviation(vector<double> *values) {

    double median = 0;
    double deviation = 0;

    // first we have to calculate the median
    for (unsigned int i = 0; i < values->size(); i += 1) {
        median += (*values)[i];
    }
    median = median / static_cast<double> (values->size());

    // with the median we can calculate the standard deviation
    for (unsigned int i = 0; i < values->size(); i += 1) {
        double helper = (*values)[i] - median;
        helper = helper * helper;
        deviation += helper;
    }
    deviation /= static_cast<double> (values->size() - 1);
    deviation = sqrt(deviation);
    return deviation;
}

/**
 *@Description Calculates the median of the values in the given Vector
 *@param  pointer to vector that contains the values(vector<double>)
 *@return  Media value (double)
 */
double globalStatistic::calcMedian(vector<double> *values) {
    if (values->size() == 0)
        return 0;

    double median = 0;
    for (unsigned int i = 0; i < values->size(); i += 1)
        median += (*values)[i];

    return ( median / static_cast<double> (values->size()));
}
