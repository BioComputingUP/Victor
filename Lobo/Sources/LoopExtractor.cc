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
 *@Class:    LoopExtractor               
 * 
 *@Description:Extracts all the loop regions (by numbers) from a given spacer.
 *      
 */

// Includes:
#include "LoopExtractor.h"

// Global constants, typedefs, etc. (to avoid):
//using namespace Biopool;
namespace Biopool {
    const int LoopExtractor::BEST_COUNT = 6;


    // CONSTRUCTORS/DESTRUCTOR:

    /**
     *@Description basic constructor
     */
    LoopExtractor::LoopExtractor() {
        position = 0;
    }

    /**
     *@Description constructor based in a spacer
     *@parameter pointer to the Spacer
     */
    LoopExtractor::LoopExtractor(Spacer* s) {
        sp = s;
        position = 1;
    }

    // PREDICATES:

    /**
     *@Description This routine extracts the next loop out of the spacer and returns the 
     *  indices of the loop in the two variables.The indices start and end contain (counting starts with 0!) the first and
     *  last of the loop: [start, end] is the loop!.If no loop is found, start and end are set to -1
     *@param  reference to the start and to the end position(int&,int&)
     *@return  changes are made internally(void)
     */
    void LoopExtractor::nextLoop(int& start, int& end) {
        bool loop_found = false;
        start = -1;
        end = -1;

        for (unsigned int i = position; i < sp->sizeAmino(); i++) {
            if (loop_found == false) {
                if ((sp->getAmino(i)).getState() == HELIX
                        || (sp->getAmino(i)).getState() == STRAND)
                    continue;
                else { // we found a loop
                    start = i;
                    loop_found = true;
                    continue;
                }
            }
            else { // extend an already found loop

                if ((sp->getAmino(i)).getState() == HELIX
                        || (sp->getAmino(i)).getState() == STRAND) {
                    end = i;
                    position = i + 1; // set new start position
                    return;
                }
                else
                    continue;
            }
        }
        // if we reach this point, we either found no loop or a loop 
        // at the end of the protein:
        if (start == -1) { // we found no loop
            position = sp->sizeAmino();
        } else { // we have a loop at the end
            end = sp->sizeAmino() - 1;
            position = sp->sizeAmino();
        }
    }

    /**
     *@Description The rms, the ranking of the solution regarded to the propensity and 
     *  the VDW-forces is written to the file specified as a parameter.
     *@param vector containing the propensity percentages and van der whals percentages(vector<double>, vector<double>),
     * reference to the set of RMS(multiset<ranking_helper2> &), the reference to the output file(ofstream &)
     *@return  changes are made internally(void)
     */
    void LoopExtractor::writeFile(vector<double> prop_percent, vector<double> vdw_percent,
            multiset<ranking_helper2> &rmsh, ofstream &statOut) {

        int count = 0; // used to count to BEST_COUNT
        set<ranking_helper2>::iterator pos; // used to iterate over the rmsh

        ASSERT(prop_percent.size() == vdw_percent.size(), exception);
        ASSERT(prop_percent.size() >= BEST_COUNT, exception);

        for (pos = rmsh.begin(); pos != rmsh.end(); ++pos) {

            if (count >= BEST_COUNT) {
                break; // we output the values only for the best solutions
            }
            statOut << "RMS der " << count + 1 << " besten Loesung: " << pos->get_value() <<
                    " Prozent VDW: " << vdw_percent[count] << " Prozent Propensities: " << prop_percent[count] << endl;
            count += 1;
        }

    }

    // MODIFIERS:

    /**
     *@Description sets the spacer in the object
     *@param  pointer to  the Spacer to set(Spacer*)
     *@return  changes are made internally(void)
     */
    void LoopExtractor::setSpacer(Spacer* s) {
        sp = s;
    }

    /**
     *@Description Calculates the rank (seen as percent) of the top RMS solutions
     *  in the rhs set (which contains the ranking of the solutions
     *  by propensities).The rank of each of the top solutions regarded to the propensities
     *  is put into the result vector. 
     *@param  reference to the resulting data(vector<double>), reference to the set of 
     * RMS(multiset<ranking_helper2> &), reference to the set of RMS(multiset<ranking_helper> &)
     *@return  changes are made internally(void)
     */
    void LoopExtractor::calcPropPercent(vector<double> &result, multiset<ranking_helper2> &rmsh,
            multiset<ranking_helper> &rhs) {

        int zaehler = 1; // used in a loop in order to calculate the propensities 
        // for the lowest rms solutions
        set<ranking_helper2>::iterator pos; // used to iterate over the rmsh

        for (pos = rmsh.begin(); pos != rmsh.end(); ++pos) {
            double helper;
            if (zaehler > BEST_COUNT) {
                break; // we calculate these values only for the best solutions
            }
            helper = givePercentProp(rhs, (*pos).get_index() + 1);
            result.push_back(helper);
            cout << "Index:" << pos->get_index() << " Value (RMS) " << pos->get_value()
                    << " Prozentsatz fuer Propensities: " << helper << endl;
            zaehler += 1;
        }
    }

    /**
     *@Description Calculates the rank (seen as percent) of the top RMS solutions
     *  in the rhs set (which contains the ranking of the solutions
     *  by rms).The rank of each of the top solutions regarded to the rms
     *  is put into the result vector. 
     *@param  reference to the resulting data(vector<double>), reference to the set of
     *  RMS(multiset<ranking_helper2> &), reference to the set of RMS(multiset<ranking_helper> &)
     *@return  changes are made internally(void)
     */
    void LoopExtractor::calcVdwPercent(vector<double> &result, multiset<ranking_helper2> &rmsh,
            multiset<ranking_helper2> &rhs) {

        int zaehler = 1; // used in a loop in order to calculate the vdw
        // for the lowest rms solutions
        set<ranking_helper2>::iterator pos; // used to iterate over the rmsh

        for (pos = rmsh.begin(); pos != rmsh.end(); ++pos) {
            double helper;
            if (zaehler > BEST_COUNT) {
                break; // we calculate these values only for the best solutions
            }
            helper = givePercentVdw(rhs, (*pos).get_index() + 1);
            result.push_back(helper);
            cout << "Index:" << pos->get_index() << " Value (RMS) " << pos->get_value()
                    << " Prozentsatz fuer VDW: " << helper << endl;
            zaehler += 1;
        }
        cout << endl;
    }


    // OPERATORS:

    // HELPERS:

    /**
     *@Description This function calculates for the solution in the given vector with
     *    the entry index, what position it has in the given sorted set. 
     *@param  reference to the set of RMS(multiset<ranking_helper> &), index(int)
     *@return  corresponding propensity percentage(double)
     */
    double LoopExtractor::givePercentProp(multiset<ranking_helper> s, int count) {

        int pos = 0; // contains the position of the current solution in the set
        int size = s.size(); // contains the size of the set
        ASSERT(size > 0, exception);
        double helper; // used to store the percent value for each of the solutions

        pos = getPositionProp(s, count);
        ASSERT(pos >= 1, exception);
        helper = static_cast<double> (pos) / static_cast<double> (size);
        return helper;
    }

    /**
     *@Description This function calculates for the solution in the given vector with
     *    the entry index, what position it has in the given sorted set.
     *
     *    Uses the function getPositionVdw()
     *@param  reference to the set of RMS(multiset<ranking_helper2> &), index(int)
     *@return  corresponding  Van der Wals percentage (double)
     */
    double LoopExtractor::givePercentVdw(multiset<ranking_helper2> s, int count) {

        int pos = 0; // contains the position of the current solution in the set
        int size = s.size(); // contains the size of the set
        ASSERT(size > 0, exception);
        double helper; // used to store the percent value for each of the solutions

        pos = getPositionVdw(s, count);
        ASSERT(pos >= 1, exception);
        helper = static_cast<double> (pos) / static_cast<double> (size);
        return helper;
    }

    /**
     *@Description This function returns the position of the entry in the set which has in
     *   its index field  the entry index. It returns -1 if the entry is not in
     *   in the set. The counting starts with 1.
     *@param  reference to the set of RMS(multiset<ranking_helper> &), index(int)
     *@return  corresponding propensity position in the multiset(int)
     */
    int LoopExtractor::getPositionProp(multiset<ranking_helper> s, int index) {

        int counter = 1; // contains the position of index and thus the return value

        set<ranking_helper>::iterator pos;
        for (pos = s.begin(); pos != s.end(); ++pos) {
            if ((*pos).get_index() == index) {
                return counter;
            }
            counter += 1;
        }
        return -1;
    }

    /**
     *@Description This function returns the position of the entry in the set which has in
     *   its index field  the entry index. It returns -1 if the entry is not in
     *   in the set. The counting starts with 1.
     *@param  reference to the set of RMS(multiset<ranking_helper2> &), index(int)
     *@return  corresponding Van Der Wals position in the multiset(int)
     */
    int LoopExtractor::getPositionVdw(multiset<ranking_helper2> s, int index) {

        int counter = 1; // contains the position of index and thus the return value

        set<ranking_helper2>::iterator pos;
        for (pos = s.begin(); pos != s.end(); ++pos) {
            if ((*pos).get_index() == index) {
                return counter;
            }
            counter += 1;
        }
        return -1;
    }
}
