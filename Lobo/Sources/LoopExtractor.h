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


#ifndef _LOOPEXTRACTOR_H_
#define _LOOPEXTRACTOR_H_

// Includes:
#include<Spacer.h>
#include<set>
#include<ranking_helper.h>
#include<ranking_helper2.h>

// Global constants, typedefs, etc. (to avoid):
using namespace Victor::Biopool;
namespace Victor { namespace Lobo {

    /**
     * @brief  Extracts all the loop regions (by numbers) from a given spacer.
     * */
    class LoopExtractor {
    public:

        // CONSTRUCTORS/DESTRUCTOR:
        LoopExtractor();
        LoopExtractor(Spacer* s);

        // PREDICATES:
        void nextLoop(int& start, int& end);
        double givePercentProp(multiset<ranking_helper> s, int count);
        double givePercentVdw(multiset<ranking_helper2> s, int count);
        void writeFile(vector<double> prop_percent, vector<double> vdw_percent,
                multiset<ranking_helper2> &rmsh, ofstream &statOut);

        // MODIFIERS:
        void setSpacer(Spacer* s);
        void calcPropPercent(vector<double> &result, multiset<ranking_helper2> &rmsh,
                multiset<ranking_helper> &rhs);
        void calcVdwPercent(vector<double> &result, multiset<ranking_helper2> &rmsh,
                multiset<ranking_helper2> &rhs);
        static const int BEST_COUNT; // determines the number of solutions which 
        // are evaluated by propensities and collisions

        // OPERATORS:

    protected:

    private:

        // HELPERS:
        int getPositionProp(multiset<ranking_helper> s, int index);
        int getPositionVdw(multiset<ranking_helper2> s, int index);

        // ATTRIBUTES:
        Spacer* sp; // contains the spacer which is examined
        int position; // current working position (from where we start 
        // searching for the next loop)
    };

}} // namespace
#endif
