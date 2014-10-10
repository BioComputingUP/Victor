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


#ifndef _RANK_HELPER2_H_
#define _RANK_HELPER2_H_

// Includes:

// Global constants, typedefs, etc. (to avoid):

namespace Victor { namespace Lobo {

    /**
     * @brief the class contains methods to manage the ranking
     * 
     * @Description  JJust contains an int and double value which gets sorted (by a STL set)
     *    in order to determine the ranking of low rms solutions with special 
     *    filters. This container differs from the ranking_helper by the fact,
     *    that it has the < operator overloaded differently 
     * */
    class ranking_helper2 {
    public:

        ranking_helper2(int ind, double val);
        ranking_helper2(const ranking_helper2& c);
        int get_index() const;
        double get_value() const;
        bool operator<(const ranking_helper2 &name) const;
        ranking_helper2& operator=(const ranking_helper2& orig);
        void copy(const ranking_helper2& c);

    private:
        int index;
        double value;
    };

    /**
     *@Description returns the rms ranking of the solution
     *@param  none
     *@return  the corresponding value
     */
    inline int ranking_helper2::get_index() const {
        return index;
    }

    /**
     *@Description returns a value from a filter (like the propensity or the collision)
     *@param  none
     *@return  the corresponding value
     */
    inline double ranking_helper2::get_value() const {
        return value;
    }


}} // namespace
#endif
