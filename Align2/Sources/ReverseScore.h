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


#ifndef __ReverseScore_H__
#define __ReverseScore_H__

#include <Align.h>
#include <Alignment.h>
#include <StatTools.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace Victor { namespace Align2{

    /** @brief  
     * 
     *   

     **/
    class ReverseScore {
    public:

        // CONSTRUCTORS:

        /// Default constructor.
        ReverseScore(Align *a);

        /// Copy constructor.
        ReverseScore(const ReverseScore &orig);

        /// Destructor
        virtual ~ReverseScore();


        // OPERATORS:

        /// Assignment operator.
        ReverseScore& operator =(const ReverseScore &orig);


        // MODIFIERS:

        /// Copy orig object to this object ("deep copy").
        virtual void copy(const ReverseScore &orig);

        /// Calculate Z-score between ali and its reversal.
        double getZScore(double &forward, double &reverse, unsigned int n = 50);


    protected:

        // ATTRIBUTES:

        Align *ali; ///< Pointer to initial Align.
        Align *inv; ///< Pointer to inverted Align.


    private:

    };

}} // namespace

#endif
