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

#ifndef KABSCHMETHOD_H
#define KABSCHMETHOD_H

#include <Rotator.h>


namespace Victor {
    namespace Biopool {

        /**@brief Implementation of the kabsch method for find the optimal 
         * rotoslation for superimpose two molecules.
         * 
         * */
        class KabschMethod:public Rotator{
        public:

            // CONSTRUCTORS/DESTRUCTOR:
            KabschMethod();
            
            Eigen::Affine3d* rotate(Eigen::Matrix3Xd set1Matrix, Eigen::Matrix3Xd set2Matrix) const;

        };

    }
} // namespace

#endif
