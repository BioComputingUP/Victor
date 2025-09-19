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
 * Author: Stefano Cecconello
 *
 */

#ifndef ROTATOR_H
#define ROTATOR_H

#include <Spacer.h>
#include <eigen3/Eigen/Geometry>

namespace Victor {
    namespace Biopool {

        /**@brief Abstract class for the different types of rotation methods.
         * 
         * */
        class Rotator {
        public:
            // CONSTRUCTORS/DESTRUCTOR:
            Rotator();

            // PREDICATES:
            virtual Eigen::Affine3d* rotate(Eigen::Matrix3Xd set1Matrix, Eigen::Matrix3Xd set2Matrix) const = 0;



        };

    }
} // namespace

#endif
