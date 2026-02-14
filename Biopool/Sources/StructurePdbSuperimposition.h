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
 * Author: Simone Piano
 *
 */

#ifndef STRUCTUREPDBSUPERIMPOSITION_H
#define	STRUCTUREPDBSUPERIMPOSITION_H

#include <Eigen/Geometry>

namespace Victor {
    namespace Biopool {
        class StructurePdbSuperimposition {            
        public:
            //CONSTRUCTOR/DESTRUCTOR
            StructurePdbSuperimposition();
            ~StructurePdbSuperimposition();
            
            //ROTOTRANSLATION METHOD
            static std::pair<Eigen::Matrix3d,Eigen::Vector3d> KabschAlgorithm(Eigen::MatrixX3d&, Eigen::MatrixX3d&);
            
            //ITERATIVE SEARCH METHODS
            static std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d> MaxSub(Eigen::MatrixX3d&, Eigen::MatrixX3d&, Spacer&, Spacer&, int, double);
            static void MaxSubExtend(std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d>&, Eigen::MatrixX3d&, Eigen::MatrixX3d&, double);
            
            //SCORE METHODS
            static double RMSD(Eigen::MatrixX3d&, Eigen::MatrixX3d&, Spacer&, Spacer&);
            static double GDTTS(Eigen::MatrixX3d&, Eigen::MatrixX3d&, Spacer&, Spacer&, int, double[]);
            static double TMScore(Eigen::MatrixX3d&, Eigen::MatrixX3d&, Spacer&, Spacer&, int);
            static double MaxSubScore(Eigen::MatrixX3d&, Eigen::MatrixX3d&, Spacer&, Spacer&, int, double, bool);
            
            //SUPPORT METHODS
            static Eigen::MatrixX3d fromSpacerToEigenMatrix(Spacer*);
            static void applyMatrixRototranslation(Eigen::MatrixX3d&, Eigen::MatrixX3d&);
            static void applyVectorRototranslation(Eigen::Vector3d&);
            static void removeRow(Eigen::MatrixX3d&, int);
            
            //OTHER METHODS
            static void savePdbFile(Spacer&, Spacer&, string);
            static Spacer updateSpacer(Spacer&);                    
            
        protected:
            
        private:
            
        };
    }
}
#endif
