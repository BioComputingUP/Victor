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

#ifndef SUPERIMPOSITOR_H
#define	SUPERIMPOSITOR_H

#include "Protein.h"
#include <KabschMethod.h>
#include <eigen3/Eigen/Geometry>
#include <vector>

namespace Victor {
    namespace Biopool {

        /**
         * @brief Do the superimposition between two proteins using different
         * methods of rotation. Also get back the value for different metrics, 
         * the ranges use for the superimposition for every metrics and the 
         * rotated protein in pdb format.
         * 
         * */
        class SuperImpositor {
        public:

            // CONSTRUCTORS/DESTRUCTOR:
            SuperImpositor(Protein* firstProtein, Protein* secondProtein, string method);
            virtual ~SuperImpositor();

            // PREDICATES:
            void calculateRMSD();
            void calculateMaxSub(double d, std::vector<std::pair<int, int> > vectorSet, char E);
            void calculateGdt(std::vector<std::pair<int, int> > vectorSet);
            void calculateTMScore(std::vector<std::pair<int, int> > vectorSet);



            Spacer* getSet2() const;
            Spacer* getSet1() const;
            Spacer getRMSDset2() const;
            Spacer getRMSDset1() const;
            Spacer getMaxSubset2() const;
            Spacer getMaxSubset1() const;
            Spacer getTMScoreset2() const;
            Spacer getTMScoreset1() const;
            Spacer getGdtset4_2() const;
            Spacer getGdtset4_1() const;
            Spacer getGdtset3_2() const;
            Spacer getGdtset3_1() const;
            Spacer getGdtset2_2() const;
            Spacer getGdtset2_1() const;
            Spacer getGdtset1_2() const;
            Spacer getGdtset1_1() const;
            double getTMScoreValue() const;
            double getGdtValue() const;
            double getMaxsubValue() const;
            double getRmsdValue() const;
            std::vector<std::pair<int, int> > getTMScoreAlignment() const;
            std::vector<std::pair<int, int> > getGdtAlignment4() const;
            std::vector<std::pair<int, int> > getGdtAlignment3() const;
            std::vector<std::pair<int, int> > getGdtAlignment2() const;
            std::vector<std::pair<int, int> > getGdtAlignment1() const;
            std::vector<std::pair<int, int> > getMaxsubAlignment() const;

            //Help function
            static void calculateRotation(Eigen::Matrix3Xd& firstSet, Eigen::Affine3d* rotoTraslation);
        private:

            // PREDICATES:

            double maxEvaluate(double d, std::vector<std::pair<int, int> > vectorSet, char E, Eigen::Matrix3Xd& modifyMatrixSet1, Eigen::Matrix3Xd& modifyMatrixSet2, Eigen::Affine3d*& rotoTraslation, std::vector<std::pair<int, int> >& range);
            std::vector<std::pair<int, int> > maxSubAlignment(Eigen::Matrix3Xd& firstSet, Eigen::Matrix3Xd& secondSet, std::vector< std::pair<int, int> > vectorSet, double d, char E, Eigen::Affine3d*& rotoTraslation);
            std::vector< std::pair<int, int> > Extend(std::vector<std::pair<int, int> > M, std::vector<std::pair<int, int> > vectorSet, Eigen::Matrix3Xd& A, Eigen::Matrix3Xd& B, double d, int L, int n, char E, Eigen::Affine3d*& rotoTraslation);


            //Help function
            vgMatrix3<double> fromMatrix3XdTovgMatrix3(Eigen::Matrix3Xd matrix3Xd) const;
            Spacer fromMatrix3XdToSpacer(Eigen::Matrix3Xd matrix3Xd, int num) const;
            Spacer rotateSpacer(Eigen::Affine3d* rotoTraslation, int num) const;
            vgVector3<double> fromVector3dTovgVector3(Eigen::Vector3d Vector3d) const;
            Eigen::Vector3d fromvgVector3ToVector3d(vgVector3<double> vgVector3) const;
            Eigen::Matrix3Xd fromSpacerToMatrix3Xd(Spacer spacerSet) const;


            // ATTRIBUTES:
            // This is the rotation algorithm chose for this superImpositor
            Rotator* rotationAlgorith;
            //This is the Spacer of the first protein in input
            Spacer* set1;
            //This is the Spacer of the second protein in input
            Spacer* set2;
            //This is the Spacer of the first protein in input for  RMSD method
            Spacer RMSDset1;
            //This is the Spacer of the second protein in input for  RMSD method
            Spacer RMSDset2;
            //This is the Spacer of the first protein in input for  RMSD method
            Spacer MaxSubset1;
            //This is the Spacer of the second protein in input for  RMSD method
            Spacer MaxSubset2;

            //This is the Spacer of the first protein in input for  RMSD method
            Spacer Gdtset1_1;
            //This is the Spacer of the second protein in input for  RMSD method
            Spacer Gdtset1_2;
            //This is the Spacer of the first protein in input for  RMSD method
            Spacer Gdtset2_1;
            //This is the Spacer of the second protein in input for  RMSD method
            Spacer Gdtset2_2;
            //This is the Spacer of the first protein in input for  RMSD method
            Spacer Gdtset3_1;
            //This is the Spacer of the second protein in input for  RMSD method
            Spacer Gdtset3_2;
            //This is the Spacer of the first protein in input for  RMSD method
            Spacer Gdtset4_1;
            //This is the Spacer of the second protein in input for  RMSD method
            Spacer Gdtset4_2;


            //This is the Spacer of the first protein in input for  RMSD method
            Spacer TMScoreset1;
            //This is the Spacer of the second protein in input for  RMSD method
            Spacer TMScoreset2;
            //This is a matrix that contain all the coords of the CA atoms in set1
            Eigen::Matrix3Xd matrixSet1;
            //This is a matrix that contain all the coords of the CA atoms in set2
            Eigen::Matrix3Xd matrixSet2;


            //This is the value of the RMSD between the sequences set1 and set2
            double rmsdValue;
            //This is the value of the maxsub between the sequences set1 and set2
            double maxsubValue;
            //This is the value of the gdt between the sequences set1 and set2
            double gdtValue;
            //This is the value of the TMSCORE between the sequences set1 and set2
            double TMScoreValue;

            //This is an array that contain the alignment between the sequences set1 and set2 using maxsub method
            std::vector<std::pair<int, int> > maxsubAlignment;
            //This is an array that contain the alignment between the sequences set1 and set2 using gdt method
            std::vector<std::pair<int, int> > gdtAlignment1;
            std::vector<std::pair<int, int> > gdtAlignment2;
            std::vector<std::pair<int, int> > gdtAlignment3;
            std::vector<std::pair<int, int> > gdtAlignment4;
            //This is an array that contain the alignment between the sequences set1 and set2 using TMScore method
            std::vector<std::pair<int, int> > TMScoreAlignment;
        };
    }
}

#endif	/* SUPERIMPOSITOR_H */

