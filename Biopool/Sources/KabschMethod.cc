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

#include <KabschMethod.h>
#include <eigen3/Eigen/Geometry>

using namespace Victor;
using namespace Victor::Biopool;

// CONSTRUCTORS:

KabschMethod::KabschMethod() {
}
/**
 * This method implement the KabschMethod method for obtein the rototraslation 
 * that minimize the rmsd value of the rotated protein and the second protein.
 * 
 * @param set1Matrix (Eigen::Matrix3Xd), the 3*N matrix with the coordinate
 * of the atom of the first protein;
 * @param set2Matrix (Eigen::Matrix3Xd), the 3*N matrix with the coordinate
 * of the atom of the second protein;
 * @return Eigen::Affine3d*, the rototraslation returned by the method.
 */
Eigen::Affine3d* KabschMethod::rotate(Eigen::Matrix3Xd set1Matrix, Eigen::Matrix3Xd set2Matrix) const{
    Eigen::Affine3d* outputCoords = new Eigen::Affine3d();
    //Rotation Matrix
    outputCoords->linear() = Eigen::Matrix3d::Identity(3, 3);
    //Translation Vector
    outputCoords->translation() = Eigen::Vector3d::Zero();

    //Check that the two set have the same length
    if (set1Matrix.cols() != set2Matrix.cols())
        throw "I due set devono presentare uguale lunghezza";


    // Calculate the distance between the consecutive points in the matrix
    double set1DistancesSum = 0, set2DistancesSum = 0;
    for (int col = 0; col < set1Matrix.cols() - 1; col++) {
        set1DistancesSum += (set1Matrix.col(col + 1) - set1Matrix.col(col)).norm();
        set2DistancesSum += (set2Matrix.col(col + 1) - set2Matrix.col(col)).norm();
    }

    //If all the point is equal simply return outputCoords
    if (set1DistancesSum <= 0 || set2DistancesSum <= 0)
        return outputCoords;


    // Find the centroids and then shift to the origin
    Eigen::Vector3d set1Centroids = Eigen::Vector3d::Zero();
    Eigen::Vector3d set2Centroids = Eigen::Vector3d::Zero();
    for (int col = 0; col < set1Matrix.cols(); col++) {
        set1Centroids += set1Matrix.col(col);
        set2Centroids += set2Matrix.col(col);
    }
    set1Centroids /= set1Matrix.cols();
    set2Centroids /= set2Matrix.cols();
    for (int col = 0; col < set1Matrix.cols(); col++) {
        set1Matrix.col(col) -= set1Centroids;
        set2Matrix.col(col) -= set2Centroids;
    }
    //FINE

    // Calculate svd decomposition
    Eigen::MatrixXd Cov = set1Matrix * set2Matrix.transpose();
    // Modification of the matrix for a more efficient
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Define the direction of rotation
    double direction = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (direction > 0)
        direction = 1.0;
    else
        direction = -1.0;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    //Change the last value of rotation matrix to use the correct direction of rotation
    I(2, 2) = direction;
    
    //Calculation of the rotation matrix
    Eigen::Matrix3d rotationMatrix = svd.matrixV() * I * svd.matrixU().transpose();
    outputCoords->linear() = rotationMatrix;
    
    //The translation is calculate has the distance between the centroids of the
    //second set and the first rototrasled
    outputCoords->translation() = set2Centroids - rotationMatrix * set1Centroids;
    return outputCoords;
}

