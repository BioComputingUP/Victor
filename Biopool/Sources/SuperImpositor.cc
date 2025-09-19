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

#include "SuperImpositor.h"
#include <math.h>
#include <vector>

using namespace Victor;
using namespace Victor::Biopool;

/**
 * This is the constructor for the superImpositor object. It sets the indicated 
 * method for the rotation and save the proteins where he has to work on.
 *  
 * @param prot1 (Protein*) , the first protein in input;
 * @param prot2 (Protein*), the second protein in input;
 * @param method (string), a string with the name of the rotation method used. 
 */
SuperImpositor::SuperImpositor(Protein* firstProtein, Protein* secondProtein, string method = "") {
    if (method == "Kabsch") {
        rotationAlgorith = new KabschMethod();
    }



    //set default rotation method
    if (method == "") {
        rotationAlgorith = new KabschMethod();
    }
    //If the input method does not exist then the default method is set
    if (rotationAlgorith) {
        rotationAlgorith = new KabschMethod();
    }

    //Get the spacers
    set1 = firstProtein->getSpacer((unsigned int) 0);
    set2 = secondProtein->getSpacer((unsigned int) 0);
    matrixSet1 = fromSpacerToMatrix3Xd(*(set1));
    matrixSet2 = fromSpacerToMatrix3Xd(*(set2));
    rmsdValue = 999;
    maxsubValue = 999;
    gdtValue = 999;
    TMScoreValue = 999;
}

SuperImpositor::~SuperImpositor() {
}

/**
 *  Calculates RMSD value between the two proteins given in the constructor.
 */
void SuperImpositor::calculateRMSD() {
    if (rmsdValue == 999) {
        Eigen::Matrix3Xd modifyMatrixSet1 = matrixSet1;
        Eigen::Matrix3Xd modifyMatrixSet2 = matrixSet2;
        Eigen::Affine3d* rotoTraslation = rotationAlgorith->rotate(modifyMatrixSet1, modifyMatrixSet2);
        calculateRotation(modifyMatrixSet1, rotoTraslation);
        RMSDset1 = rotateSpacer(rotoTraslation, 1);
        //RMSDset2 = rotateSpacer(rotoTraslation, 2);
        //RMSDset1 = fromMatrix3XdToSpacer(modifyMatrixSet1, 1);
        RMSDset2 = fromMatrix3XdToSpacer(modifyMatrixSet2, 2);
        //Calculate rmsd
        double distanceSquared = 0;
        double msd;
        int columns = modifyMatrixSet1.cols();
        for (int col = 0; col < columns; col++) {
            distanceSquared += (modifyMatrixSet1.col(col) - modifyMatrixSet2.col(col)).squaredNorm();
        }
        msd = distanceSquared / columns;
        rmsdValue = sqrt(msd);
    }
}

/**
 * Calculates MaxSub value between the two proteins given in the constructor.
 * 
 * @param d (double), the threshold for the maximum distance, between two atoms,
 * accepted during the research;
 * @param vectorSet (std::vector<std::pair<int, int> >), the vector that contain
 * the couples of aligned position of the two structures in input;
 * @param E (char), a parameter for decide if the atoms over the distance d have 
 * to be deleted in the research of the best structure.
 */
void SuperImpositor::calculateMaxSub(double d, std::vector<std::pair<int, int> > vectorSet, char E) {
    if (maxsubValue == 999) {
        Eigen::Matrix3Xd modifyMatrixSet1 = matrixSet1;
        Eigen::Matrix3Xd modifyMatrixSet2 = matrixSet2;
        double maxSub;
        Eigen::Affine3d* rotoTraslation;
        std::vector<std::pair<int, int> > range;
        maxSub = maxEvaluate(d, vectorSet, E, modifyMatrixSet1, modifyMatrixSet2, rotoTraslation, range);
        maxsubAlignment = range;
        MaxSubset1 = rotateSpacer(rotoTraslation, 1);
        MaxSubset2 = fromMatrix3XdToSpacer(modifyMatrixSet2, 2);
        maxsubValue = maxSub;
    }
}

/**
 * Calculates Gdt value between the two proteins given in the constructor.
 * 
 * @param vectorSet (std::vector<std::pair<int, int> >), the vector that contain
 * the couples of aligned position of the two structures in input.
 */
void SuperImpositor::calculateGdt(std::vector<std::pair<int, int> > vectorSet) {
    if (gdtValue == 999) {
        //0.6135 con maxsub programmino per 4 e 5
        Eigen::Matrix3Xd modifyMatrixSet1 = matrixSet1;
        Eigen::Matrix3Xd modifyMatrixSet2 = matrixSet2;
        Eigen::Matrix3Xd modifyMatrixSet3 = matrixSet1;
        Eigen::Matrix3Xd modifyMatrixSet4 = matrixSet2;
        Eigen::Matrix3Xd modifyMatrixSet5 = matrixSet1;
        Eigen::Matrix3Xd modifyMatrixSet6 = matrixSet2;
        Eigen::Matrix3Xd modifyMatrixSet7 = matrixSet1;
        Eigen::Matrix3Xd modifyMatrixSet8 = matrixSet2;
        Eigen::Affine3d* rotoTraslation1;
        Eigen::Affine3d* rotoTraslation2;
        Eigen::Affine3d* rotoTraslation3;
        Eigen::Affine3d* rotoTraslation4;
        double sum = 0;
        std::vector<std::pair<int, int> > range1;
        std::vector<std::pair<int, int> > range2;
        std::vector<std::pair<int, int> > range3;
        std::vector<std::pair<int, int> > range4;
        sum = sum + maxEvaluate(1, vectorSet, 'n', modifyMatrixSet1, modifyMatrixSet2, rotoTraslation1, range1);
        sum = sum + maxEvaluate(2, vectorSet, 'n', modifyMatrixSet3, modifyMatrixSet4, rotoTraslation2, range2);
        sum = sum + maxEvaluate(4, vectorSet, 'n', modifyMatrixSet5, modifyMatrixSet6, rotoTraslation3, range3);
        sum = sum + maxEvaluate(8, vectorSet, 'n', modifyMatrixSet7, modifyMatrixSet8, rotoTraslation4, range4);
        gdtAlignment1 = range1;
        gdtAlignment2 = range2;
        gdtAlignment3 = range3;
        gdtAlignment4 = range4;

        sum = sum / 4;

        Gdtset1_1 = rotateSpacer(rotoTraslation1, 1);
        Gdtset1_2 = fromMatrix3XdToSpacer(modifyMatrixSet2, 2);
        Gdtset2_1 = rotateSpacer(rotoTraslation2, 1);
        Gdtset2_2 = fromMatrix3XdToSpacer(modifyMatrixSet4, 2);
        Gdtset3_1 = rotateSpacer(rotoTraslation3, 1);
        Gdtset3_2 = fromMatrix3XdToSpacer(modifyMatrixSet6, 2);
        Gdtset4_1 = rotateSpacer(rotoTraslation4, 1);
        Gdtset4_2 = fromMatrix3XdToSpacer(modifyMatrixSet8, 2);
        gdtValue = sum;
    }
}

/**
 * Calculates TMScore value between the two proteins given in the constructor.
 * 
 * @param vectorSet (std::vector<std::pair<int, int> >), the vector that contain
 * the couples of aligned position of the two structures in input.
 */
void SuperImpositor::calculateTMScore(std::vector<std::pair<int, int> > vectorSet) {
    if (TMScoreValue == 999) {
        Eigen::Matrix3Xd modifyMatrixSet1 = matrixSet1;
        Eigen::Matrix3Xd modifyMatrixSet2 = matrixSet2;
        double d0;
        d0 = (1.24 * cbrt(modifyMatrixSet1.cols() - 15)) - 1.8;
        double TMScore;
        Eigen::Affine3d* rotoTraslation;
        std::vector<std::pair<int, int> > range;
        TMScore = maxEvaluate(d0, vectorSet, 'n', modifyMatrixSet1, modifyMatrixSet2, rotoTraslation, range);
        TMScoreAlignment = range;
        TMScoreset1 = rotateSpacer(rotoTraslation, 1);
        TMScoreset2 = fromMatrix3XdToSpacer(modifyMatrixSet2, 2);
        TMScoreValue = TMScore;
    }
}

/**
 * Calculates a value for indicate the quality of the superimposition found.
 * 
 * @param d (double), the threshold for the maximum distance, between two atoms,
 * accepted during the research
 * @param vectorSet (std::vector<std::pair<int, int> >), the vector that contain
 * the couples of aligned position of the two structures in input;
 * @param E (char), a parameter for decide if the atoms over the distance d have 
 * to be deleted in the research of the best structure;
 * @param modifyMatrixSet1 (Eigen::Matrix3Xd&), the 3*N matrix with the coordinate
 * of the atom of the first protein;
 * @param modifyMatrixSet2 (Eigen::Matrix3Xd&), the 3*N matrix with the coordinate
 * of the atom of the second protein;
 * @param rotoTraslation (Eigen::Affine3d*&), a reference for return the final 
 * rotation matrix and translation vector;
 * @param range (std::vector<std::pair<int, int> >&) a reference for return the 
 * couples of the best superimposition.
 * 
 * @return double a value between 0 and 1 for the best superimposition found 
 * with the given parameter.
 */
double SuperImpositor::maxEvaluate(double d, std::vector<std::pair<int, int> > vectorSet, char E, Eigen::Matrix3Xd& modifyMatrixSet1, Eigen::Matrix3Xd& modifyMatrixSet2, Eigen::Affine3d*& rotoTraslation, std::vector<std::pair<int, int> >& range) {
    double distance;
    double sum = 0;


    range = maxSubAlignment(modifyMatrixSet1, modifyMatrixSet2, vectorSet, d, E, rotoTraslation);
    for (unsigned int i = 0; i < range.size(); i++) {
        distance = (modifyMatrixSet1.col(range[i].first) - modifyMatrixSet2.col(range[i].second)).squaredNorm();
        distance = sqrt(distance);
        sum = sum + (1 / (1 + pow((distance / d), 2)));
    }
    return sum / modifyMatrixSet1.cols();
}

/**
 * The algorithm for find the best superimposition between he two input protein.
 * 
 * @param firstSet (Eigen::Matrix3Xd&), the 3*N matrix with the coordinate
 * of the atom of the first protein;
 * @param secondSet (Eigen::Matrix3Xd&), the 3*N matrix with the coordinate
 * of the atom of the second protein;
 * @param vectorSet (std::vector<std::pair<int, int> >), the vector that contain
 * the couples of aligned position of the two structures in input;
 * @param d (double), the threshold for the maximum distance, between two atoms,
 * accepted during the research
 * @param E (char), a parameter for decide if the atoms over the distance d have 
 * to be deleted in the research of the best structure;
 * @param rotoTraslation (Eigen::Affine3d*&), a reference for return the final 
 * rotation matrix and translation vector.
 * 
 * @return std::vector<std::pair<int, int> > the couples representing the finally best
 * superimposition found.
 */
std::vector<std::pair<int, int> > SuperImpositor::maxSubAlignment(Eigen::Matrix3Xd& firstSet, Eigen::Matrix3Xd& secondSet, std::vector<std::pair<int, int> > vectorSet, double d, char E, Eigen::Affine3d*& rotoTraslation) {
    long unsigned int sMax = 0;
    int n = firstSet.cols();
    int L = 4;
    Eigen::Matrix3Xd modifiedFirstSet;
    Eigen::Matrix3Xd modifiedSecondSet;
    Eigen::Matrix3Xd optFirstSet;
    Eigen::Matrix3Xd optSecondSet;
    Eigen::Affine3d* optRotoTraslation;
    Eigen::Affine3d* modifyRotoTraslation;

    std::vector<std::pair<int, int> > mMax;
    std::vector<std::pair<int, int> > M;
    for (int i = 0; i < n - L + 1; i++) {
        modifiedFirstSet = firstSet;
        modifiedSecondSet = secondSet;
        M.clear();
        //Add the L elements
        for (int j = 0; j < L; j++) {
            M.push_back(vectorSet[i + j]);
        }
        M = Extend(M, vectorSet, modifiedFirstSet, modifiedSecondSet, d, L, n, E, modifyRotoTraslation);


        if (M.size() > sMax) {
            optFirstSet = modifiedFirstSet;
            optSecondSet = modifiedSecondSet;
            optRotoTraslation = modifyRotoTraslation;
            sMax = M.size();
            mMax = M;
        }
    }
    firstSet = optFirstSet;
    secondSet = optSecondSet;
    rotoTraslation = optRotoTraslation;

    return mMax;
}

/**
 * The algorithm for find the best superimposition between the two input protein.
 * 
 * @param M (std::vector<std::pair<int, int> >), a reference for return the actually
 * vector containing a subset of vectorSet;
 * @param vectorSet (std::vector<std::pair<int, int> >), the vector that contain
 * the couples of aligned position of the two structures in input;
 * @param A (Eigen::Matrix3Xd&), the 3*N matrix with the coordinate
 * of the atom of the first protein;
 * @param B (Eigen::Matrix3Xd&), the 3*N matrix with the coordinate
 * of the atom of the second protein;
 * @param d (double), the threshold for the maximum distance, between two atoms,
 * accepted during the research
 * @param L (double), the number of minimum amino acids that the algorithm need to
 * found;
 * @param n (double), the number of atoms in the proteins;
 * @param E (char), a parameter for decide if the atoms over the distance d have 
 * to be deleted in the research of the best structure;
 * @param rotoTraslation (Eigen::Affine3d*&), a reference for return the final 
 * rotation matrix and translation vector.
 * 
 * @return std::vector<std::pair<int, int> > the couples representing the actually best
 * superimposition found.
 */
std::vector< std::pair<int, int> > SuperImpositor::Extend(std::vector<std::pair<int, int> > M, std::vector<std::pair<int, int> > vectorSet, Eigen::Matrix3Xd& A, Eigen::Matrix3Xd& B, double d, int L, int n, char E, Eigen::Affine3d*& rotoTraslation) {
    std::vector<std::pair<int, int> > N;
    double threshold;
    int k = 4;
    double distance;
    Eigen::Matrix3Xd M1(3, L);
    Eigen::Matrix3Xd M2(3, L);

    for (int i = 0; i < L; i++) {
        M1.col(i) = A.col(M[i].first);
        M2.col(i) = B.col(M[i].second);
    }

    Eigen::Matrix3Xd ARototrasled;
    Eigen::Matrix3Xd BRototrasled;


    for (int j = 1; j <= k; j++) {
        ARototrasled = A;
        BRototrasled = B;
        rotoTraslation = rotationAlgorith->rotate(M1, M2);
        calculateRotation(ARototrasled, rotoTraslation);

        //Calculate the distance between the points

        N.clear();

        for (int i = 0; i < n; i++) {
            distance = (ARototrasled.col(vectorSet[i].first) - BRototrasled.col(vectorSet[i].second)).squaredNorm();
            distance = sqrt(distance);
            threshold = ((j * d) / k);
            if (distance <= threshold) {
                N.push_back(vectorSet[i]);
            }

        }


        Eigen::Matrix3Xd N1(3, N.size());
        Eigen::Matrix3Xd N2(3, N.size());
        for (unsigned int i = 0; i < N.size(); i++) {
            N1.col(i) = A.col(N[i].first);
            N2.col(i) = B.col(N[i].second);
        }
        M1 = N1;
        M2 = N2;
    }

    rotoTraslation = rotationAlgorith->rotate(M1, M2);

    calculateRotation(A, rotoTraslation);

    //Calculate the distance between the points
    M = N;
    for (std::vector<std::pair<int, int> >::iterator it = M.begin(); it != M.end(); ++it) {
        distance = (A.col(it->first) - B.col(it->second)).squaredNorm();
        distance = sqrt(distance);
        if (distance > d && E == 'y') {
            M.erase(it);
            if (it == M.end()) {
                break;
            }
        }
    }


    return M;
}

/**
 * This method modifies the input matrix applying on this the rotation and the 
 * translation.
 * 
 * @param firstSet (Eigen::Matrix3Xd&), the 3*N matrix with the coordinate
 * of the atom of the first protein;
 * @param rotoTraslation (Eigen::Affine3d*&), a reference for return the final 
 * rotation matrix and translation vector.
 * 
 */
void SuperImpositor::calculateRotation(Eigen::Matrix3Xd& firstSet, Eigen::Affine3d* rotoTraslation) {
    Eigen::Matrix3Xd R = rotoTraslation->linear();
    Eigen::Vector3d S = rotoTraslation->translation();
    //Apply the rototraslation
    int NumAmino = firstSet.cols();
    Eigen::Matrix3Xd rotoTraslSet(3, NumAmino);
    for (int col = 0; col < NumAmino; col++) {
        rotoTraslSet.col(col) = R * firstSet.col(col) + S;
    }
    //Save rototraslation
    firstSet = rotoTraslSet;
}

/**
 * This is a converting method for get a Matrix3Xd from a spacer. The final matrix
 * content is the 3D coordinates originally present in the spacer.
 * 
 * @param spacerSet (spacer), the spacer that need to be converted.
 * 
 * @return Eigen::Matrix3Xd, the output Matrix3Xd.
 */
Eigen::Matrix3Xd SuperImpositor::fromSpacerToMatrix3Xd(Spacer spacerSet) const {
    int NumAmino = (int) spacerSet.sizeAmino();
    Eigen::Matrix3Xd matrixSet(3, NumAmino);


    for (int col = 0; col < NumAmino; col++) {
        vgVector3<double> coords = spacerSet.getAmino(col)[CA].getCoords();
        matrixSet(0, col) = coords[0];
        matrixSet(1, col) = coords[1];
        matrixSet(2, col) = coords[2];
    }

    return matrixSet;
}

/**
 * This is a converting method for get a vgMatrix3 from a Matrix3Xd. The final matrix
 * content is the 3D coordinates originally present in the input matrix. The method
 * take in input 3*3 matrices.
 * 
 * @param matrix3Xd (Eigen::Matrix3Xd), the matrix that need to be converted.
 * 
 * @return vgMatrix3<double>, the output vgMatrix3.
 */
vgMatrix3<double> SuperImpositor::fromMatrix3XdTovgMatrix3(Eigen::Matrix3Xd matrix3Xd) const {
    //Conversion from Matrix3Xd to vgMatrix3
    vgMatrix3<double> rotationMatrix;
    rotationMatrix.x.x = matrix3Xd(0, 0);
    rotationMatrix.x.y = matrix3Xd(0, 1);
    rotationMatrix.x.z = matrix3Xd(0, 2);
    rotationMatrix.y.x = matrix3Xd(1, 0);
    rotationMatrix.y.y = matrix3Xd(1, 1);
    rotationMatrix.y.z = matrix3Xd(1, 2);
    rotationMatrix.z.x = matrix3Xd(2, 0);
    rotationMatrix.z.y = matrix3Xd(2, 1);
    rotationMatrix.z.z = matrix3Xd(2, 2);

    return rotationMatrix;
}

/**
 * This is a converting method for get a spacer from a Matrix3Xd. The final spacer
 * content is the 3D coordinates originally present in the Matrix3Xd. The dimension
 * of the input matrix and the output matrix are the same of the spacer given in input 
 * to the constructor.
 * 
 * @param matrix3Xd (Eigen::Matrix3Xd matrix3Xd), the matrix that need to be converted;
 * @param num (int), the number of the spacer given in input to the constructor 
 * that need to be used how base for the new spacer.
 * 
 * @return spacer, the output spacer.
 */
Spacer SuperImpositor::fromMatrix3XdToSpacer(Eigen::Matrix3Xd matrix3Xd, int num) const {
    //Change atoms coordinates of set1
    Spacer newSpacer;
    if (num == 1)
        newSpacer = *(set1);
    else
        newSpacer = *(set2);
    int NumAmino = matrix3Xd.cols();
    for (int i = 0; i < NumAmino; i++) {
        newSpacer.getAmino(i)[CA].setCoords(fromVector3dTovgVector3(matrix3Xd.col(i)));
    }
    return newSpacer;
}

/**
 * This method apply the rototraslation given in input to the original spacers given in
 * input to the constructor. The spacer is selected by the second parameter.
 * 
 * @param otoTraslation (Eigen::Affine3d*d), the rotation matrix and the translation 
 * vector;
 * @param num (int), the number of the spacer given in input to the constructor 
 * that need to be used how base for the new spacer.
 * 
 * @return spacer, the output spacer.
 */
Spacer SuperImpositor::rotateSpacer(Eigen::Affine3d* rotoTraslation, int num) const {
    //Change atoms coordinates of set1
    Spacer newSpacer1;
    Spacer newSpacer2;
    Spacer newModifySpacer;
    if (num == 1) {
        newSpacer1 = *(set1);
        newSpacer2 = *(set1);
        newModifySpacer = *(set1);
    } else {
        newSpacer1 = *(set2);
        newSpacer2 = *(set2);
        newModifySpacer = *(set2);
    }
    int NumAmino = newSpacer1.sizeAmino();
    Eigen::Matrix3Xd R = rotoTraslation->linear();
    Eigen::Vector3d S = rotoTraslation->translation();


    vector<Atom> atoms;
    unsigned int contAtom;
    Eigen::Vector3d coords;
    Eigen::Vector3d newCoords;



    for (int i = 0; i < NumAmino; i++) {
        atoms = newSpacer1.getAmino(i).giveAtoms();
        contAtom = atoms.size();
        AminoAcid& AA = newModifySpacer.getAmino(i);
        for (unsigned int j = 0; j < contAtom; j++) {
            coords = fromvgVector3ToVector3d(atoms[j].getCoords());
            newCoords = R * coords + S;
            AA[j].setCoords(fromVector3dTovgVector3(newCoords));
        }

        atoms = newSpacer2.getAmino(i).getSideChain().giveAtoms();

        for (unsigned int j = 0; j < atoms.size(); j++) {
            coords = fromvgVector3ToVector3d(atoms[j].getCoords());
            newCoords = R * coords + S;
            AA[j + contAtom].setCoords(fromVector3dTovgVector3(newCoords));
        }
    }
    return newModifySpacer;
}

/**
 * This is a converting method for get a vgVector3 from a Eigen::Vector3d<double>. 
 * These are 3*1 vectors.
 * 
 * @param Vector3d (Eigen::Vector3d), the vector that need to be converted.
 * 
 * @return vgVector3, the output vector.
 */
vgVector3<double> SuperImpositor::fromVector3dTovgVector3(Eigen::Vector3d Vector3d) const {
    vgVector3<double> traslationVector;
    for (int r = 0; r < 3; r++) {
        traslationVector[r] = Vector3d(r);
    }
    return traslationVector;
}

/**
 * This is a converting method for get a Eigen::Vector3d from a vgVector3<double>. 
 * These are 3*1 vectors.
 * 
 * @param vgVector3 (vgVector3<double>), the vector that need to be converted.
 * 
 * @return Eigen::Vector3d, the output vector.
 */
Eigen::Vector3d SuperImpositor::fromvgVector3ToVector3d(vgVector3<double> vgVector3) const {
    Eigen::Vector3d newVector;
    newVector(0) = vgVector3[0];
    newVector(1) = vgVector3[1];
    newVector(2) = vgVector3[2];
    return newVector;
}

Spacer* SuperImpositor::getSet2() const {

    return set2;
}

Spacer* SuperImpositor::getSet1() const {

    return set1;
}

Spacer SuperImpositor::getRMSDset2() const {

    return RMSDset2;
}

Spacer SuperImpositor::getRMSDset1() const {
    return RMSDset1;
}

Spacer SuperImpositor::getMaxSubset2() const {
    return MaxSubset2;
}

Spacer SuperImpositor::getMaxSubset1() const {
    return MaxSubset1;
}

Spacer SuperImpositor::getTMScoreset2() const {
    return TMScoreset2;
}

Spacer SuperImpositor::getTMScoreset1() const {
    return TMScoreset1;
}

Spacer SuperImpositor::getGdtset4_2() const {
    return Gdtset4_2;
}

Spacer SuperImpositor::getGdtset4_1() const {
    return Gdtset4_1;
}

Spacer SuperImpositor::getGdtset3_2() const {
    return Gdtset3_2;
}

Spacer SuperImpositor::getGdtset3_1() const {
    return Gdtset3_1;
}

Spacer SuperImpositor::getGdtset2_2() const {
    return Gdtset2_2;
}

Spacer SuperImpositor::getGdtset2_1() const {
    return Gdtset2_1;
}

Spacer SuperImpositor::getGdtset1_2() const {
    return Gdtset1_2;
}

Spacer SuperImpositor::getGdtset1_1() const {
    return Gdtset1_1;
}

double SuperImpositor::getTMScoreValue() const {
    return TMScoreValue;
}

double SuperImpositor::getGdtValue() const {
    return gdtValue;
}

double SuperImpositor::getMaxsubValue() const {
    return maxsubValue;
}

double SuperImpositor::getRmsdValue() const {
    return rmsdValue;
}

std::vector<std::pair<int, int> > SuperImpositor::getTMScoreAlignment() const {
    return TMScoreAlignment;
}

std::vector<std::pair<int, int> > SuperImpositor::getGdtAlignment4() const {
    return gdtAlignment4;
}

std::vector<std::pair<int, int> > SuperImpositor::getGdtAlignment3() const {
    return gdtAlignment3;
}

std::vector<std::pair<int, int> > SuperImpositor::getGdtAlignment2() const {
    return gdtAlignment2;
}

std::vector<std::pair<int, int> > SuperImpositor::getGdtAlignment1() const {
    return gdtAlignment1;
}

std::vector<std::pair<int, int> > SuperImpositor::getMaxsubAlignment() const {
    return maxsubAlignment;
}