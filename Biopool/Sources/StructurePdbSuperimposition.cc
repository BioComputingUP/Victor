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

#include <Spacer.h>
#include <PdbSaver.h>
#include <Eigen/Geometry>
#include <StructurePdbSuperimposition.h>

using namespace Victor; using namespace Victor::Biopool;
using namespace Eigen;

//Rt hold the last rotation matrix and translation vector
std::pair<Eigen::Matrix3d, Eigen::Vector3d> Rt;

/**
 * Constructor
 * @param none
 */
StructurePdbSuperimposition::StructurePdbSuperimposition() {

}

/**
 * Destructor
 * @param none
 */
StructurePdbSuperimposition::~StructurePdbSuperimposition() {

}

/**
 * support method that return an Eigen matrix from the coordinates of Calpha,
 * in order to simplify the work with matrices 
 * @param Spacer*
 * @return Eigen::MatrixX3d
 */
Eigen::MatrixX3d StructurePdbSuperimposition::fromSpacerToEigenMatrix(Spacer *spacer) {
    int num = (int)spacer->sizeAmino();
    Eigen::MatrixX3d matrix;
    matrix.resize(num, 3);
    vgVector3<double> vectorCACoords;
    
    for (int i=0; i<num; i++) {
        vectorCACoords = spacer->getAmino(i)[CA].getCoords();
        matrix(i, 0) = vectorCACoords.x;
        matrix(i, 1) = vectorCACoords.y;
        matrix(i, 2) = vectorCACoords.z;
    }
    
    return matrix;
}

/**
 * 
 * @param Eigen::MatrixX3d&
 * @param Eigen::MatrixX3d&
 * @return std::pair<Eigen::Matrix3d,Eigen::Vector3d>
 */
std::pair<Eigen::Matrix3d,Eigen::Vector3d> StructurePdbSuperimposition::KabschAlgorithm(Eigen::MatrixX3d& A, Eigen::MatrixX3d& B) {
    Eigen::MatrixX3d P = A;
    Eigen::MatrixX3d Q = B;
    if ((P.rows() == Q.rows()) && (P.cols() == Q.cols()) && (P.rows() != 0)) {
        int size = P.rows();

        // STEP 1 - Translation
        Eigen::Vector3d centroidP = P.colwise().mean();
        Eigen::Vector3d centroidQ = Q.colwise().mean();

        for (int i = 0; i < size; i++) {
            P.row(i) -= centroidP;
            Q.row(i) -= centroidQ;
        }

        // STEP 2 - Covariance matrix
        Eigen::MatrixXd cov = P.transpose() * Q;

        // STEP 3 - Optimal rotation matrix
        //compute SVD
        JacobiSVD<Eigen::MatrixXd> svd(cov, ComputeThinU | ComputeThinV);

        //calculate determinant for the rotation sign  
        double det = (svd.matrixV() * svd.matrixU().transpose()).determinant();
        Eigen::Vector3d e(1, 1, (det < 0)? -1 : 1);

        //compute optimal rotation matrix    
        Eigen::Matrix3d U = (svd.matrixV() * e.asDiagonal() * svd.matrixU().transpose());

        //calculate optimal translation vector
        Eigen::Vector3d t = centroidQ - (U * centroidP);
        
        //return the best rotation matrix and translation vector
        return std::make_pair(U,t);
    }
    else {
        cout << "KabschAlgorithm(): Error matrix size" << endl;
    }
}

/**
 * calculate the RMSD score after optimal rototranslation returned by Kabsch method
 * @param Eigen::MatrixX3d&
 * @param Eigen::MatrixX3d&
 * @param Spacer&
 * @param Spacer&
 * @return double
 */
double StructurePdbSuperimposition::RMSD(Eigen::MatrixX3d& A, Eigen::MatrixX3d& B, Spacer& sp1, Spacer& sp2) {
    Eigen::MatrixX3d newA;
    
    Rt = KabschAlgorithm(A, B);
    
    applyMatrixRototranslation(newA, A);
    
    savePdbFile(sp1, sp2, "SuperimpositionRMSD.pdb");
    
    double rmsdScore=0.0;
    if((newA.rows() == B.rows()) && (newA.cols() == B.cols())) {
        for (int i=0; i<newA.rows(); i++)
            for (int j=0; j<newA.cols(); j++)
                rmsdScore += pow(newA(i,j) - B(i,j),2);
        rmsdScore /= newA.rows();
        rmsdScore = sqrt(rmsdScore);
    }
    else
        cout << "RMSD() error: the number of points are different" << endl;
    
    return rmsdScore;
}

/**
 * calculate and return the GDT_TS score with four distance using the iterative algorithm MaxSub
 * @param Eigen::MatrixX3d&
 * @param Eigen::MatrixX3d&
 * @param Spacer&
 * @param Spacer&
 * @param int (seed length for MaxSub)
 * @param double[] (scores of four distance)
 * @return double
 */
double StructurePdbSuperimposition::GDTTS(Eigen::MatrixX3d& A, Eigen::MatrixX3d& B, Spacer& sp1, Spacer& sp2, int L, double p[]) {
    p[0] = MaxSubScore(A, B, sp1, sp2, L, 1, false);
    savePdbFile(sp1, sp2, "SuperimpositionGDTTS1A.pdb");
    
    p[1] = MaxSubScore(A, B, sp1, sp2, L, 2, false);
    savePdbFile(sp1, sp2, "SuperimpositionGDTTS2A.pdb");
    
    p[2] = MaxSubScore(A, B, sp1, sp2, L, 4, false);
    savePdbFile(sp1, sp2, "SuperimpositionGDTTS4A.pdb");
    
    p[3] = MaxSubScore(A, B, sp1, sp2, L, 8, false);
    savePdbFile(sp1, sp2, "SuperimpositionGDTTS8A.pdb");
    
    return (p[0] + p[1] + p[2] + p[3])/4.0;
}

/**
 * calculate and return the TM score using the iterative algorithm MaxSub
 * @param Eigen::MatrixX3d&
 * @param Eigen::MatrixX3d&
 * @param Spacer&
 * @param Spacer&
 * @param int (seed length for MaxSub)
 * @return double
 */
double StructurePdbSuperimposition::TMScore(Eigen::MatrixX3d& A, Eigen::MatrixX3d& B, Spacer& sp1, Spacer& sp2, int L) {
    std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d> M_max;    
    int Lt, Ln = A.rows();
    double sum = 0, d0 = 1.24 * cbrt(Ln-15) - 1.8;
    if (d0 <= 0)
        d0 = 0.17; //approximated costant value
    
    M_max = MaxSub(A, B, sp1, sp2, L, d0);
    
    savePdbFile(sp1, sp2, "SuperimpositionTM.pdb");
    
    Lt = M_max.first.rows();
    
    for(int i=0; i<Lt; i++) {
	sum += ( 1.0 / (1.0 + pow((M_max.first.row(i)-M_max.second.row(i)).norm()/d0, 2)) );
    }
    return sum/(double)Ln;
}

/**
 * calculate and return the MaxSub score using iterative algorithm MaxSub
 * @param Eigen::MatrixX3d&
 * @param Eigen::MatrixX3d&
 * @param Spacer&
 * @param Spacer&
 * @param int (seed length for MaxSub)
 * @param double (threshold for MaxSub)
 * @param bool
 * @return double
 */
double StructurePdbSuperimposition::MaxSubScore(Eigen::MatrixX3d& A, Eigen::MatrixX3d& B, Spacer& sp1, Spacer& sp2, int L, double d, bool b) {
    std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d> M_max;
    double sizeA = A.rows();
    
    M_max = MaxSub(A, B, sp1, sp2, L, d);
    
    // If b is true then save pdb file for maxsub, otherwise for GDTTS
    if (b)
        savePdbFile(sp1, sp2, "SuperimpositionMaxSub.pdb");

    double sum = 0;
    for (int i=0; i<M_max.first.rows(); i++) {
        sum += 1.0 / (1.0 + pow((M_max.first.row(i)-M_max.second.row(i)).norm()/d,2));
    }
    
    //return the normalized score
    return sum/sizeA;    
}

/**
 * iterative algorithm that return the largest subset of atoms lying below the input distance
 * @param Eigen::MatrixX3d&
 * @param Eigen::MatrixX3d&
 * @param Spacer&
 * @param Spacer&
 * @param int (seed length)
 * @param double (threshold)
 * @return std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d>
 */
std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d> StructurePdbSuperimposition::MaxSub(Eigen::MatrixX3d& A, Eigen::MatrixX3d& B, Spacer& sp1, Spacer& sp2, int L, double d) {
    std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d> M, M_max;
    
    if((A.rows() == B.rows()) && (A.cols() == B.cols())) {
        int n = B.rows();
        int s_max = 0;
        for (int i=0; i<n-L; i++) {
            M.first = A.block(i,0,L,A.cols());
            M.second = B.block(i,0,L,B.cols());
                
            MaxSubExtend(M, A, B, d);
                
            if (M.first.rows() > s_max) {
                s_max = M.first.rows();
                M_max = M;                    
            }
        }
    }
    else
        cout << "MaxSub() error: the number of points are different" << endl;
    
    if(M_max.first.rows() == M_max.second.rows() && M_max.first.rows() > 0) {
        //update Rt with the best rototranslation, after MaxSub extend
        Rt = KabschAlgorithm(M_max.first, M_max.second);    
        applyMatrixRototranslation(M_max.first, M_max.first);
    }
    return M_max;
}

/**
 * support method for 'MaxSub' that try to extend the input set
 * @param std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d>&
 * @param Eigen::MatrixX3d&
 * @param Eigen::MatrixX3d&
 * @param double (threshold)
 * @return void
 */
void StructurePdbSuperimposition::MaxSubExtend(std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d>& M, Eigen::MatrixX3d& A, Eigen::MatrixX3d& B, double d) {
    std::pair<Eigen::MatrixX3d, Eigen::MatrixX3d> N;
    int k = 4, position;
    double threshold;
    
    for(int j=1; j<=k; j++) {
		if(M.first.rows() == M.second.rows() && M.first.rows() > 0)
        	Rt = KabschAlgorithm(M.first, M.second);
        
        Eigen::MatrixX3d Arototraslated;
        applyMatrixRototranslation(Arototraslated, A);
        
        N.first.resize(A.rows(),3);
        N.second.resize(A.rows(),3);
        
        threshold = (double(j)*d)/(double)k;
        position = 0;
        
        for(int i=0; i<A.rows(); i++) {
            if((Arototraslated.row(i) - B.row(i)).norm() <= threshold) {
                N.first.row(position) = A.row(i);
                N.second.row(position) = B.row(i);
                position++;
            }
        }
            
        N.first.conservativeResize(position,3);
        N.second.conservativeResize(position,3);
        
        M = N;
        
        //If no pair has been added, don't try to extend this seed
        if (position == 0)
            break;
    }
    
    if (position > 0) {
	Eigen::MatrixX3d M1rototraslated;
	applyMatrixRototranslation(M1rototraslated, M.first);
		        
	for(int i=0; i<M.first.rows(); i++) {
            if((M1rototraslated.row(i) - M.second.row(i)).norm() > d) {
		removeRow(M.first, i);
		removeRow(M.second, i);
		i--;
            }
	}
    }
}

/**
 * apply the global rototranslation Rt to the first input Eigen matrix
 * @param Eigen::MatrixX3d&
 * @param Eigen::MatrixX3d&
 * @return void
 */
void StructurePdbSuperimposition::applyMatrixRototranslation(Eigen::MatrixX3d& matrixToModify, Eigen::MatrixX3d& M) {
    matrixToModify = (Rt.first * M.transpose()).transpose();
    for(int i=0; i<M.rows(); i++) {
        matrixToModify.row(i) += Rt.second;
    }
}

/**
 * delete a row from the Eigen matrix in input
 * @param Eigen::MatrixX3d&
 * @param int
 * @return void
 */
void StructurePdbSuperimposition::removeRow(Eigen::MatrixX3d& M, int rowToRemove) {
    int numRows = M.rows()-1;
    int numCols = M.cols();

    if(rowToRemove < numRows)
        M.block(rowToRemove,0,numRows-rowToRemove,numCols) = M.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    M.conservativeResize(numRows,numCols);
}

/**
 * save the pdb file with the first spacer modify and the second spacer
 * @param Spacer&
 * @param Spacer&
 * @param string
 */
void StructurePdbSuperimposition::savePdbFile(Spacer& sp1, Spacer& sp2, string name) {
    ofstream outFile(name.c_str());
    if (!outFile)
	ERROR("File not found.", exception);
    Spacer sp = updateSpacer(sp1);
    PdbSaver ps = PdbSaver(outFile);
    ps.saveSpacer(sp);
    ps.saveSpacer(sp2);
    outFile.close();
}

/**
 * update the coordinates of the spacer in input with the last rototranslation
 * compute
 * @param Spacer&
 * @return Spacer
 */
Spacer StructurePdbSuperimposition::updateSpacer(Spacer& sp1) {    
    int numAmino = sp1.sizeAmino();
    Spacer updatedSp = Spacer(sp1);
    vector<Atom> aaAtoms, scAtoms;
    
    for(int i=0; i<numAmino; i++) {
        aaAtoms = updatedSp.getAmino(i).giveAtoms();
            
        for(unsigned int j=0; j<aaAtoms.size(); j++) {
            vgVector3<double> vectorVgCoords = aaAtoms[j].getCoords();
            Eigen::Vector3d vectorEigenCoords(vectorVgCoords.x, vectorVgCoords.y, vectorVgCoords.z);
            applyVectorRototranslation(vectorEigenCoords);
            updatedSp.getAmino(i).getAtom(j).setCoords(vectorEigenCoords(0),vectorEigenCoords(1),vectorEigenCoords(2));
        }
            
        scAtoms = updatedSp.getAmino(i).getSideChain().giveAtoms();
            
        for(unsigned k=0; k<scAtoms.size(); k++) {
            vgVector3<double> vectorVgCoords = scAtoms[k].getCoords();
            Eigen::Vector3d vectorEigenCoords(vectorVgCoords.x, vectorVgCoords.y, vectorVgCoords.z);
            applyVectorRototranslation(vectorEigenCoords);
            updatedSp.getAmino(i).getSideChain().getAtom(k).setCoords(vectorEigenCoords(0),vectorEigenCoords(1),vectorEigenCoords(2));
        }
    }
    
    return updatedSp;
}

/**
 * support method for 'updateSpacer', to apply the rototranslation on the Eigen input vector
 * @param Eigen::Vector3d&
 * @return void
 */
void StructurePdbSuperimposition::applyVectorRototranslation(Eigen::Vector3d& vectorToModify) {
    vectorToModify = Rt.first * vectorToModify;
    vectorToModify += Rt.second;
}
