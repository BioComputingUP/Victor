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
 * Author: Filippo Gamberoni
 *
 */

#include <math.h>
#include <iostream>
#include <Eigen/SVD>
#include <Eigen/Geometry> 
#include <KabschAlgorithm.h>


//using namespace Eigen;
using namespace std;
using namespace Victor;
using namespace Victor::Biopool;

/**
 *  Constructor
 *@param void
 */
KabschAlgorithm::KabschAlgorithm()
{}

/**
 *  DESTRUCTOR
 *@param none
 */
KabschAlgorithm::~KabschAlgorithm()
{}


/**
 * Get the stored roto-translation matrix as a Eigen::Transform object
 * @param void 
 * @return Eigen::Affine3d
 */
Eigen::Affine3d KabschAlgorithm::getRotoTranslation()
{
	return rotoTranslationMatrix;
}
	
/**
 * Get the stored translation vector 
 * @param void 
 * @return Eigen::Vector3d
 */		
Eigen::Vector3d KabschAlgorithm::getTranslation()
{
	return rotoTranslationMatrix.translation();
}

/**
 * Get the stored rotation matrix 
 * @param void 
 * @return Eigen::Matrix3d (double 3x3 matrix)
 */
Eigen::Matrix3d KabschAlgorithm::getRotation()
{
	return rotoTranslationMatrix.linear();
}

/**
 * Reset the current roto-translation matrix
 * @param void 
 * @return void
 */
void KabschAlgorithm::resetRotoTranslation()
{
	rotoTranslationMatrix.linear() = Eigen::Matrix3d::Identity(3, 3);
	rotoTranslationMatrix.translation() = Eigen::Vector3d::Zero();
}

/**
 * Applies the current optimal roto-translation matrix to a copy of the input matrix
 * @param Eigen::MatrixX3d (double nx3 matrix) 
 * @return Eigen::MatrixX3d (double nx3 matrix)
 */
Eigen::Matrix<double,Eigen::Dynamic,3> KabschAlgorithm::applyRotoTranslation(Eigen::Matrix<double,Eigen::Dynamic,3> mat)
{
	Eigen::Matrix<double,Eigen::Dynamic,3> result;
	Eigen::Matrix<double,1,3> translationVector= rotoTranslationMatrix.translation();

	result=(rotoTranslationMatrix.linear() *(mat.transpose())).transpose();

	for(int i=0;i<mat.rows();i++)
	{
		result.row(i)+= translationVector;
	}
	return result;
}

/**
 * Wrapper to the actual Kabsch method
 * @param Eigen::MatrixX3d (double nx3 matrix)  
 * @param Eigen::MatrixX3d (double nx3 matrix)  
 * @return bool
 */
bool KabschAlgorithm::calculateRotoTranslation(Eigen::Matrix<double,Eigen::Dynamic,3> P, Eigen::Matrix<double,Eigen::Dynamic,3> Q)
{
	return KabschAlgorithm::kabsch(P,Q);
}

/**
 * Alternative Kabsch method to assess the correctness of the one used 
 * (found here: https://en.wikipedia.org/wiki/Kabsch_algorithm#External_links)
 * @param Eigen::MatrixX3d (double nx3 matrix)  
 * @param Eigen::MatrixX3d (double nx3 matrix)  
 * @return Eigen::Affine3d
 */
Eigen::Affine3d KabschAlgorithm::Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out) 
{
	// Default output
	Eigen::Affine3d A;
	A.linear() = Eigen::Matrix3d::Identity(3, 3);
	A.translation() = Eigen::Vector3d::Zero();
	
	if (in.cols() != out.cols())
		cout<<"Find3DAffineTransform(): input data mis-match"<<endl;
	
	// Find the centroids then shift to the origin
	Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
	Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
	for (int col = 0; col < in.cols(); col++) 
	{
		in_ctr += in.col(col);
		out_ctr += out.col(col);
	}
	in_ctr /= in.cols();
	out_ctr /= out.cols();
	for (int col = 0; col < in.cols(); col++) 
	{
		in.col(col) -= in_ctr;
		out.col(col) -= out_ctr;
	}
	
	// SVD
	Eigen::MatrixXd Cov = in * out.transpose();
	
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

	// Find the rotation
	double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
	if (d > 0)
		d = 1.0;
	else
		d = -1.0;
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
	I(2, 2) = d;
	Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();
	
	// The final transform
	A.linear() =  R;
	
	A.translation() = (out_ctr - R*in_ctr);
	
	return A;
}
	
/**
 * A function to test Find3DAffineTransform()
 * @param void 
 * @return void
 */
void KabschAlgorithm::TestFind3DAffineTransform(){
	// Create datasets with known transform
	Eigen::Matrix3Xd in(3, 100), out(3, 100);
	Eigen::Quaternion<double> Q(1, 3, 5, 2);
	Q.normalize();
	Eigen::Matrix3d R = Q.toRotationMatrix();
	double scale = 2.0;
	for (int row = 0; row < in.rows(); row++) 
	{
		for (int col = 0; col < in.cols(); col++) 
		{
			in(row, col) = log(2*row + 10.0)/sqrt(1.0*col + 4.0) + sqrt(col*1.0)/(row + 1.0);
		}
	}
	Eigen::Vector3d S;
	S << -5, 6, -27;
	for (int col = 0; col < in.cols(); col++)
		out.col(col) = scale*R*in.col(col) + S;
	Eigen::Affine3d A = Find3DAffineTransform(in, out);
	// See if we got the transform we expected
	if ( (scale*R-A.linear()).cwiseAbs().maxCoeff() > 1e-13 || (S-A.translation()).cwiseAbs().maxCoeff() > 1e-13)
		throw "Could not determine the affine transform accurately enough";
}

/**
 * Tests the roto-translation matrices get by kabsch() and Find3DAffineTransform(), with custom tolerance
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @param double
 * @return bool
 */
bool KabschAlgorithm::kabschTest(Eigen::Matrix<double,Eigen::Dynamic,3> A, Eigen::Matrix<double,Eigen::Dynamic,3> B, double threshold)
{
	//push
	Eigen::Affine3d rt=getRotoTranslation();

	Eigen::Affine3d result1=Find3DAffineTransform(A.transpose(), B.transpose());
	kabsch(A,B);
	Eigen::Affine3d result2= getRotoTranslation();
	
	double fail=0.0;
	for(int i=0;i<result1.linear().rows();i++)
	{
		for(int j=0;j<3;j++)
		{
			if( fabs(result1.linear()(i,j)+fail - result2.linear()(i,j) ) > threshold )
			{
				cout<<"("<<i<<","<<j<<"): "<<fabs(result1.linear()(i,j) +fail - result2.linear()(i,j) )<<endl;
				return false;
			}
		}
	}

	for(int i=0;i<3;i++)
	{
		if( fabs( result1.translation()(i)-result2.translation()(i) ) > threshold )
			return false;
	}

	//pop
	rotoTranslationMatrix= rt;
	return true;
}

/**
 * Calculates the optimal roto-translation matrix between the two input matrices. 
 * Matrix B is superimposed onto A, and the roto-translation is kept in the object member
 * @param Eigen::MatrixX3d (double nx3 matrix) A 
 * @param Eigen::MatrixX3d (double nx3 matrix) B
 * @return bool
 */
bool KabschAlgorithm::kabsch(Eigen::Matrix<double,Eigen::Dynamic,3> A, Eigen::Matrix<double,Eigen::Dynamic,3> B)
{
	//centroids
	resetRotoTranslation();
	
	Eigen::Matrix<double,Eigen::Dynamic,3> P=A;
	Eigen::Matrix<double,Eigen::Dynamic,3> Q=B;

	bool success=false;

	if(P.rows() != Q.rows())
	{
		cout<<"kabsch(): the two protein structures differ in the number of points"<<endl;
	}
	else
	{
		if(P.rows()==0 || Q.rows()==0)
		{
			success=true; 
		}
		else
		{
			double det=0.0;
			
			Eigen::Vector3d centroidP= Eigen::Vector3d::Zero();
			Eigen::Vector3d centroidQ= Eigen::Vector3d::Zero();
			Eigen::Vector3d t;
			Eigen::MatrixXd H(P.rows(),P.rows());
			Eigen::MatrixXd S,U,V,R,Id;

			Eigen::MatrixXd Pt,Qt;
			
			//Accumulo distintamente le coordinate x,y,z di ciascun punto
			for(int i=0;i<P.rows();i++)
			{
				//Somma di vettori: sommo  il vettore centroid con la i-esima riga della matrice
				centroidP+=P.row(i);
				centroidQ+=Q.row(i);
			}

			//calcolo della media di ciascuna coordinata
			//divisione del vettore per uno scalare
			centroidP/= P.rows();
			centroidQ/= Q.rows();

			
			//coordinates traslation
			//Sottraggo dalle coordinate di ciascun punto le coordinate del rispettivo centroide 
			Pt=P;
			Qt=Q;
			for(int i=0;i<P.rows();i++)
			{
				P.row(i)-=centroidP;
				Q.row(i)-=centroidQ;
			}

			
			H=(P.transpose())*Q;
			
			//SVD
			Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
			S=svd.singularValues();
			U= svd.matrixU();
			V= svd.matrixV();
			
			det= (V*U.transpose()).determinant();
			
			if(det>0)
				det=1.0;
			else
				det=-1.0;

			Id= Eigen::Matrix3d::Identity(H.rows(),H.rows());
			Id(H.rows()-1,H.rows()-1)=det;
			
			R= V*Id*U.transpose();
			
			rotoTranslationMatrix.linear()=R;
			//Traslation  t= -R * centroidA + centroidB
			t= -R*centroidP + centroidQ;
			
			rotoTranslationMatrix.translation()=t;


	
			success=true;
		}
	}
	return success;
}
