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

#include <iostream>
#include <utility>
#include <Eigen/Geometry> 
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <Protein.h>
#include <KabschAlgorithm.h>
#include <StructuralSuperimposition.h>



using namespace std;
using namespace Victor; 
using namespace Victor::Biopool;

/**
 *  Constructor
 *@param void
 */
StructuralSuperimposition::StructuralSuperimposition()
{}

/**
 *  DESTRUCTOR
 *@param void
 */
StructuralSuperimposition::~StructuralSuperimposition()
{}


/**
 * MaxSub heuristic iterative algorithm for finding the largest subset M such that, for each (ai,bi) in M,
 * the euclidean distance between ai and T(bi) is below the distance threshold d
 * @param Eigen::MatrixX3d (double nx3 matrix)  
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @param int
 * @param double
 * @return std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d>
 */
pair<Eigen::MatrixX3d,Eigen::MatrixX3d> StructuralSuperimposition::maxSub(Eigen::Matrix<double,Eigen::Dynamic,3> A, Eigen::Matrix<double,Eigen::Dynamic,3> B, int L, double d)
{
	int sMax=0; //size of the largest subset found so far

	pair<Eigen::MatrixX3d,Eigen::MatrixX3d> M, Mmax;

	if(A.rows() != B.rows())
	{
		cout<<"maxSub(): the two protein structures differ in the number of elements"<<endl;
	}
	else
	{
		if(A.rows()==0 || B.rows()==0)
		{
			cout<<"maxSub(): one or more data structures are empty"<<endl;
		}
		else
		{
			Eigen::Affine3d currentBestRotoTranslation;
			for(int i=0; i<A.rows()-L+1 ;i++)
			{
				
				M.first.resize(L,3);
				M.second.resize(L,3);
				for(int j=0;j<L;j++)
				{
					M.first.row(j)=A.row(j+i);
					M.second.row(j)=B.row(j+i);

				}
			
				currentBestRotoTranslation=extend(M,A,B,d);
				
				if(M.first.rows() > sMax)
				{
					
					sMax=M.first.rows();
					optimalRotoTranslation=currentBestRotoTranslation;
					
					Mmax.first=M.first;
					
					Mmax.second=M.second;
					
				}
				

			}
		}
	}
	return Mmax;
}

/**
 * Step of the iterative algorithm which tries to extend the initial seed set to include additional pairs
 * @param std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d> reference
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @param double
 * @return Eigen::Affine3d
 */
Eigen::Affine3d StructuralSuperimposition::extend(pair<Eigen::MatrixX3d,Eigen::MatrixX3d> &M, Eigen::Matrix<double,Eigen::Dynamic,3> A, Eigen::Matrix<double,Eigen::Dynamic,3> B, double d)
{
	double threshold;
	KabschAlgorithm kb= KabschAlgorithm();
	Eigen::Matrix<double,Eigen::Dynamic,3> transformedA;
	pair<Eigen::MatrixX3d,Eigen::MatrixX3d> newM,N;
	
	int k=4;
	
	for (int j=1; j<=k; j++)
	{
		//init N
		N.first.resize(A.rows(),3);
		N.second.resize(A.rows(),3);
		
		if(!kb.calculateRotoTranslation(M.first, M.second))
		{
			cout<<"extend(): an error occurred while calculating the roto-translation"<<endl;
		}
		else
		{
			if(!kb.kabschTest(M.first, M.second, 0.001))
				exit(-1);
			
			transformedA= kb.applyRotoTranslation(A);
			
			threshold=(j*d)/(double)k;
			
			int it=0;
			
			for(int i=0; i<A.rows();i++)
			{	
				if( (transformedA.row(i)-B.row(i)).norm()  <= threshold)
				{
					N.first.row(it)=A.row(i);
					N.second.row(it)=B.row(i);
					it++;
				}
			}
				N.first.conservativeResize(it,3);
				N.second.conservativeResize(it,3);
	
				M.first=N.first;
				M.second=N.second;
		}
	}

	if(!kb.calculateRotoTranslation(M.first, M.second))
	{
		cout<<"extend(): an error occurred while calculating the roto-translation"<<endl;
	}
	else
	{
		transformedA= kb.applyRotoTranslation(M.first);


		newM.first.resize(M.first.rows(),3);
		newM.second.resize(M.first.rows(),3);
		int j=0;
		
		for(int i=0;i<M.first.rows();i++)
		{
			if( (transformedA.row(i)-M.second.row(i)).norm() <= d)
			{
				newM.first.row(j)=M.first.row(i);
				newM.second.row(j)=M.second.row(i);
				j++;
			}
		}


		newM.first.conservativeResize(j,3);
		newM.second.conservativeResize(j,3);

		M.first=newM.first;
		M.second=newM.second;
	}
	return kb.getRotoTranslation();
}

/**
 * This function applies the optimal roto-translation (found by the iterative search algorithm) on the input matrix
 * and loads the protein atoms coordinates in the output Spacer
 * @param Spacer reference 
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @return Spacer
 */
Spacer StructuralSuperimposition::updateSpacerCoords(Spacer &spA, Eigen::Matrix<double,Eigen::Dynamic,3> A)
{
	unsigned int AANumber= spA.sizeAmino();
	Spacer updatedSp= Spacer(spA);
	vector<Atom> aaAtoms, sideChainAtoms;
	vgVector3<double> spatialAtomCoords;
	Eigen::Vector3d eigenAtomCoords;
	if(A.rows() == 0)
	{
		cout<<"updateSpacerCoords(): input data structure is empty. Spacer hasn't been updated\n";
	} 
	else
	{
		if(AANumber != A.rows())
		{
			cout<<"updateSpacerCoords(): spacer and matrix sizes don't match. Spacer hasn't been updated"<<endl;
		}
		else
		{ 

			Eigen::MatrixX3d tA= applyOptRotoTranslation(A);

			for(unsigned int i=0;i<AANumber;i++)
			{
				//get all the atoms of the i-th AA, using the function Group::getAtoms()
				aaAtoms = updatedSp.getAmino(i).giveAtoms();

				for(unsigned j=0;j<aaAtoms.size();j++)
				{
					spatialAtomCoords = aaAtoms[j].getCoords();

					eigenAtomCoords(0)= spatialAtomCoords[0];
					eigenAtomCoords(1)= spatialAtomCoords[1];
					eigenAtomCoords(2)= spatialAtomCoords[2];

					eigenAtomCoords = optimalRotoTranslation.linear() * eigenAtomCoords;
					eigenAtomCoords += optimalRotoTranslation.translation();

					updatedSp.getAmino(i).getAtom(j).setCoords(eigenAtomCoords(0),eigenAtomCoords(1),eigenAtomCoords(2));
				}


				//for each amino, get its side chain and update its atoms coords in the same way
				sideChainAtoms = updatedSp.getAmino(i).getSideChain().giveAtoms();
				
				for(unsigned k=0;k<sideChainAtoms.size();k++)
				{
					spatialAtomCoords = sideChainAtoms[k].getCoords();
					eigenAtomCoords(0)= spatialAtomCoords[0];
					eigenAtomCoords(1)= spatialAtomCoords[1];
					eigenAtomCoords(2)= spatialAtomCoords[2];
					eigenAtomCoords = optimalRotoTranslation.linear() * eigenAtomCoords;
					eigenAtomCoords += optimalRotoTranslation.translation();

					updatedSp.getAmino(i).getSideChain().getAtom(k).setCoords(eigenAtomCoords(0),eigenAtomCoords(1),eigenAtomCoords(2));
				}
			}
		}
	}
	return updatedSp;
}


/**
 * This function print a PDB file with the resulting optimally roto-translated structures 
 * @param string
 * @param Spacer reference 
 * @param Spacer reference
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @return void
 */
void StructuralSuperimposition::printPDBFile(string fileName, Spacer &spA, Spacer &spB, Eigen::Matrix<double,Eigen::Dynamic,3> A)
{
	ofstream outputFile;
	outputFile.open(fileName.c_str(), ios::out);
	
	if(outputFile.is_open())
	{
		Spacer updatedSp=updateSpacerCoords(spA, A);
		
		PdbSaver pdbS= PdbSaver(outputFile);
		pdbS.saveSpacer(updatedSp);
		pdbS.saveSpacer(spB);
		
		outputFile.close();
	}
	else
		cout<<"printPDBFile(): can't open the files: "<<fileName<<" write operation aborted."<<endl;
}

/**
 * This function copies the spacer's CA atoms coordinates inside an Eigen matrix
 * @param Spacer reference 
 * @return Eigen::MatrixX3d (double nx3 matrix)
 */
Eigen::Matrix<double,Eigen::Dynamic,3> StructuralSuperimposition::convertSpacerToEigen(Spacer *sp)
{
	Eigen::MatrixXd coordMatrix;
	unsigned int AANumber= sp->sizeAmino();
	if(AANumber!=0)
	{
		vgVector3<double> spatialCAlphaCoords;
		//cout<<"numero AA nello spacer: "<<AANumber<<endl;
		coordMatrix.resize(AANumber,3);
		
		for(unsigned int i=0;i<AANumber;i++)
		{
			spatialCAlphaCoords=sp->getAmino(i)[CA].getCoords();
			coordMatrix(i,0)=spatialCAlphaCoords[0];
			coordMatrix(i,1)=spatialCAlphaCoords[1];
			coordMatrix(i,2)=spatialCAlphaCoords[2];
		}
	}
	else
		cout<<"convertSpacerToEigen(): warning, the input matrix is empty.\n";
	return coordMatrix;
}


/**
 * Calculation of RMSD score formula
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @return double
 */
double StructuralSuperimposition::rmsd(Eigen::Matrix<double,Eigen::Dynamic,3> P, Eigen::Matrix<double,Eigen::Dynamic,3> Q)
{
	double rmsd=0.0;
	if(P.rows()!=Q.rows())
	{
		cout<<"rmsd(): the two protein structures differ in the number of points."<<endl; 
		return -1.0;
	}
	else
	{
		if(P.rows()!=0 && Q.rows()!=0)
		{
			for(int i=0;i<P.rows();i++)
			{
				for(int j=0;j<P.cols();j++)
				{
					rmsd+=pow( (P(i,j)-Q(i,j)), 2);
				}
			}
			rmsd= rmsd/(double)P.rows();
			return sqrt(rmsd);
		}
		else
		{
			cout<<"rmsd(): the input matrices are empty."<<endl;
			return -1.0;
		}
	}
}

/**
 * Calculation of the RMSD-Score 
 * @param Eigen::MatrixX3d (double nx3 matrix) 
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @param int
 * @param double
 * @return double
 */
double StructuralSuperimposition::getRmsdScore(Eigen::Matrix<double,Eigen::Dynamic,3> A, Eigen::Matrix<double,Eigen::Dynamic,3> B, int L, double d)
{
	pair<Eigen::MatrixX3d,Eigen::MatrixX3d> M = maxSub(A,B,L,d);
	
	Eigen::MatrixX3d tA= applyOptRotoTranslation(A);
	
	//M.first= applyOptRotoTranslation(M.first);
	//return rmsd(M.first,M.second);
	return rmsd(tA,B);
}

/**
 * Calculation of the TM-Score 
 * @param Eigen::MatrixX3d (double nx3 matrix) 
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @param int
 * @return double
 */
double StructuralSuperimposition::getTMScore(Eigen::Matrix<double,Eigen::Dynamic,3> A, Eigen::Matrix<double,Eigen::Dynamic,3> B, int L)
{
	double d0= getTmScoreScaleParam(A.rows());
	pair<Eigen::MatrixX3d,Eigen::MatrixX3d> M = maxSub(A,B,L,d0);

	return tmScore(M, A.rows());
}

/**
 * Calculation of the TM-Score scale parameter, based on the number N of protein residues
 * @param long int
 * @return double
 */
double StructuralSuperimposition::getTmScoreScaleParam(long int Ln)
{
	double cubicRootExp=1.0/3.0;
	//Since in C++ versions prior to C++11 there isn't a cubic root function, a pow() funct. with exp=1/3 must be used
	//An internal error raises if the base is finite negative and the exponent is finite but not an integer value, causing a domain error.
	if(Ln-15 <= 0) 
	{
		cout<<"getTmScoreScaleParam(): warning, couldn't evaluate d0 appropriately, using approximated constant value of 0.17\n ";
		return 0.17;
	}
	else	
		return (1.24 * pow((Ln-15),cubicRootExp) - 1.8);
}

/**
 * Calculation of TM-Score formula
 * @param std::pair<Eigen::MatrixX3d>,Eigen::MatrixX3d>
 * @param long int
 * @return double
 */
double StructuralSuperimposition::tmScore(pair<Eigen::MatrixX3d,Eigen::MatrixX3d> M, long int Ln)
{
	long int Lt=M.first.rows();
	double di, partialSum=0;
	double d0;
	if(Ln<=0)
	{
		cout<<"tmScore(): can't evaluate a structure with 0 residues."<<endl;
		return -1;
	}
	else
	{
		d0=getTmScoreScaleParam(Ln);
		cout<<"TM-Score calculated with d0 threshold: "<<d0<<endl;
		M.first= applyOptRotoTranslation(M.first);
		for(int i=0;i<Lt;i++)
		{
			di=(M.first.row(i)-M.second.row(i)).norm();
			partialSum+= ( 1.0/ (1.0+pow(di/d0,2)) );
		}
		return partialSum/Ln;
	}

}

/**
 * Calculation of MaxSub score from the maximum subset obtained from the search algorithm
 * @param std::pair<Eigen::MatrixX3d,Eigen::MatrixX3d>
 * @param double
 * @param long int
 * @return double
 */
double StructuralSuperimposition::maxSubScore(pair<Eigen::MatrixX3d,Eigen::MatrixX3d> M, double d, long int q)
{
	double sum=0, di=0;
	if(q==0)
	{
		cout<<"maxSubScore(): can't evaluate a structure with 0 residues."<<endl;
		return -1;
	}
	else
	{
		for(int i=0;i<M.first.rows();i++)
		{
			di= (M.first.row(i)-M.second.row(i)).norm();
			sum= sum+ ( 1.0/(1.0 + pow(di/d,2)) );
		}
		return sum/(double)q;
	}
}

/**
 * Applies the current optimal roto-translation matrix to a copy of the input matrix
 * @param Eigen::MatrixX3d (double nx3 matrix) 
 * @return Eigen::MatrixX3d (double nx3 matrix)
 */
Eigen::Matrix<double,Eigen::Dynamic,3> StructuralSuperimposition::applyOptRotoTranslation(Eigen::Matrix<double,Eigen::Dynamic,3> mat)
{
	Eigen::Matrix<double,Eigen::Dynamic,3> result;
	Eigen::Matrix<double,1,3> translationVector= optimalRotoTranslation.translation();

	result=(optimalRotoTranslation.linear() *(mat.transpose())).transpose();

	for(int i=0;i<mat.rows();i++)
	{
		result.row(i)+= translationVector;
	}
	return result;
}

/**
 * Calculation of the MaxSub score 
 * @param Eigen::MatrixX3d (double nx3 matrix) 
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @param int
 * @param double
 * @return double
 */
double StructuralSuperimposition::getMaxSubScore(Eigen::Matrix<double,Eigen::Dynamic,3> A, Eigen::Matrix<double,Eigen::Dynamic,3> B, int L, double d)
{
	resetOptRTMatrix();
	pair<Eigen::MatrixX3d,Eigen::MatrixX3d> M = maxSub(A,B,L,d);
	M.first= applyOptRotoTranslation(M.first);
	return maxSubScore(M,d,A.rows());
}


/**
 * Calculation of GTD-TS score 
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @param Eigen::MatrixX3d (double nx3 matrix)
 * @return double
 */
double StructuralSuperimposition::getGDTTSScore(Eigen::Matrix<double,Eigen::Dynamic,3> A, Eigen::Matrix<double,Eigen::Dynamic,3> B, int L)
{
	double gdt1A, gdt2A, gdt4A, gdt8A;

	gdt1A= getMaxSubScore(A,B,L,1.0);
	
	gdt2A= getMaxSubScore(A,B,L,2.0);
	
	gdt4A= getMaxSubScore(A,B,L,4.0);
	
	gdt8A= getMaxSubScore(A,B,L,8.0);

	cout<<"gdt-ts score (d<1A): "<<gdt1A<<endl;
	cout<<"gdt-ts score (d<2A): "<<gdt2A<<endl;
	cout<<"gdt-ts score (d<4A): "<<gdt4A<<endl;
	cout<<"gdt-ts score (d<8A): "<<gdt8A<<endl;

	return (gdt1A + gdt2A + gdt4A + gdt8A)/4.0;
}

/**
 * Reset the current optimal roto-translation matrix
 * @param void 
 * @return void
 */
void StructuralSuperimposition::resetOptRTMatrix()
{
	optimalRotoTranslation.linear() = Eigen::Matrix3d::Identity(3, 3);
	optimalRotoTranslation.translation() = Eigen::Vector3d::Zero();
}
