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

#ifndef _STRUCTURAL_SUPERIMPOSITION_H_
#define _STRUCTURAL_SUPERIMPOSITION_H_

//Includes
#include <Eigen/Geometry> 
#include <PdbLoader.h>
#include <Protein.h>

// Global constants, typedefs, etc. (to avoid):


namespace Victor 
{
    namespace Biopool 
    {
    	/**@brief Implements the search algorithms and the corresponding score metrics
    	* 
    	* */
		class StructuralSuperimposition
		{

			private:
				
				//ATTRIBUTES

				/** 
     			* Eigen's Transform object, holding the optimal rotation matrix and translation vector
     			* resulting from applying the iterative search algorithm on the two protein structures
     			*/
				Eigen::Affine3d optimalRotoTranslation;


				//METHODS

				Eigen::Affine3d extend(pair<Eigen::MatrixX3d,Eigen::MatrixX3d>&, Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>, double);
				
				double getTmScoreScaleParam(long int);

			public:
				
				// CONSTRUCTORS/DESTRUCTOR:

				StructuralSuperimposition();
				~StructuralSuperimposition();


				//METHODS:

				Eigen::Matrix<double,Eigen::Dynamic,3> applyOptRotoTranslation(Eigen::Matrix<double,Eigen::Dynamic,3>);

				//INTERFACE METHODS BETWEEN Eigen and Victor FORMATS

				static Eigen::Matrix<double,Eigen::Dynamic,3> convertSpacerToEigen(Spacer *);
				Spacer updateSpacerCoords(Spacer &, Eigen::Matrix<double,Eigen::Dynamic,3>);
				

				void printPDBFile(string,Spacer &, Spacer &, Eigen::Matrix<double,Eigen::Dynamic,3>);
				

				//ITERATIVE SEARCH ALGORITHM 

				pair<Eigen::MatrixX3d,Eigen::MatrixX3d> maxSub(Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>, int, double);
				

				//SCORE METHODS: 

				double getMaxSubScore(Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>, int, double);
				double getTMScore(Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>, int);
				double getGDTTSScore(Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>, int);
				double getRmsdScore(Eigen::Matrix<double,Eigen::Dynamic,3>,Eigen::Matrix<double,Eigen::Dynamic,3>,int, double);
				//void tmScoreSearchEngine (Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>, int, double);
				
				double tmScore(pair<Eigen::MatrixX3d,Eigen::MatrixX3d>, long int);
				double maxSubScore(pair<Eigen::MatrixX3d,Eigen::MatrixX3d>, double, long int);
				double rmsd(Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>);

				void resetOptRTMatrix();
		};
	}
}


#endif //_STRUCTURAL_SUPERIMPOSITION_H_