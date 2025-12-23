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

#ifndef _KABSCH_H_
#define _KABSCH_H_

//Includes
#include <Eigen/Geometry> 
#include <RotoTranslation.h>

// Global constants, typedefs, etc. (to avoid):

namespace Victor 
{
    namespace Biopool 
    {
    /**@brief Implements the optimal roto-translation method required for superimposition
    * 
    * */
		class KabschAlgorithm : public RotoTranslation
		{

			
			private:

				//ATTRIBUTES:

				/** 
	     		* Eigen's Transform object holding the rotation matrix and the translation vector
	     		* Every time the method is invoked, the resulting roto-translation matrix is held here
	     		*/
				Eigen::Affine3d rotoTranslationMatrix;


				//METHODS: 

				bool kabsch(Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>);

			public:

				// CONSTRUCTORS/DESTRUCTOR:
				
				KabschAlgorithm();
				~KabschAlgorithm();


				//METHODS:


				//TEST FUNCTIONS: 
				
				Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd);
				void TestFind3DAffineTransform();
				bool kabschTest(Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>, double);
		

				//GET METHODS:

				Eigen::Affine3d getRotoTranslation();
				Eigen::Vector3d getTranslation();
				Eigen::Matrix3d getRotation();


				bool calculateRotoTranslation(Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>);
				Eigen::Matrix<double,Eigen::Dynamic,3> applyRotoTranslation(Eigen::Matrix<double,Eigen::Dynamic,3>);
				void resetRotoTranslation();	
		};	
	}
}

#endif //_KABSCH_H_