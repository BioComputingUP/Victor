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

#ifndef _ROTO_TRANSLATION_H_
#define _ROTO_TRANSLATION_H_

//Includes
#include <Eigen/Geometry> 

// Global constants, typedefs, etc. (to avoid):

namespace Victor 
{
    namespace Biopool 
    {

		/**@brief Abstract class for the implementation of different roto-translation methods
	    * 
	    * */
		class RotoTranslation
		{
			public:
				
				// CONSTRUCTORS/DESTRUCTOR:

				RotoTranslation();
				~RotoTranslation();
				

				//VIRTUAL METHOD:

				virtual bool calculateRotoTranslation(Eigen::Matrix<double,Eigen::Dynamic,3>, Eigen::Matrix<double,Eigen::Dynamic,3>) =0;
		};	
	}
}

#endif //_ROTO_TRANSLATION_H_