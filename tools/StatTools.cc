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
/**
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "StatTools.h"

/**
 * @Description get Random Number
 * @param double low, double high 
* @return double
*/
double getRandomNumber( double low, double high )
{
	double interval=max<double>( low, high )-min<double>(low, high );
	double delta=min<double>( low, high );
	double randomNumber=((double)rand())/(((double)RAND_MAX)+((double)1));
	return randomNumber*interval+delta;
}

/**
 * @Description get Gaussian Random Number
 * @param double low, double high double seed
* @return double
*/
double getGaussianRandomNumber( double low, double high, double seed )
{
	double sum=0.0;
	for ( long i=0; i<static_cast<long>(seed); ++i )
	{
		double rand=getRandomNumber(0,1);
		sum += rand;
	}
	double interval=max<double>( low, high )-min<double>(low, high );
	return interval*(sum-seed/2)/(seed/2);
}



