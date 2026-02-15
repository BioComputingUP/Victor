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

#ifndef __ClusterMatrix_H__
#define __ClusterMatrix_H__

#include <map>

#include <Cluster.h>
#include <DistanceMatrix.h>

namespace Victor { namespace Phylogenesis{

	/** @brief    Implements a triangular matrix usefull to store distances between Cluster without wasting memory.
	 * 			  Access and update operations are performed in constant time, erase in linear time.
	 *  @author   Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class ClusterDistanceMatrix{
	public:
		/// Default constructor.
		ClusterDistanceMatrix();

		/// Access element
		double at(Cluster* c1, Cluster* c2);

		/// Erase all (c1,*), (c2,*), (*,c1), (*,c2) elements in the matrix.
		void erase(Cluster* c1, Cluster* c2);

		/// Add a new value to the matrix at (c1,c2)
		void insert(Cluster* c1, Cluster* c2, double distance);

		/// Clear the matrix
		void clear();

	private:

		// ATTRIBUTES

		std::map<Cluster*, std::map<Cluster*,double> > distanceClusters;
	};

}}

#endif
