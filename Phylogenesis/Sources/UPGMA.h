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

#ifndef __UPGMA_H_
#define __UPGMA_H_

#include <Cluster.h>
#include <HierarchicalClusterBuilder.h>
#include <DistanceMatrix.h>
#include <ClusterDistanceMatrix.h>
#include <vector>
#include "Sequence.h"

namespace Victor { namespace Phylogenesis{

	/** @brief   Hierarchical clustering by UPGMA (Unweighted Pair Group Method with Arithmetic Mean).
	 *  @author  Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class UPGMA : public HierarchicalClusterBuilder{
		public:
			/// Constructor.
			UPGMA(std::vector<Sequence*> species, DistanceMatrix* distanceSpecies);

			/// Perform the hierarchical clustering with UPGMA algorithm.
			Cluster* performClustering();

		private:

			/// Update the ClusterDistanceMatrix with the new Cluster
			void updateClusterDistanceMatrix(Cluster* newCluster);

			// ATTRIBUTES:

			/// Initial DistanceMatrix for given sequences
			DistanceMatrix* distanceSpecies;

			/// Clusters still not clustered
			std::vector<Cluster*> clusters;

			/// ClusterDistanceMatrix
			ClusterDistanceMatrix distanceCluster;
	};

}}

#endif
