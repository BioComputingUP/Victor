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

#ifndef __NeighborJoining_H__
#define __NeighborJoining_H__

#include <Cluster.h>
#include <ClusterDistanceMatrix.h>
#include <HierarchicalClusterBuilder.h>
#include <DistanceMatrix.h>
#include "Sequence.h"

namespace Victor { namespace Phylogenesis{

	/** @brief   Hierarchical clustering by NeighborJoining.
	 *  @author  Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class NeighborJoining : public HierarchicalClusterBuilder
	{
		public:
			/// Constructor.
			NeighborJoining(std::vector<Sequence*> species, DistanceMatrix* distanceSpecies);

			/// Perform the hierarchical clustering with NJ algorithm.
			Cluster* performClustering();

		private:
			// ATTRIBUTES:

			/// Initial DistanceMatrix for given sequences
			DistanceMatrix* distanceSpecies;

			/// NJ MinimalDistanceCluster
			ClusterDistanceMatrix minimalDistanceCluster;

			/// ClusterDistanceMatrix
			ClusterDistanceMatrix distanceCluster;

			/// Clusters still not clustered
			std::vector<Cluster*> clusters;

			/// Mean distance between other clusters
			std::vector<double> meanDistances;

			void updateMinimalClusterDistance();
			void updateDistanceCluster(Cluster* newCluster);
			double meanDistance(Cluster* c);
	};

}}

#endif
