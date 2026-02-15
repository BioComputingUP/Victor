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

#include <UPGMA.h>
#include <SingleCluster.h>
#include <iostream>
#include <limits>

namespace Victor { namespace Phylogenesis{

	/**
	 *
	 * @param species initial sequences
	 * @param distanceSpecies distance between sequences
	 */
	UPGMA::UPGMA(std::vector<Sequence*> species, DistanceMatrix* distanceSpecies) : distanceSpecies(distanceSpecies)
	{
		for (unsigned int i = 0; i < species.size(); i++)
		{
			Cluster* cluster = new SingleCluster(species[i]);
			(cluster->indices).push_back(species[i]->id);
			clusters.push_back(cluster);
		}
	}

	Cluster* UPGMA::performClustering()
	{
		//Init distanceCluster
		for (unsigned int i = 0; i < clusters.size() - 1; i++)
		{
			for (unsigned int j = i+1; j < clusters.size(); j++)
			{
				distanceCluster.insert(clusters[i],clusters[j], distanceSpecies->at(i,j));
			}
		}

		std::cout<<"\nPerforming UPGMA clustering\n";

		while(clusters.size() != 1)
		{
			// Find the two clusters at minimum distance
			double minDistance = std::numeric_limits<double>::max();
			Cluster* firstToMerge = NULL;
			Cluster* secondToMerge = NULL;

			for (unsigned int i = 0; i < clusters.size() - 1; i++)
			{
				for (unsigned int j = i+1; j < clusters.size(); j++)
				{
					Cluster* first = clusters[i];
					Cluster* second = clusters[j];
					double pairDistance =  distanceCluster.at(first, second);

					if(pairDistance < minDistance){
						minDistance = pairDistance;
						firstToMerge = first;
						secondToMerge = second;
					}
				}
			}

			// Build the new cluster
			Cluster* merged = new Cluster();
			merged->branches.push_back(firstToMerge);
			merged->branches.push_back(secondToMerge);

			// Current cluster distance
			double distance = minDistance / 2.0;
			merged->distance = distance;

			// Distance to clusters
			merged->lengths.push_back(distance - firstToMerge->distance);
			merged->lengths.push_back(distance - secondToMerge->distance);

			// Update species in the current cluster
			std::vector<int> ms = std::vector<int>(firstToMerge->indices);
			ms.insert(ms.end(),secondToMerge->indices.begin(),secondToMerge->indices.end());
			merged->indices = ms;

			// Remove merged data distances
			distanceCluster.erase(firstToMerge, secondToMerge);

			// Remove merged clusters from the initial list
			std::vector<Cluster*> clustersUpdated;
			for (unsigned int i = 0; i < clusters.size(); i++)
			{
				if(clusters[i] == firstToMerge || clusters[i] == secondToMerge)
					continue;

				clustersUpdated.push_back(clusters[i]);
			}
			clusters = clustersUpdated;

			// Add new distance from other clusters
			updateClusterDistanceMatrix(merged);

			// Add the new cluster
			//clusters.push_back(merged);
			//To keep original ClustalW order
			clusters.insert(clusters.begin(), merged);
		}

		Cluster* root = clusters.front();
		root->setRooted(true);
		return root;
	}

	void UPGMA::updateClusterDistanceMatrix(Cluster* newCluster)
	{
		for (unsigned int i = 0; i < clusters.size(); i++)
		{
			Cluster* c = clusters[i];

			// Calculates the distance based on the present species
			double distance = 0.0;
			std::vector<int> leftS(newCluster->indices);
			std::vector<int> rightS(c->indices);
			
			for (unsigned int ii = 0; ii < leftS.size(); ii++)
			{
				for (unsigned int jj = 0; jj < rightS.size(); jj++)
				{
					distance += distanceSpecies->at(leftS[ii],rightS[jj]);
				}
			}
			distance = distance / (leftS.size() * rightS.size());

			// Add it
			distanceCluster.insert(c, newCluster, distance);
		}
	}

}}
