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

#include <NeighborJoining.h>
#include <SingleCluster.h>
#include <iostream>
#include <limits>

namespace Victor { namespace Phylogenesis{

	/**
	 *
	 * @param species initial sequences
	 * @param distanceSpecies distance between sequences
	 */
	NeighborJoining::NeighborJoining(std::vector<Sequence*> species, DistanceMatrix* distanceSpecies) : distanceSpecies(distanceSpecies)
	{
		for (unsigned int i = 0; i < species.size(); i++)
		{
			Cluster* cluster = new SingleCluster(species[i]);
			(cluster->indices).push_back(species[i]->id);
			clusters.push_back(cluster);
		}
	}

	double NeighborJoining::meanDistance(Cluster* c)
	{
		double mean = 0.0;
		for (unsigned int i = 0; i < clusters.size(); i++)
		{
			Cluster* current = clusters[i];
			if(current != c)
			{
				mean+= distanceCluster.at(current,c);
			}
		}
		mean = mean / double((clusters.size() - 2));

		return mean;
	}

	void NeighborJoining::updateMinimalClusterDistance()
	{
		meanDistances.clear();

		for (unsigned int i = 0; i < clusters.size(); i++)
		{
			meanDistances.push_back(meanDistance(clusters[i]));
		}
	
		minimalDistanceCluster.clear();

		for (unsigned int i = 0; i < clusters.size() - 1; i++)
		{
			for (unsigned int j = i+1; j < clusters.size(); j++)
			{
				// M(i,j);
				double distance_ij = distanceCluster.at(clusters[i],clusters[j]);
				double r_i = meanDistances[i];
				double r_j = meanDistances[j];
				double distance = distance_ij - r_i - r_j;

				minimalDistanceCluster.insert(clusters[i],clusters[j], distance);
			}
		}
	}
	
	void NeighborJoining::updateDistanceCluster(Cluster* newCluster)
	{
		Cluster* a = newCluster->branches[0];
		Cluster* b = newCluster->branches[1];

		//Note: a and b was removed from clusters list
		for (unsigned int i = 0; i < clusters.size(); i++)
		{
			Cluster* c = clusters[i];
			double distance = (distanceCluster.at(a,c) + distanceCluster.at(b,c) - distanceCluster.at(a,b)) / 2.0;
			distanceCluster.insert(newCluster, c, distance);
		}
	}

	Cluster* NeighborJoining::performClustering()
	{
		std::cout<<"\nPerforming NJ clustering...\n";

		//Init distanceCluster
		for (unsigned int i = 0; i < clusters.size() - 1; i++)
		{
			for (unsigned int j = i+1; j < clusters.size(); j++)
			{
				distanceCluster.insert(clusters[i],clusters[j], distanceSpecies->at(i,j));
			}
		}

		if(clusters.size() == 1){
			std::cout<<"done\n";
			return clusters[0];
		}

		while(clusters.size() > 2)
		{
			updateMinimalClusterDistance();

			// Find the two clusters to 'smaller' distance according to NJ (minimalDistanceCluster)
			double minDistance = std::numeric_limits<double>::max();
			Cluster* firstToMerge = NULL;
			Cluster* secondToMerge = NULL;
			int firstIndex = -1;
			int secondIndex = -2;

			for (unsigned int i = 0; i < clusters.size() - 1; i++)
			{
				for (unsigned int j = i+1; j < clusters.size(); j++)
				{
					Cluster* first = clusters[i];
					Cluster* second = clusters[j];
					// NJ distance
					double pairDistance =  minimalDistanceCluster.at(first, second);

					if(pairDistance < minDistance){
						minDistance = pairDistance;
						firstToMerge = first;
						secondToMerge = second;
						firstIndex = i;
						secondIndex = j;
					}
				}
			}

			// Build the new cluster
			Cluster* merged = new Cluster();
			merged->branches.push_back(firstToMerge);
			merged->branches.push_back(secondToMerge);

			// Update species in the current cluster
			std::vector<int> ms = std::vector<int>(firstToMerge->indices);
			ms.insert(ms.end(),secondToMerge->indices.begin(),secondToMerge->indices.end());
			merged->indices = ms;

			//branchLength a,merged = d(A,B) / 2 + [r(A)-r(B)] / 2
			//branchLength b,merged = d(A,B) / 2 + [r(B)-r(A)] / 2
			double lenghtToFirst = (distanceCluster.at(firstToMerge,secondToMerge) / 2.0) + ((meanDistances[firstIndex] - meanDistances[secondIndex])/2.0);
			double lenghtToSecond = distanceCluster.at(firstToMerge,secondToMerge) - lenghtToFirst;
			merged->lengths.push_back(lenghtToFirst);
			merged->lengths.push_back(lenghtToSecond);

			// Remove merged clusters from the initial list
			std::vector<Cluster*> clustersUpdated;
			for (unsigned int i = 0; i < clusters.size(); i++)
			{
				if(clusters[i] == firstToMerge || clusters[i] == secondToMerge)
					continue;

				clustersUpdated.push_back(clusters[i]);
			}
			clusters = clustersUpdated;

			// Update distanceCluster
			updateDistanceCluster(merged);

			// Remove merged data distances
			distanceCluster.erase(firstToMerge, secondToMerge);

			// Add the new cluster
			clusters.push_back(merged);
		}

		// Merge last 2 clusters with an edge
		clusters[1]->branches.push_back(clusters[0]);
		clusters[1]->lengths.push_back(distanceCluster.at(clusters[0],clusters[1]));

		return new Cluster(*clusters[1]);
	}

}}
