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

#include "ClusterDistanceMatrix.h"

namespace Victor { namespace Phylogenesis{

	ClusterDistanceMatrix::ClusterDistanceMatrix(){}

	double ClusterDistanceMatrix::at(Cluster* c1, Cluster* c2)
	{
		if(distanceClusters.find(c1) != distanceClusters.end() && distanceClusters.at(c1).find(c2) != distanceClusters.at(c1).end()){
			return distanceClusters.at(c1).at(c2);
		}
		return distanceClusters.at(c2).at(c1);
	}

	void ClusterDistanceMatrix::erase(Cluster* c1, Cluster* c2)
	{
		std::map<Cluster*, std::map<Cluster*,double> > updatedDistanceClusters;
		for (std::map<Cluster*, std::map<Cluster*,double> >::iterator it=distanceClusters.begin(); it!=distanceClusters.end(); ++it)
		{
			if(it->first == c1 || it->first  == c2)
				continue;

			std::map<Cluster*,double> row = it->second;
			std::map<Cluster*,double> updatedRow;
			for(std::map<Cluster*,double>::iterator itr=row.begin(); itr!=row.end(); ++itr)
			{
				if(itr->first == c1 || itr->first  == c2)
					continue;

				updatedRow.insert(*itr);
			}
			updatedDistanceClusters.insert(std::pair<Cluster*, std::map<Cluster*,double> >(it->first, updatedRow));
		}
		distanceClusters = updatedDistanceClusters;
	}

	void  ClusterDistanceMatrix::insert(Cluster* c1, Cluster* c2, double distance)
	{
		if(distanceClusters.find(c1) == distanceClusters.end() && distanceClusters.find(c2) != distanceClusters.end()){
			distanceClusters[c2].insert(std::pair<Cluster*,double>(c1, distance));
		}
		else if(distanceClusters.find(c2) == distanceClusters.end() && distanceClusters.find(c1) != distanceClusters.end()){
			distanceClusters[c1].insert(std::pair<Cluster*,double>(c2, distance));
		}
		else{
			std::map<Cluster*,double> row;
			distanceClusters.insert(std::pair<Cluster*, std::map<Cluster*,double> >(c1, row));
			distanceClusters[c1].insert(std::pair<Cluster*,double>(c2, distance));
		}
	}

	void ClusterDistanceMatrix::clear(){
		distanceClusters.clear();
	}
}}
