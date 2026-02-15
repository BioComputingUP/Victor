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

#include "Cluster.h"
#include "SingleCluster.h"
#include <limits>
#include <iostream>

namespace Victor { namespace Phylogenesis{

vector<Cluster*> Cluster::leafs;

void Cluster::initUnrootedTree(Cluster* cluster, Cluster* parent, double dist, int nseqs){
	cluster->dirs.resize(nseqs, -1);
	if(parent != NULL){
		cluster->branches.push_back(parent);
		cluster->lengths.push_back(dist);
	}
	if(cluster->isSingle()){
		leafs.push_back(cluster);
		return;
	}
	for(unsigned int i=0;i<cluster->branches.size();i++){
		if(cluster->branches[i] != parent){
			initUnrootedTree(cluster->branches[i], cluster, cluster->lengths[i], nseqs);
		}
	}
}

double Cluster::maxFarthestLeaf(Cluster* cluster, Cluster* parent, double dist, int id){

	if(cluster->isSingle() && parent != NULL){
		return dist;
	}

	double max = std::numeric_limits<double>::min();
	int dir = -1;

	for(unsigned int i=0;i<cluster->branches.size();i++){
		if(cluster->branches[i] != parent){
			double d = maxFarthestLeaf(cluster->branches[i], cluster, dist + cluster->lengths[i], id);
			if(d > max){
				max = d;
				dir = i;
			}
		}
	}

	cluster->dirs[id] = dir;
	return max;
}


void Cluster::removeParent(Cluster* cluster, Cluster* parent){
	if(parent != NULL){
		int parentId = -1;
		for(unsigned int i=0;i<cluster->branches.size();i++){
			if(cluster->branches[i] == parent){
				parentId = i;
			}
		}
		cluster->branches.erase(cluster->branches.begin()+parentId);
		cluster->lengths.erase(cluster->lengths.begin()+parentId);
		cluster->dirs.clear();
	}
}

void Cluster::removeParentRec(Cluster* cluster, Cluster* parent){
	removeParent(cluster, parent);

	if(cluster->isSingle()){
		return;
	}

	for(unsigned int i=0;i<cluster->branches.size();i++){
		removeParentRec(cluster->branches[i], cluster);
	}
}

Cluster* Cluster::addHalfwayRoot(Cluster* cluster, Cluster* parent, double dist, int id, double avg, double lastDist){

	if(dist < avg){
		Cluster* next =  cluster->branches[cluster->dirs[id]];
		double nextDist = cluster->lengths[cluster->dirs[id]];
		return addHalfwayRoot(next, cluster, dist + nextDist, id, avg, nextDist);
	}

	Cluster* root = new Cluster();
	root->branches.push_back(parent);
	root->branches.push_back(cluster);
	root->lengths.push_back(avg - (dist - lastDist));
	root->lengths.push_back(dist - avg);

	//Remove reverse edges
	removeParent(cluster, parent);
	removeParent(parent, cluster);
	removeParentRec(cluster, NULL);
	removeParentRec(parent, NULL);

	return root;
}

Cluster* Cluster::midpointRooting(Cluster* tree, int leafsnum){

	// Init the current cluster with parameters useful for making the rooted tree
	// like reverse edges and extract the leafs

	leafs.clear();
	initUnrootedTree(tree, NULL, 0.0, leafsnum);

	// For each leafs find the other most farthest leaf and keeps the max

	int maxId = -1;
	double maxDist = std::numeric_limits<double>::min();

	for(unsigned int i=0;i<leafs.size();i++){
		SingleCluster* leaf = static_cast<SingleCluster*>(leafs[i]);
		double dist = maxFarthestLeaf(leafs[i], NULL, 0.0, leaf->sequence->id);
		if(dist > maxDist){
			maxDist = dist;
			maxId = i;
		}
		//cout<<"Leaf id "<<leaf->sequence->id<<" with farthest leaf at "<<dist<<endl;
	}

	// Put the root halfway between the two leaves
	SingleCluster* l = static_cast<SingleCluster*>(leafs[maxId]);
	Cluster* rootedTree = addHalfwayRoot(l, NULL, 0.0, l->sequence->id, maxDist/2.0, 0.0);
	return rootedTree;
}

}}
