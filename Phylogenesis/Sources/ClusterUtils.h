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

#ifndef __ClusterUtils_h__
#define __ClusterUtils_h__

#include <Cluster.h>
#include <SingleCluster.h>

namespace Victor { namespace Phylogenesis{

	/// Count number of sequences in the current cluster
	static int countSequences(Cluster* c){
		if(c->isSingle()){
			return 1;
		}
		int tot = 0;
		for(unsigned int i = 0; i < c->branches.size(); i++){
			tot+= countSequences(c->branches[i]);
		}
		return tot;
	}

	/// Calculate weight sequences in the current cluster
	static void calculatesWeights(Cluster* c, double weight){
		if(c->isSingle()){
			SingleCluster* cs = static_cast<SingleCluster*>(c);
			cs->sequence->weight = weight;
			return;
		}
		for(unsigned int i = 0; i < c->branches.size(); i++){
			calculatesWeights(c->branches[i], weight + (c->lengths[i]/countSequences(c->branches[i])));
		}
	}

	/// Get max weight in the cluster
	static double getMaxWeight(Cluster* c, double max){
		if(c->isSingle()){
			SingleCluster* cs = static_cast<SingleCluster*>(c);
			if(cs->sequence->weight  > max){
				max = cs->sequence->weight;
			}
		}
		for(unsigned int i = 0; i < c->branches.size(); i++){
			double m = getMaxWeight(c->branches[i], max);
			if(m > max){
				max = m;
			}
		}
		return max;
	}

	/// Normalizes weights in the current cluster
	static void normalizesWeights(Cluster* c, double norm){
		if(c->isSingle()){
			SingleCluster* cs = static_cast<SingleCluster*>(c);
			cs->sequence->normalizedWeight = cs->sequence->weight / norm;
			return;
		}
		for(unsigned int i = 0; i < c->branches.size(); i++){
			normalizesWeights(c->branches[i], norm);
		}
	}

	/// Returns the sum of all weights in the current cluster
	static double countWeights(Cluster* c){
		if(c->isSingle()){
			SingleCluster* cs = static_cast<SingleCluster*>(c);
			return cs->sequence->weight;
		}
		double tot = 0;
		for(unsigned int i = 0; i < c->branches.size(); i++){
			tot+= countWeights(c->branches[i]);
		}
		return tot;
	}

	/// Prints weights and normalized weights in the current cluster
	static void print(Cluster* c){
		if(c->isSingle()){
			SingleCluster* cs = static_cast<SingleCluster*>(c);
			string name = cs->sequence->name;
			name = name.substr(0,name.find(' '));
			cout<<name<<" "<<cs->sequence->weight<<" "<<cs->sequence->normalizedWeight<<endl;
			return;
		}
		for(unsigned int i = 0; i < c->branches.size(); i++){
			print(c->branches[i]);
		}
	}

}}

#endif
