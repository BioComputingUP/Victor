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

#ifndef __Cluster_H__
#define __Cluster_H__

#include <vector>

using namespace std;

namespace Victor { namespace Phylogenesis{

	/** @brief    It represents the group of one or more clusters (recursive).
	 *			  Given a clustering, each Cluster is a node of the clustering tree.
	 *	@author   Andrea Nasari (andrea.nasari\@studenti.unipd.it)
	 *
	 **/
	class Cluster{
	public:

		/// Default constructor.
		Cluster() : distance(0), rooted(false){
		}

		/// Destructor.
		virtual ~Cluster(){
		}

		/// Check if the current cluster not contains no other clusters.
		virtual bool isSingle();

		/// Check if the cluster is rooted
		bool isRooted();

		/// Set if the cluster is rooted
		void setRooted(bool rt);

		// ATTRIBUTES:

		/// Pointers to the successors of of the current cluster.
		std::vector<Cluster*> branches;

		/// Length of the branches from the current cluster to the successors.
		std::vector<double> lengths;

		/// Distance of the current cluster.
		double distance;

		/// Indices of the sequences contained in the current cluster.
		std::vector<int> indices;

		/// Reroot the current tree with midpoint method
		static Cluster* midpointRooting(Cluster* tree,  int leafsnum);

	private:

		// ATTRIBUTES

		/// Cluster is rooted
		bool rooted;

		/// Directions used for root the tree
		std::vector<int> dirs;

		/// Current cluster's leafs
		static vector<Cluster*> leafs;

		// Static methods for making the rooted tree by midpoint method

		static void initUnrootedTree(Cluster* cluster, Cluster* parent, double dist, int nseqs);

		static double maxFarthestLeaf(Cluster* cluster, Cluster* parent, double dist, int id);

		static void removeParent(Cluster* cluster, Cluster* parent);

		static void removeParentRec(Cluster* cluster, Cluster* parent);

		static Cluster* addHalfwayRoot(Cluster* cluster, Cluster* parent, double dist, int id, double avg, double lastDist);
	};

	inline bool Cluster::isSingle(){
		return false;
	}

	inline bool Cluster::isRooted(){
		return rooted;
	}

	inline void Cluster::setRooted(bool rt){
		rooted = rt;
	}

}}

#endif
