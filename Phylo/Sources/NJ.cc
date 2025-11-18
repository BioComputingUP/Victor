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

#include <NJ.h>

using std::max;

namespace Victor {
namespace Phylo {
    typedef UnrootedTree::Node Node;


    NJ::~NJ() {
    }

    
    
    UnrootedTree &
    NJ::apply(const DistanceMatrix &matrix) const {
        UnrootedTree *tree = new UnrootedTree();
        DistanceMatrix d = matrix;
        set<string> OTUs = matrix.getOTUs();
        map<string, Node *> node_pool;


        // Initializes leaves
        for (auto OTU : OTUs) {
            Node *leaf = new Node(OTU);
            node_pool[OTU] = leaf;
            tree->addNode(*leaf);
        }


        // Repeats until every element in the matrix is consumed
        while (!d.isEmpty()) {
            // Calculates matrix Q, based on current distance matrix
            DistanceMatrix Q;
            OTUs = d.getOTUs();
            const size_t n = OTUs.size();
            for (auto i : OTUs) {
                Q.addOTU(i);
                for (auto j : OTUs) {
                    if (i == j) continue;
                    double q_ij, d_ik = 0.0, d_jk = 0.0;
                    for (auto k : OTUs){
                        d_ik += d(i, k);
                        d_jk += d(j, k);
                    }

                    q_ij = (n - 2) * d(i, j) - d_ik - d_jk;
                    q_ij = q_ij;
                    Q(i, j, q_ij);
                }
            }


            // Finds f, g such that Q(f, g) is minimum (i != j), then
            // joins them to a new node u
            const pair<string, string> min_pos = Q.getMinimumPosition();
            const string f = min_pos.first,
                         g = min_pos.second,
                         u = f  + "+" + g;
            Node *node = new Node();
            node_pool[u] = node;
            tree->addNode(*node);


            // Calculates distances between new node and merged ones. If
            // there are only two nodes left, operation is trivial
            if (2 == n) {
                node->addNeighbor(*node_pool[f], d(f, g) / 2.0);
                node->addNeighbor(*node_pool[g], d(f, g) / 2.0);
            } else {
                double delta_fu, delta_gu, d_fk = 0.0, d_gk = 0.0;;
                for (auto k : OTUs) {
                    d_fk += d(f, k);
                    d_gk += d(g, k);
                }
                delta_fu = max(0.5 * d(f, g) + (1.0 / (2.0 * (n - 2))) * (d_fk - d_gk), 0.0);
                delta_gu = max(d(f, g) - delta_fu, 0.0);
                node->addNeighbor(*node_pool[f], delta_fu);
                node->addNeighbor(*node_pool[g], delta_gu);
            }


            // Calculates distances from each OTU to the new node
            for (auto k : OTUs) {
                double d_uk = 0.5 * (d(f, k) + d(g, k) - d(f, g));
                d(u, k, d_uk);
            }


            // Removes "consumed" OTUs and adds newly generated one
            d.removeOTU(f)
             .removeOTU(g)
             .addOTU(u);
        }

        return *tree;
    }
}  // namespace Phylo
}  // namespace Victor