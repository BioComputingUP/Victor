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

/** Returns union of sets.
 * @param[in]       A         First set
 * @param[in]       B         Second set
 * @return Set union between A and B: A U B
 */
template <class T>
set<T> set_union(const set<T> &A, const set<T> &B) {
    set<T> C = A;
    C.insert(B.begin(), B.end());
    return C;
}



namespace Victor {
namespace Phylo {
    UPGMA::~UPGMA() {
    }


    
    RootedTree &
    UPGMA::apply(const DistanceMatrix &matrix) const {
        map<string, RootedTree *> node_pool;
        map<string, double> cumulative;
        map<string, set<string>> components;
        RootedTree *tree;
        DistanceMatrix DM = matrix;
        set<string> OTUs = matrix.getOTUs();


        // Initializes auxilary data structures
        for (auto &label : OTUs) {
            node_pool[label] = new RootedTree(label);
            cumulative[label] = 0.0;
            components[label].insert(label);
        }

        while (!DM.isEmpty()) {
            // Takes minimum from distance matrix and computes new label and
            // distance
            pair<string, string> min = DM.getMinimumPosition();
            string i = min.first,
                   j = min.second,
                   new_label = i + "+" + j;
            double new_distance = DM(i, j) / 2.0;

            // Creates a new node for this joining
            tree = new RootedTree();
            tree->addChild(*node_pool[i]);
            tree->addChild(*node_pool[j]);

            // Sets distance between new joining node and children
            node_pool[i]->setDistance(new_distance - cumulative[i]);
            node_pool[j]->setDistance(new_distance - cumulative[j]);

            // Updates auxiliary data structures
            node_pool[new_label] = tree;
            cumulative[new_label] = new_distance;
            components[new_label] = set_union(components[i], components[j]);


            // Updates distance matrix with new distances
            OTUs = DM.getOTUs();
            for (auto &l : OTUs) {
                // Does not update entries about consumed labels
                if (l == i || l == j) {
                    continue;
                }

                // Computes average distance
                double average_distance = 0.0;
                for (auto &a : components[new_label]) {
                    for (auto &b : components[l]) {
                        average_distance += matrix(a, b);
                    }
                }
                average_distance /= (components[new_label].size() * components[l].size());

                // Updates matrix
                DM(l, new_label, average_distance);
            }

            // Adds new OTU to the distance matrix, and removes consumed ones
            DM.addOTU(new_label)
              .removeOTU(i)
              .removeOTU(j);
        }
        return *tree;
    }

}  // namespace Phylo
}  // namespace Victor