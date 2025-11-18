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

#include <FitchMargoliash.h>
#include <algorithm>
#include <iterator>

using Victor::Phylo::DistanceMatrix;

/** @brief Returns distance between two clusters of nodes.
 * @param[in]           in_X      Set of (names of) nodes in X
 * @param[in]           in_O      Set of (names of) nodes in O
 * @param[in]           d         Reference distance matrix
 */
static double d_XO(
    const set<string> &in_X,
    const set<string> &in_O,
    const DistanceMatrix &d);

/** @brief Returns union of sets.
 * @param[in]       A         First set
 * @param[in]       B         Second set
 * @return Set union between A and B: A U B
 */
template <class T>
set<T> set_union(const set<T> &A, const set<T> &B);

/** @brief Returns difference of sets.
 *
 * This is a more comfortable inteface to the std set difference function
 * which, unfortunately, does not seems a set difference function at all.
 * @param[in]       A         First set
 * @param[in]       B         Second set
 * @return Set difference between A and B: A - B
 */
template <class T>
set<T> set_difference(const set<T> &A, const set<T> &B);



namespace Victor {
namespace Phylo {
    typedef UnrootedTree::Node Node;


    FitchMargoliash::~FitchMargoliash() {
    }



    UnrootedTree &
    FitchMargoliash::apply(const DistanceMatrix &matrix) const {
        UnrootedTree *tree = new UnrootedTree();
        DistanceMatrix d   = matrix;
        set<string> OTUs   = matrix.getOTUs();
        map<string, Node *> node_pool;
        map<string, size_t> size_pool;
        map<string, set<string>> sets_pool;


        // Initializes leaves
        for (auto &OTU : OTUs) {
            Node *leaf = new Node(OTU);
            node_pool[OTU] = leaf;
            size_pool[OTU] = 1;
            sets_pool[OTU].insert(OTU);
            tree->addNode(*leaf);
        }


        // Repeats until only three nodes are left
        while (d.getSize() > 3) {
            OTUs = d.getOTUs();
            set<string> others;

            // Finds closest pair A and B, then join them into a new
            // node R
            const pair<string, string> min_pos = d.getMinimumPosition();
            const string A = min_pos.first,
                         B = min_pos.second,
                         R = A + "+" + B;
            Node *r      = new Node();
            node_pool[R] = r;
            size_pool[R] = size_pool[A] + size_pool[B];
            sets_pool[R] = set_union(sets_pool[A], sets_pool[B]);
            others       = set_difference(matrix.getOTUs(), sets_pool[R]);
            tree->addNode(*r);


            // Calculates distances among A, B and the other nodes (O)
            const double d_AB = d(A, B),
                         d_AO = d_XO(sets_pool[A], others, matrix),
                         d_BO = d_XO(sets_pool[B], others, matrix);


            // Calculates distances betwen R and A (a), R and B (b), R
            // and the other nodes (c)
            const double a = std::max(0.5 * (d_AO + d_AB - d_BO), 0.0),
                         b = std::max(0.5 * (d_BO + d_AB - d_AO), 0.0);
                         //c = abs(0.5 * (d_AO + d_BO - d_AB));


            // Set distances among A, B, R
            r->addNeighbor(*node_pool[A], a);
            r->addNeighbor(*node_pool[B], b);


            // Updates distance matrix
            for (auto &K : OTUs) {
                if (K == A || K == B) continue;
                const double d_RK = 0.5 * (d(A, K) + d(B, K));
                d(R, K, d_RK);
            }


            // Removes "consumed" OTUs and adds newly generated one
            d.unsetDistance(A)
             .unsetDistance(B)
             .removeOTU(A)
             .removeOTU(B)
             .addOTU(R);
        }


        // If there are only three nodes left, the operation is trivial
        if (d.getSize() == 3) {
            vector<string> OTUs;
            for (auto &OTU : d.getOTUs()) {
                OTUs.push_back(OTU);
            }

            const string A = OTUs[0],
                         B = OTUs[1],
                         C = OTUs[2];
            const double a = std::max(0.5 * (d(A, B) + d(A, C) - d(B, C)), 0.0),
                         b = std::max(0.5 * (d(A, B) + d(B, C) - d(A, C)), 0.0),
                         c = std::max(0.5 * (d(A, C) + d(B, C) - d(A, B)), 0.0);

            Node *r = new Node();
            tree->addNode(*r);
            r->addNeighbor(*node_pool[A], a);
            r->addNeighbor(*node_pool[B], b);
            r->addNeighbor(*node_pool[C], c);
        }

        // If there were only two nodes in the matrix, operation is trivial
        if (d.getSize() == 2) {
            vector<string> OTUs;
            for (auto &OTU : d.getOTUs()) {
                OTUs.push_back(OTU);
            }
            
            const string A = OTUs[0],
                         B = OTUs[1];

            Node *r = new Node();
            tree->addNode(*r);
            r->addNeighbor(*node_pool[A], d(A, B) / 2.0);
            r->addNeighbor(*node_pool[B], d(A, B) / 2.0);
        }

        return *tree;
    }

}  // namespace Phylo
}  // namespace Victor



////////////////////////////////////////////////////////////////////////
// Support (local) functions
static double d_XO(
    const set<string> &in_X,
    const set<string> &in_O,
    const DistanceMatrix &d) {
    const size_t N_X = in_X.size(),
                 N_O = in_O.size();
    double d_xo = 0.0;

    for (auto &x : in_X) {
        for (auto &o : in_O) {
            d_xo += d(x, o);
        }
    }
    d_xo /= (N_X * N_O);

    return d_xo;
}



template <class T>
set<T> set_union(const set<T> &A, const set<T> &B) {
    set<T> C = A;
    C.insert(B.begin(), B.end());
    return C;
}



template <class T>
set<T> set_difference(const set<T> &A, const set<T> &B) {
    set<T> result;

    std::set_difference(
        A.begin(), A.end(), B.begin(), B.end(),
        std::inserter(result, result.end())
    );

    return result;
}