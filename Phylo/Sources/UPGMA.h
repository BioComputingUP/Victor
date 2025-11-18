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

#ifndef _VICTOR_PHYLO_UPGMA_H_
#define _VICTOR_PHYLO_UPGMA_H_

#include <ClusteringAlgorithm.h>
#include <RootedTree.h>

namespace Victor {
namespace Phylo {

/** @brief Implements the UPGMA algorithm.
 *
 * Implements UPGMA as a functor.
 * UPGMA (Unweighted Pair Group Method with Arithmetic Mean) is an
 * agglomerative, bottom-up hierarchical clustering method.
 * The UPGMA algorithm constructs a rooted tree (dendrogram) that
 * reflects the structure present in a pairwise similarity matrix
 * (or a dissimilarity matrix).
 * 
 * For more information see:
 * - A quantitative approach to a problem in classification
 *   (Michanel and Sockal)
 *
 * This class follows the Strategy Design Pattern.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class UPGMA: public ClusteringAlgorithm {
 public:
    /** @brief Destructor. */
    virtual ~UPGMA();


    
    /** @brief Builds a phylogenetic tree from given distance matrix.
     * @param[in]       matrix    Distance matrix
     * @return Phylogenetic tree
     */
    virtual RootedTree &apply(const DistanceMatrix &matrix) const;
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_UPGMA_