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

#ifndef _VICTOR_PHYLO_CLUSTERINGALGORITHM_H_
#define _VICTOR_PHYLO_CLUSTERINGALGORITHM_H_

#include <PhylogeneticAlgorithm.h>
#include <PhylogeneticTree.h>
#include <DistanceMatrix.h>

namespace Victor {
namespace Phylo {

/** @brief A clustering algorithm.
 *
 * Clustering algorithms are phylogenetic algorithms which produce
 * phylogenetic trees by using information about distances among
 * alignments.
 *
 * Clustering algorithms are implemented as functors from a distance matrix to
 * a phylogenetic tree.
 *
 * This class follows the Strategy Design Pattern.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class ClusteringAlgorithm: public PhylogeneticAlgorithm {
 public:
    /** @brief Builds a phylogenetic tree from given distance matrix.
     * @param[in]       matrix    Distance matrix
     * @return Phylogenetic tree
     */
    inline PhylogeneticTree &operator()(const DistanceMatrix &matrix) const {
        return apply(matrix);
    }

    /** @brief Builds a phylogenetic tree from given distance matrix.
     * @param[in]       matrix    Distance matrix
     * @return Phylogenetic tree
     */
    virtual PhylogeneticTree &apply(const DistanceMatrix &matrix) const = 0;
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_CLUSTERINGALGORITHM_