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

#ifndef _VICTOR_PHYLO_FITCHMARGOLIASH_H_
#define _VICTOR_PHYLO_FITCHMARGOLIASH_H_

#include <ClusteringAlgorithm.h>
#include <UnrootedTree.h>

namespace Victor {
namespace Phylo {

/** @brief Implements the Fitch-Margoliash clustering algorithm.
 * 
 * Implements Fitch-Margoliash as a functor.
 * Fitch-Margoliash algorithm is a bottom-up agglomerative method for
 * the creation of phylogenetic trees. This algorithm requires knowledge
 * about the distance between each pair of taxa (e.g. species or
 * sequences) in order to build the tree.
 *
 * For more information see:
 * - Contruction of Phylogenetic Trees (Fitch and Margoliash)
 *
 * This class follows the Strategy Design Pattern.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class FitchMargoliash: public ClusteringAlgorithm {
 public:
    /** @brief Destructor. */
    ~FitchMargoliash();


    
    /** @brief Builds a phylogenetic tree from given distance matrix.
     * @param[in]       matrix    Distance matrix
     * @return Phylogenetic tree
     */
    virtual UnrootedTree &apply(const DistanceMatrix &matrix) const;
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_FITCHMARGOLIASH_H_