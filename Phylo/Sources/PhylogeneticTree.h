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

#ifndef _VICTOR_PHYLO_PHYLOGENETICTREE_H_
#define _VICTOR_PHYLO_PHYLOGENETICTREE_H_

#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <DistanceMatrix.h>
#include <Visitor.h>

using std::istream;
using std::ostream;
using std::string;
using std::vector;

namespace Victor {
namespace Phylo {

// Forward declaration
class RootedTree;
class UnrootedTree;

/** @brief A phylogenetic tree.
 *
 * A phylogenetic tree is a tree where every node may have a label, and
 * arcs have a length.
 * 
 * This class follows the Visitor Design Pattern and uses Method
 * Cascading (through Method Chaining).
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class PhylogeneticTree {
 public:
    /** @brief Destructor. */
    virtual ~PhylogeneticTree();


    
    /** @brief Returns matrix of distances among leaves.
     * @return Distance matrix
     */
    virtual DistanceMatrix &getDistanceMatrix() const = 0;

    /** @brief Returns a rooted version of this tree.
     * @return A rooted version of this tree
     */
    virtual RootedTree &asRootedTree() const = 0;

    /** @brief Returns an unrooted version of this tree.
     * @return An unrooted version of this tree
     */
    virtual UnrootedTree &asUnrootedTree() const = 0;

    /** @brief Parses a string in Newick format.
     * @param[in]       input     String in Newick format
     * @return This tree itself
     */
    virtual PhylogeneticTree &parseNewick(const string &input) = 0;

    /** @brief Returns a Newick version of this tree.
     * @return String representing this tree in Newick format
     */
    virtual string asNewick() const = 0;

    /** @brief Loads a phylogenetic tree in Newick format.
     * @param[in]       input     Stream to read input from
     */
    void loadNewick(const istream &input);

    /** @brief Saves a phylogenetic tree in Newick format.
     * @param[out]      output    Stream to write output to
     */
    void saveNewick(ostream &output) const;



    /** @brief Accepts a visitor.
     * @param[in, out]  visitor   Visitor to accept
     * @return This tree itself
     */
    virtual PhylogeneticTree &accept(Visitor &visitor);
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_PHYLOGENETICTREE_