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

#ifndef _VICTOR_PHYLO_VISITOR_H_
#define _VICTOR_PHYLO_VISITOR_H_

namespace Victor {
namespace Phylo {

// Forward class declaration
class PhylogeneticTree;
class RootedTree;
class UnrootedTree;

/** @brief A visitor for the phylogenetic tree hierarchy.
 * 
 * This class defines the interface of a visitor for the phylogenetic
 * tree hierarchy: PhylogeneticTree, RootedTree and UnrootedTree.
 *
 * This class follows the Visitor Design Pattern.
 *
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class Visitor {
 public:
    /** @brief Destructor. */
    virtual ~Visitor() = 0;



    /** @brief Visits a phylogenetic tree.
     * @param[in, out]  tree      Phylogenetic tree to visit
     */
    virtual void visit(PhylogeneticTree &tree) = 0;

    /** @brief Visits a rooted tree.
     * @param[in, out]  tree      Rooted tree to visit
     */
    virtual void visit(RootedTree &tree) = 0;

    /** @brief Visits an unrooted tree.
     * @param[in, out]  tree      Unrooted tree to visit
     */
    virtual void visit(UnrootedTree &tree) = 0;
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_VISITOR_H_