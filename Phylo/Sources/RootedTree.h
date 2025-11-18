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

#ifndef _VICTOR_PHYLO_ROOTEDTREE_H_
#define _VICTOR_PHYLO_ROOTEDTREE_H_

#include <iostream>
#include <string>
#include <vector>
#include <PhylogeneticTree.h>

using std::string;
using std::vector;

namespace Victor {
namespace Phylo {

/** @brief Implements a rooted phylogenetic tree
 *
 * A rooted tree is a phylogenetic tree where every node may have a
 * label, arcs have a length and every node has a parent (except the
 * root node).
 *
 * A rooted tree is implemented using a recursive structure, and a
 * (sub)tree is represented by its root node.
 *
 * This class follows the Visitor Design Pattern and uses Method
 * Cascading (through Method Chaining).
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class RootedTree: public PhylogeneticTree {
 public:
    /** Subtree is an alias for RootedTree. */
    typedef RootedTree Subtree;

    /** Node is an alias for RootedTree. */
    typedef RootedTree Node;


    /** Distance from parent node Not Available. */
    static const double NA;


    /** @brief Default constructor.
     * @param[in,out]   parent    Pointer to parent node
     * @param[in]       distance  Distance from parent node
     * @param[in]       label     Label of this node
     */
    RootedTree(
        Node *parent = nullptr,
        const double distance = NA,
        const string &label = "");

    /** @brief Constructor.
     * @param[in]       distance  Distance from parent node
     * @param[in]       label     Label of this node
     */
    explicit RootedTree(const double distance, const string &label = "");

    /** @brief Constructor.
     * @param[in]       label     Label of this node
     */
    explicit RootedTree(const string &label);

    /** @brief Destructor. */
    virtual ~RootedTree();



    /** @brief Equality test for rooted trees.
     *
     * Two rooted trees are equal if (they point to the same
     * memory area, or) they have identical labels, distance and children.
     * @param[in]       other     Other tree
     * @retval true     Trees are equal
     * @retval false    Trees are not equal
     */
    bool operator==(const Subtree &other) const;

    /** @brief Inequality test for rooted trees.
     *
     * Two rooted trees are equal if (they point to the same
     * memory area, or) they have identical labels, distance and children.
     * @param[in]       other     Other tree
     * @retval true     Trees are not equal
     * @retval false    Trees are equal
     */
    bool operator!=(const Subtree &other) const;

    /** @brief Tells whether this node is an ancestor of given node.
     * @param[in]       other     Given node
     * @retval true     This node is an ancestor of given one
     * @retval false    This node is not an ancestor of given one
     */
    bool operator<(const Node &other) const;

    /** @brief Tells whether this node is a descendant of given node.
     * @param[in]       other     Given node
     * @retval true     This node is a descendant of given one
     * @retval false    This node is not a descendant of given one
     */
    bool operator>(const Node &other) const;

    /** @brief Adds a child to this node.
     * @param[in, out]  child     Child to be added
     * @return This node itself
     */
    RootedTree &operator+(Node &child);

    /** @brief Returns child at given index.
     *
     * If index is invalid, returns an empty tree.
     * @param[in]       index     Index of the child
     * @return Child at given index
     */
    Subtree &operator[](const size_t index) const;



    /** @brief Tells whether this node is a leaf.
     * @retval true     This node is a leaf
     * @retval false    This node has at least one child
     */
    bool isLeaf() const;

    /** @brief Tells whether this node has siblings.
     * @retval true     This node has at least one sibling
     * @retval false    This node does not have any siblings
     */
    bool hasSiblings() const;

    /** @brief Tells whether this node has children.
     * @retval true     This node has at least one child
     * @retval false    This node does not have any child
     */
    bool hasChildren() const;

    /** @brief Tells whether this node is a root node.
     * @retval true     This node is a root node
     * @retval false    This node is not a root node
     */
    bool isRoot() const;

    /** @brief Tells whether this node is an ancestor of given node.
     * @param[in]       node      Given node
     * @retval true     This node is an ancestor of giiven node
     * @retval false    This node is not an ancestor of given node
     */
    bool isAncestor(const Node &node) const;

    /** @brief Tells whether this node is a descendant of given node.
     * @param[in]       node      Given node
     * @retval true     This node is a descendant of given node
     * @retval false    This node is not a descendant of given node
     */
    bool isDescendant(const Node &node) const;

    /** @brief Tells whether this node has the distance from parent node set.
     * @retval true     A distance was set for this node
     * @retval false    A distance was not set for this node
     */
    bool hasDistance() const;



    /** @brief Returns number of nodes in this (sub)tree.
     * @return Number of nodes in this (sub)tree
     */
    size_t getSize() const;

    /** @brief Returns leaves of this (sub)tree.
     * @return Leaves as a vector
     */
    vector<Node *> getLeaves() const;

    /** @brief Returns number of leaves in this (sub)tree.
     * @return Number of leaves in this (sub)tree
     */
    size_t getLeavesNumber() const;

    /** @brief Returns number of siblings of this node.
     * @return Number of siblings
     */
    size_t getSiblingsNumber() const;

    /** @brief Return number of children of this node.
     * @return Number of children
     */
    size_t getChildrenNumber() const;

    /** @brief Returns child at given index.
     *
     * By default, returns first child.
     * If index is invalid, returns an empty tree.
     * @param[in]       index     Index of the child to be returned
     * @return Child node
     */
    Subtree &getChild(const size_t index = 0) const;

    /** @brief Returns root node of the tree this node belongs to.
     * @return Root node
     */
    Node &getRoot() const;

    /** @brief Returns depth of this node.
     * @return Depth of this node
     */
    size_t getDepth() const;

    /** @brief Returns height of this node.
     * @return Height of this node
     */
    size_t getHeight() const;

    /** @brief Returns parent of this node.
     *
     * If this node is a root node, returns an empty tree.
     * @return Parent node
     */
    Node &getParent() const;

    /** @brief Returns distance from this node to its parent.
     *
     * Equivalently, returns the length of the arc.
     * @return Distance from this node to its parent
     */
    double getDistance() const;

    /** @brief Returns distance from this node to given one.
     *
     * @param[in]       node      Node
     * @return Distance from this node to given one.
     */
    double getDistance(const Node &node) const;

    /** @brief Returns total distance from this node to its root.
     * @return Cumulative distance from this node to the root
     */
    double getTotalDistance() const;

    /** @brief Returns maximum distance in this tree.
     * @return Maximum cumulative distance
     */
    double getMaxDistance() const;

    /** @brief Returns label of this node.
     * @return Label of this node
     */
    string getLabel() const;

    /** @brief Returns previous sibling of this node.
     *
     * If this node has no previous sibling, returns an empty tree.
     * @return Previous sibling
     */
    Node &getPreviousSibling() const;

    /** @brief Returns next sibling of this node.
     *
     * If this node has no next sibling, returns an empty tree.
     * @return Next sibling
     */
    Node &getNextSibling() const;

    /** @brief Returns matrix of distances among leaves.
     * @return Distance matrix
     */
    virtual DistanceMatrix &getDistanceMatrix() const;



    /** @brief Sets distance from this node to its parent.
     *
     * Equivalently, sets the length of the arc.
     * @param[in]       distance  New distance from this node to its parent
     * @return This node
     */
    RootedTree &setDistance(const double distance = NA);

    /** @brief Sets label of this node.
     * @param[in]       label     New label of this node
     * @return This node
     */
    RootedTree &setLabel(const string &label = "");



    /** @brief Adds a new child to this node.
     *
     * New child become the last child of this node.
     * @param[in, out]  child     New child
     * @return This node
     */
    Subtree &addChild(Subtree &child);

    /** @brief Returns a rooted version of this tree.
     * @return A rooted version of this tree
     */
    virtual RootedTree &asRootedTree() const;

    /** @brief Returns an unrooted version of this tree.
     * @return An unrooted version of this tree
     */
    virtual UnrootedTree &asUnrootedTree() const;

    /** @brief Parses a string in Newick format.
     * @param[in]       input     String in Newick format
     * @return This tree itself
     */
    virtual RootedTree &parseNewick(const string &input);

    /** @brief Returns a Newick version of this tree.
     * @return String representing this tree in Newick format
     */
    virtual string asNewick() const;



    /** @brief Accepts a visitor.
     * @param[in, out]  visitor   Visitor to accept
     * @return This tree itself
     */
    virtual RootedTree &accept(Visitor &visitor);

 private:
    static RootedTree emptyTree;  ///< Empty tree

    string label;                 ///< Label of this node
    double distance;              ///< Distance from parent node
    Node *parent;                 ///< Parent of this node
    vector<Subtree *> children;   ///< Children of this node
    Node *previous;               ///< Pointer to Previous sibling
    Node *next;                   ///< Pointer to next sibling
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_ROOTEDTREE_H_