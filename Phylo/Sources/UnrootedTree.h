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

#ifndef _VICTOR_PHYLO_UNROOTEDTREE_H_
#define _VICTOR_PHYLO_UNROOTEDTREE_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <PhylogeneticTree.h>

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::make_pair;

namespace Victor {
namespace Phylo {

/** @brief Implements an unrooted tree.
 *
 * An unrooted tree is a phylogenetic tree where every node may have a
 * label, arcs have a length and nodes may have neighbors.
 *
 * This class follows the Visitor Design Pattern and uses Method
 * Cascading (through Method Chaining).
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class UnrootedTree: public PhylogeneticTree {
 public:
    /** @brief Implements a node in an unrooted tree. */
    class Node {
        friend class UnrootedTree;
     public:
        /** @brief Builds a node and inserts it into a tree.
         * @param[in, out] tree      Tree this node belongs to
         * @param[in]      label     Label of this node
         */
        explicit Node(UnrootedTree *tree, const string &label = "");

        /** @brief Default constructor.
         * @param[in]      label     Label of this node
         */
        Node(const string &label = "");



        /** @brief Tells whether two nodes are equal.
         *
         * Two nodes are equal when they have the same identifier.
         * @param[in]      other     Other node
         * @retval true    Nodes are equal
         * @retval false   Nodes are not equal
         */
        bool operator==(const Node &other) const;

        /** @brief Tells whether two nodes are not equal.
         *
         * Two nodes are equal when they have the same identifier.
         * @param[in]      other     Other node
         * @retval true    Nodes are not equal
         * @retval false   Nodes are equal
         */
        bool operator!=(const Node &other) const;



        /** @brief Tells whether this node is a leaf.
         * @retval true    This node is a leaf
         * @retval false   This node has at least two neighbors
         */
        bool isLeaf() const;

        /** @brief Tells whether this node has neighbors.
         * @retval true    This node has at least one neighbor
         * @retval false   This node does not have any neighbors
         */
        bool hasNeighbors() const;

        /** @brief Tells whether this node is a neighbor of given node.
         * @param[in]      identifier Identifier of given node
         * @retval true    This node is a neighbor of giiven node
         * @retval false   This node is not a neighbor of given node
         */
        bool isNeighbor(const size_t identifier) const;

        /** @brief Tells whether this node is a neighbor of given node.
         * @param[in]      node      Given node
         * @retval true    This node is a neighbor of given node
         * @retval false   This node is not a neighbor of given node
         */
        bool isNeighbor(const Node &node) const;

        /** @brief Tells whether a path exists from this node to given node.
         * @param[in]      identifier Identifier of given node
         * @param[in]      previous  Previously visited node
         * @retval true    A path between this node and given one exists
         * @retval false   A path between this node and given one does not exist
         */
        bool hasPath(const size_t identifier, const Node *previous = nullptr) const;

        /** @brief Tells whether a path exists from this node to given node.
         * @param[in]      node      Given node
         * @retval true    A path between this node and given one exists
         * @retval false   A path between this node and given one does not exist
         */
        bool hasPath(const Node &node) const;

        /** @brief Tells whether this node belongs to given tree
         * @param[in]      tree      Given tree
         * @retval true    This node belongs to given tree
         * @retval false   This node does not belong to given tree
         */
        bool belongsTo(const UnrootedTree &tree) const;



        /** @brief Returns identifier of this node.
         * @return Identifier of this node
         */
        size_t getIdentifier() const;

        /** @brief Returns identifier of this node.
         * @return Identifier of this node
         */
        size_t getID() const;

        /** @brief Returns number of neighbors of this node.
         * @return Number of neighbors of this node
         */
        size_t getSize() const;

        /** @brief Returns number of neighbors of this node.
         * @return Number of neighbors of this node
         */
        size_t getNeighborsNumber() const;

        /** @brief Returns neighbor at given index.
         *
         * If node is not a neighbor, returns an empty node.
         * @param[in]      identifier Identifier of the neighbor
         * @return Child node
         */
        Node &getNeighbor(const size_t identifier) const;

        /** @brief Returns neighbor closest to this node.
         * 
         * If more than one neighbor share minimal distance, choses one
         * in non-deterministic way.
         *
         * If this node has no neighbors, returns this node itself.
         * @return Closest neighbor
         */
        Node &getClosestNeighbor() const;

        /** @brief Returns neighbor closest to target node.
         *
         * If more than one neighbor share minimal distance, choses one
         * in non-deterministic way.
         *
         * If this node has no neighbors, or target node is not reachable,
         * returns this node.
         * @param[in]   target    Target node
         * @param[in]   avoid     Node to be avoided 
         * @return Neighbor closest to target
         */
        Node &getNeighborClosestTo(const Node &target, const Node *avoid = nullptr) const;

        /** @brief Returns neighbors of this node.
         * @return Neighbors as a vector
         */
        vector<Node *> getNeighbors() const;

        /** @brief Returns distance from this node to given one.
         * @param[in]      identifier Identifier of the node
         * @param[in]      previous   Neighbor which has been already visited
         * @return Distance from this node to given one.
         */
        double getDistance(const size_t identifier, const Node *previous = nullptr) const;

        /** @brief Returns distance from this node to given one.
         * @param[in]      node      Node
         * @return Distance from this node to given one.
         */
        double getDistance(const Node &node) const;

        /** @brief Returns unrooted tree this node belongs to.
         * @return Unrooted tree this node belongs to
         */
        UnrootedTree &getTree() const;

        /** @brief Returns label of this node.
         * @return Label of this node
         */
        string getLabel() const;



        /** @brief Sets label of this node.
         * @param[in]      label     New label of this node
         * @return This node
         */
        Node &setLabel(const string &label = "");

        /** @brief Sets distance between this node and given one.
         *
         * If given node is not a neighbor of this node, nothing happens.
         * @param[in, out] identifier Identifier of given node
         * @param[in]      distance  New distance
         * @return This node
         */
        Node &setDistance(const size_t identifier, const double distance = NA);

        /** @brief Sets distance between this node and given one.
         *
         * If given node is not a neighbor of this node, nothing happens.
         * @param[in, out] node      Given node
         * @param[in]      distance  New distance
         * @return This node
         */
        Node &setDistance(Node &node, const double distance = NA);



        /** @brief Adds a neighbor to this node.
         *
         * Neighbor node must belong to the same tree. If it does not
         * belong to the same tree, nothing happens.
         * @param[in, out] node      New neighbor
         * @param[in]      distance  Distance of new neighbor
         * @return This node
         */
        Node &addNeighbor(Node &node, const double distance = NA);

        /** @brief Adds a node between this node and one neighbor of its.
         *
         * If new node does not belong to the tree, or given node is
         * not a neighbor, nothing happens.
         * @param[in]      neighbor  Old neighbor
         * @param[in]      node      New neighbor
         * @param[in]      distance  Distance from this node
         * @return This node itself
         */
        Node &addBetween(Node &neighbor, Node &node, const double distance);

        /** @brief Removes a neighbor from this node.
         *
         * If node is not a neighbor of this node, nothing happens.
         * @param[in, out] node      Neighbor to remove
         * @return This node
         */
        Node &removeNeighbor(Node &node);


        /** @brief Returns a rooted tree rooted in this node.
         * @return A rooted tree version of the tree this node belongs to
         */
        RootedTree &asRootedTree() const;



     private:
        static Node emptyNode;    ///< Empty node
        static size_t next_ID;    ///< Next available identifier for a node
        static const double NA;   ///< Distance Not Available

        size_t ID;                ///< Identifier of this node
        string label;             ///< Label of this node
        map<size_t, pair<Node *, double>> neighbors;  ///< Neighbors of this node
        UnrootedTree *tree;       ///< Unrooted tree this node belongs to


        /** Returns next available identifier for a node.
         * @return Next available identifier
         */
        static size_t getNextIdentifier();
    };



    /** @brief Default constructor.
     *
     * Produces an empty unrooted tree.
     */
    UnrootedTree();

    /** @brief Destructor. */
    virtual ~UnrootedTree();



    /** @brief Equality test for unrooted trees.
     *
     * Two unrooted trees are equal if (they point to the same
     * memory area, or) they have the same nodes.
     * @param[in]       other     Other tree
     * @retval true     Trees are equal
     * @retval false    Trees are not equal
     */
    bool operator==(const UnrootedTree &other) const;

    /** @brief Inequality test for unrooted trees.
     *
     * Two unrooted trees are equal if (they point to the same
     * memory area, or) they have the same nodes.
     * @param[in]       other     Other tree
     * @retval true     Trees are not equal
     * @retval false    Trees are equal
     */
    bool operator!=(const UnrootedTree &other) const;

    /** @brief Adds a node to this tree.
     * @param[in, out]  node      Node to be added
     * @return This tree itself
     */
    UnrootedTree &operator+(Node &node);

    /** @brief Returns node with given identifier.
     *
     * If indentifier is invalid, returns an empty node.
     * @param[in]       identifier Index of the node
     * @return Node with given identifier
     */
    Node &operator[](const size_t identifier) const;



    /** @brief Tells whether this tree is empty.
     * @retval true     This tree is empty
     * @retval false    There is at least one node in this tree
     */
    bool isEmpty() const;

    /** @brief Tells whether a node belongs to this tree.
     * @param[in]       ID        Identifier of given node
     * @retval true     Node belongs to this tree
     * @retval false    Node does not belong to this tree
     */
    bool hasNode(const size_t ID) const;

    /** @brief Tells whether a node belongs to this tree.
     * @param[in]       node      Given node
     * @retval true     Node belongs to this tree
     * @retval false    Node does not belong to this tree
     */
    bool hasNode(const Node &node) const;



    /** @brief Returns number of nodes in this unrooted tree
     * @return Number of nodes in this tree
     */
    size_t getSize() const;

    /** @brief Returns node with given identifier.
     *
     * If identifier is invalid, returns an empty node.
     * @param[in]       identifier Identifier of the node
     * @return Node with given identifier
     */
    Node &getNode(const size_t identifier) const;

    /** @brief Returns every node in this tree.
     * @return Nodes as a vector
     */
    vector<Node *> getNodes() const;

    /** @brief Returns leaves in this tree.
     * @return Leaves as a vector
     */
    vector<Node *> getLeaves() const;

    /** @brief Returns matrix of distances among leaves.
     * @return Distance matrix
     */
    virtual DistanceMatrix &getDistanceMatrix() const;



    /** @brief Adds a node into this tree.
     * @param[in]      node      Node to insert
     * @return This tree itself
     */
    UnrootedTree &addNode(Node &node);



    /** @brief Returns a rooted version of this tree.
     * 
     * Root node is chosen using the midpoint rooting strategy.
     * @return A rooted version of this tree
     */
    virtual RootedTree &midpointRoot() const;

    /** @brief Returns a rooted version of this tree.
     *
     * Root node is chosen with a midpoint strategy.
     * @return A rooted version of this tree
     */
    virtual RootedTree &asRootedTree() const;

    /** @brief Returns a rooted version of this tree.
     *
     * Tree is rooted in one node of its.
     * @param[in]       root      Identifier of the root node
     * @return A rooted version of this tree
     */
    virtual RootedTree &asRootedTree(const size_t root) const;

    /** @brief Returns a rooted version of this tree.
     *
     * Tree is rooted in one node of its.
     * @param[in]       root      Node that will become root
     * @param[in]       parent    Pointer to parent node
     * @return A rooted version of this tree
     */
    virtual RootedTree &asRootedTree(Node &root, const Node *parent = nullptr) const;

    /** @brief Returns an unrooted version of this tree.
     * @return An unrooted version of this tree
     */
    virtual UnrootedTree &asUnrootedTree() const;

    /** @brief Parses a string in Newick format.
     * @param[in]       input     String in Newick format
     * @return This tree itself
     */
    virtual UnrootedTree &parseNewick(const string &input);

    /** @brief Returns a Newick version of this tree.
     * @return String representing this tree in Newick format
     */
    virtual string asNewick() const;



    /** @brief Accepts a visitor.
     * @param[in, out]  visitor   Visitor to accept
     * @return This tree itself
     */
    virtual UnrootedTree &accept(Visitor &visitor);



 private:
    /** Empty tree. */
    static UnrootedTree emptyTree;

    /** Nodes in this tree. */
    map<const size_t, Node *> nodes;
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_UNROOTEDTREE_H_