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

#include <sstream>
#include <limits>
#include <RootedTree.h>
#include <UnrootedTree.h>
#include <NewickParser.h>

using std::stringstream;

namespace Victor {
namespace Phylo {
    typedef RootedTree::Node Node;
    typedef RootedTree::Subtree Subtree;

    RootedTree RootedTree::emptyTree;
    const double RootedTree::NA = -1.0;

    RootedTree::RootedTree(Node *parent, const double distance, const string &label) 
    : parent(parent), label(label), distance(distance) {
        this->previous = nullptr;
        this->next     = nullptr;
        if (nullptr != parent) {
            parent->addChild(*this);
        }
    }

    RootedTree::RootedTree(const double distance, const string &label):
        RootedTree(nullptr, distance, label) {}

    RootedTree::RootedTree(const string &label):
        RootedTree(nullptr, NA, label) {}



    RootedTree::~RootedTree() {
    }



    bool
    RootedTree::operator==(const Subtree &other) const {
        bool equal = true;
        size_t i = 0;

        // If trees point the same memory area, they are equal.
        if (this == &other) {
            return true;
        }

        // If trees have different labels, number of children or distance,
        // they are not equal.
        if (label != other.label
            || getChildrenNumber() != other.getChildrenNumber()
            || getDistance() != other.getDistance()) {
            return false;
        }

        // If children of trees are equal, trees are equal.
        while (equal && i < getChildrenNumber()) {
            equal &= getChild(i) == other.getChild(i);
            i++;
        }

        return equal;
    }



    bool
    RootedTree::operator!=(const Subtree &other) const {
        return !(*this == other);
    }



    bool
    RootedTree::operator<(const Node &other) const {
        return isAncestor(other);
    }



    bool
    RootedTree::operator>(const Node &other) const {
        return isDescendant(other);
    }



    RootedTree &
    RootedTree::operator+(Node &child) {
        return addChild(child);
    }



    Subtree &
    RootedTree::operator[](const size_t index) const {
        return getChild(index);
    }



    bool
    RootedTree::isLeaf() const {
        return children.size() == 0;
    }



    bool
    RootedTree::hasSiblings() const {
        return nullptr != previous || nullptr != next;
    }



    bool
    RootedTree::hasChildren() const {
        return !isLeaf();
    }



    bool
    RootedTree::isRoot() const {
        return nullptr == parent;
    }



    bool
    RootedTree::isAncestor(const Node &node) const {
        if (node.isRoot()) {
            return false;
        } else if (node.parent == this) {
            return true;
        } else {
            return isAncestor(node.getParent());
        }
    }



    bool
    RootedTree::isDescendant(const Node &node) const {
        return node.isAncestor(*this);
    }



    bool
    RootedTree::hasDistance() const {
        return distance != NA;
    }




    size_t
    RootedTree::getSize() const {
        size_t size = 0;

        for (Node *child : children) {
            size += child->getSize();
        }

        return size + 1;
    }



    vector<Node *>
    RootedTree::getLeaves() const {
        vector<Node *> leaves;

        if (isLeaf()) {
            leaves.push_back(const_cast<Node *>(this));
        }
        else {
            for (Node *child : children) {
                vector<Node *> child_leaves = child->getLeaves();
                leaves.reserve(leaves.size() + child_leaves.size());
                leaves.insert(leaves.end(), child_leaves.begin(), child_leaves.end());
            }
        }

        return leaves;
    }



    size_t
    RootedTree::getLeavesNumber() const {
        size_t leaves = (isLeaf()) ? 1 : 0;
        
        for (Node *child : children) {
            leaves += child->getLeavesNumber();
        }

        return leaves;
    }



    size_t
    RootedTree::getSiblingsNumber() const {
        return isRoot() ? 0 : getParent().getChildrenNumber() - 1;
    }



    size_t
    RootedTree::getChildrenNumber() const {
        return children.size();
    }



    Node &
    RootedTree::getChild(const size_t index) const {
        return (index < getChildrenNumber())
             ? *(children[index]) : emptyTree;
    }



    Node &
    RootedTree::getRoot() const {
        return isRoot()
             ? *const_cast<RootedTree *>(this)
             : getParent().getRoot();
    }



    size_t
    RootedTree::getDepth() const {
        return isRoot() ? 0 : getParent().getDepth() + 1;
    }



    size_t
    RootedTree::getHeight() const {
        if (isLeaf()) {
            return 0;
        } else {
            size_t height = 0, max_height = 0;
            for (Node *child : children) {
                height = child->getHeight();
                if (height > max_height) {
                    max_height = height;
                }
            }
            return height + 1;
        }
    }



    Node &
    RootedTree::getParent() const {
        return isRoot() ? emptyTree : *parent;
    }



    double
    RootedTree::getDistance() const {
        return distance;
    }



    double
    RootedTree::getDistance(const Node &node) const {
        if (isAncestor(node)) {
            return node.getTotalDistance() - getTotalDistance();
        } else if (isDescendant(node)) {
            return getTotalDistance() - node.getTotalDistance();
        } else {
            Node *ancestor = parent;
            while (!ancestor->isAncestor(node)) {
                ancestor = ancestor->parent;
            }
            return getDistance(*ancestor) + node.getDistance(*ancestor);
        }
    }



    double
    RootedTree::getTotalDistance() const {
        return isRoot()
             ? getDistance()
             : getDistance() + getParent().getTotalDistance();
    }



    double
    RootedTree::getMaxDistance() const {
        double max_distance = 0.0;

        for (auto &leaf : getLeaves()) {
            max_distance = std::max(max_distance, leaf->getTotalDistance());
        }

        return max_distance;
    }



    string
    RootedTree::getLabel() const {
        return label;
    }



    Node &
    RootedTree::getPreviousSibling() const {
        return (nullptr != previous) ? *previous : emptyTree;
    }



    Node &
    RootedTree::getNextSibling() const {
        return (nullptr != next) ? *next : emptyTree;
    }



    DistanceMatrix &
    RootedTree::getDistanceMatrix() const {
        DistanceMatrix *d = new DistanceMatrix();
        vector<Node *> leaves = getLeaves();

        for (auto i : leaves) {
            d->addOTU(i->getLabel());
            for (auto j : leaves) {
                if (i == j) {
                    continue;
                }
                (*d)(i->getLabel(), j->getLabel(), i->getDistance(*j));
            }
        }

        return *d;
    }



    RootedTree &
    RootedTree::setDistance(const double distance) {
        this->distance = distance;
        return *this;
    }



    RootedTree &
    RootedTree::setLabel(const string &label) {
        this->label = label;
        return *this;
    }



    Subtree &
    RootedTree::addChild(RootedTree::Subtree &child) {
        if (hasChildren()) {
            children[getChildrenNumber() - 1]->next = &child;
            child.previous = children[getChildrenNumber() - 1];
        }
        children.push_back(&child);
        child.parent = this;
        return *this;
    }



    RootedTree &
    RootedTree::asRootedTree() const {
        return *(const_cast<RootedTree *>(this));
    }



    /** @brief Recursively builds an unrooted tree.
     * @param[in]       root      Current node
     * @param[in]       tree      Pointer to unrooted tree built so far
     * @return Unrooted version of this node
     * @note This local utility function must be implemented out of the
     * RootedTree class, because UnrootedTree::Node is not available
     * in the class definition. This is due to a missing feature in C++:
     * no forward declarations for nested class.
     */
    static UnrootedTree::Node &
    _buildNode(const RootedTree &root, UnrootedTree *tree = nullptr) {
        // Creates tree if null
        if (nullptr == tree) {
            tree = new UnrootedTree();
        }

        // Adds current node to the tree
        UnrootedTree::Node *node =
            new UnrootedTree::Node(tree, root.getLabel());

        // Recursively adds children of current node to the tree as neighbors
        for (size_t i = 0; i < root.getChildrenNumber(); i++) {
            const RootedTree child = root.getChild(i);
            node->addNeighbor(_buildNode(child, tree), child.getDistance());
        }

        return *node;
    }

    UnrootedTree &
    RootedTree::asUnrootedTree() const {
        return _buildNode(*this).getTree();
    }



    RootedTree &
    RootedTree::parseNewick(const string &input) {
        NewickParser parser;

        *this = parser.parse(input);

        return *this;
    }

    string
    RootedTree::asNewick() const {
        stringstream output;

        // Saves data about children
        if (hasChildren()) {
            output << "(" << getChild(0).asNewick() << ")";
        }

        // Saves label of this node
        output << getLabel();

        // Saves data about distance
        if (hasDistance()) {
            output << ":" << getDistance();
        }

        // Saves data about sibling
        if (next != nullptr) {
            output << "," << next->asNewick();
        }

        return output.str();
    }



    RootedTree &
    RootedTree::accept(Visitor &visitor) {
        visitor.visit(*this);
        return *this;
    }
}  // namespace Phylo
}  // namespace Victor