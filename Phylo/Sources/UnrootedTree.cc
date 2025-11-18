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
#include <UnrootedTree.h>
#include <RootedTree.h>

using std::stringstream;

namespace Victor {
namespace Phylo {
    typedef UnrootedTree::Node Node;

    UnrootedTree UnrootedTree::emptyTree;
    Node Node::emptyNode;
    const double Node::NA = -1.0;
    size_t Node::next_ID = 0;

    ////////////////////////////////////////////////////////////////////
    // UnrootedTree::Node
    Node::Node(UnrootedTree *tree, const string &label)
    : tree(tree), label(label) {
        this->ID = getNextIdentifier();

        if (nullptr != tree) {
            tree->addNode(*this);
        }
    }



    Node::Node(const string &label) : Node(nullptr, label) {
    }



    bool
    Node::operator==(const Node &other) const {
        return getIdentifier() == other.getIdentifier();
    }


    bool
    Node::operator!=(const Node &other) const {
        return !(*this == other);
    }



    bool
    Node::isLeaf() const {
        return getSize() < 2;
    }



    bool
    Node::hasNeighbors() const {
        return !neighbors.empty();
    }



    bool
    Node::isNeighbor(const size_t identifier) const {
        return neighbors.find(identifier) != neighbors.end();
    }



    bool
    Node::isNeighbor(const Node &node) const {
        return isNeighbor(node.getIdentifier());
    }



    bool
    Node::hasPath(const size_t identifier, const Node *previous) const {
        if (isNeighbor(identifier)) {
            return true;
        }

        bool has_path = false;
        for (auto neighbor : getNeighbors()) {
            if (neighbor == previous) {
                continue;
            }

            has_path |= neighbor->hasPath(identifier, this);
        }

        return has_path;
    }



    bool
    Node::hasPath(const Node &node) const {
        return hasPath(node.getIdentifier());
    }



    bool
    Node::belongsTo(const UnrootedTree &tree) const {
        return this->tree == &tree || getTree() == tree;
    }



    size_t
    Node::getIdentifier() const {
        return ID;
    }



    size_t
    Node::getID() const {
        return getIdentifier();
    }



    size_t
    Node::getSize() const {
        return neighbors.size();
    }



    size_t
    Node::getNeighborsNumber() const {
        return getSize();
    }



    Node &
    Node::getNeighbor(const size_t identifier) const {
        return (isNeighbor(identifier))
             ? *(neighbors.find(identifier)->second.first)
             : emptyNode;
    }



    Node &
    Node::getClosestNeighbor() const {
        Node *closest   = const_cast<Node *>(this);
        double distance = std::numeric_limits<double>::max();

        for (auto neighbor : getNeighbors()) {
            if (getDistance(*neighbor) < distance) {
                distance = getDistance(*neighbor);
                closest  = neighbor;
            }
        }

        return *closest;
    }



    Node &
    Node::getNeighborClosestTo(const Node &target, const Node *avoid) const {
        Node *closest   = const_cast<Node *>(this);
        double distance = std::numeric_limits<double>::max();

        for (auto neighbor : getNeighbors()) {
            if (avoid != nullptr && *avoid == *neighbor) continue;
            if (neighbor->getDistance(target) < distance) {
                distance = neighbor->getDistance(target);
                closest  = neighbor;
            }
        }

        return *closest;
    }



    vector<Node *>
    Node::getNeighbors() const {
        vector<Node *> neighbors;

        for (auto entry : this->neighbors) {
            neighbors.push_back(entry.second.first);
        }

        return neighbors;
    }



    double
    UnrootedTree::Node::getDistance(const size_t identifier, const Node *previous) const {
        // If nodes are neighbors, returns distance directly
        if (isNeighbor(identifier)) {
            return neighbors.find(identifier)->second.second;
        }
        
        // If nodes are not neighbors, finds a path
        double distance = NA;
        vector<Node *> neighbors = getNeighbors();

        for (auto neighbor : neighbors) {
            double to_neighbor = getDistance(*neighbor),
                   path_length;

            // Avoids infinite loops
            if (previous == neighbor || NA == to_neighbor) {
                continue;
            }

            path_length = neighbor->getDistance(identifier, this);

            // Exits the loop if the path was found
            if (path_length != NA) {
                distance = to_neighbor + path_length;
                break;
            }
        }

        return distance;
    }



    double
    UnrootedTree::Node::getDistance(const Node &node) const {
        return getDistance(node.getIdentifier());
    }



    UnrootedTree &
    Node::getTree() const {
        return (nullptr != tree) ? *tree : emptyTree;
    }



    string
    Node::getLabel() const {
        return label;
    }



    Node &
    Node::setLabel(const string &label) {
        this->label = label;
        return *this;
    }


    Node &
    Node::setDistance(const size_t identifier, const double distance) {
        if (isNeighbor(identifier)) {
             neighbors.find(identifier)->second.second = distance;
        }
        return *this;
    }


    Node &
    Node::setDistance(Node &node, const double distance) {
        return setDistance(node.getIdentifier(), distance);
    }



    Node &
    Node::addNeighbor(Node &node, const double distance) {
        // If nodes do not belong to the same tree, nothing happens
        if (!node.belongsTo(getTree())) {
            return *this;
        }

        neighbors[node.getIdentifier()] = make_pair(&node, distance);
        node.neighbors[getIdentifier()] = make_pair(this, distance);

        return *this;
    }



    Node &
    Node::addBetween(Node &neighbor, Node &node, const double distance) {
        // If node do not belong to the same tree, or are not neighbors,
        // nothing happens
        if (!node.belongsTo(getTree()) || !isNeighbor(neighbor)) {
            return *this;
        }

        const double d_ac = distance,
                     d_bc = getDistance(neighbor) - d_ac;

        addNeighbor(node, d_ac);
        neighbor.addNeighbor(node, d_bc);
        removeNeighbor(neighbor);

        return *this;
    }



    Node &
    Node::removeNeighbor(Node &node) {
        neighbors.erase(node.getIdentifier());
        node.neighbors.erase(getIdentifier());

        return *this;
    }


    RootedTree &
    Node::asRootedTree() const {
        return getTree().asRootedTree(getIdentifier());
    }



    size_t
    Node::getNextIdentifier() {
        return next_ID++;
    }





    ////////////////////////////////////////////////////////////////////
    // UnrootedTree
    UnrootedTree::UnrootedTree() {
    }



    UnrootedTree::~UnrootedTree() {
    }



    bool
    UnrootedTree::operator==(const UnrootedTree &other) const {
        // If trees point the same memory area, they are trivially equal
        if (this == &other) {
            return true;
        }

        // If trees have a different number of nodes, they cannot be equal
        if (getSize() != other.getSize()) {
            return false;
        }

        // Trees are equal if their nodes are equal
        bool equal = true;
        for (auto &entry : nodes) {
            equal &= getNode(entry.first) == other.getNode(entry.first);
        }

        return equal;
    }



    bool
    UnrootedTree::operator!=(const UnrootedTree &other) const {
        return !(*this == other);
    }



    UnrootedTree &
    UnrootedTree::operator+(Node &node) {
        return addNode(node);
    }



    Node &
    UnrootedTree::operator[](const size_t identifier) const {
        return getNode(identifier);
    }



    bool
    UnrootedTree::isEmpty() const {
        return getSize() == 0;
    }



    bool
    UnrootedTree::hasNode(const size_t ID) const {
        return nodes.find(ID) != nodes.end();
    }



    bool
    UnrootedTree::hasNode(const Node &node) const {
        return hasNode(node.getIdentifier());
    }



    size_t
    UnrootedTree::getSize() const {
        return nodes.size();
    }



    Node &
    UnrootedTree::getNode(const size_t identifier) const {
        return hasNode(identifier)
             ? *(nodes.find(identifier)->second)
             : Node::emptyNode;
    }



    vector<Node *>
    UnrootedTree::getNodes() const {
        vector<Node *> nodes;

        for (auto entry : this->nodes) {
            nodes.push_back(entry.second);
        }

        return nodes;
    }



    vector<Node *>
    UnrootedTree::getLeaves() const {
        vector<Node *> leaves;

        for (auto entry : nodes) {
            Node *node = entry.second;
            if (node->isLeaf()) {
                leaves.push_back(node);
            }
        }

        return leaves;
    }



    DistanceMatrix &
    UnrootedTree::getDistanceMatrix() const {
        vector<Node *> leaves = getLeaves();
        DistanceMatrix *d = new DistanceMatrix();

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



    UnrootedTree &
    UnrootedTree::addNode(Node &node) {
        nodes[node.getIdentifier()] = &node;
        return *this;
    }



    RootedTree &
    UnrootedTree::midpointRoot() const {
        map<string, Node *> node_pool;

        UnrootedTree tree = *this;
        DistanceMatrix d  = tree.getDistanceMatrix();

        // Initializes node pool
        for (auto node : tree.getLeaves()) {
            node_pool[node->getLabel()] = node;
        }

        // Finds most distant nodes A and B
        pair<string, string> max = d.getMaximumPosition();
        Node *A = node_pool[max.first],
             *B = node_pool[max.second];
        double max_distance = d(max.first, max.second);

        // Finds adjacent nodes between which new root will be placed
        Node *current  = A,
             *previous = B;
        double path_distance = max_distance;

        while (path_distance >= max_distance / 2.0 && !current->isNeighbor(*B)) {
            Node *avoid   = previous;
            previous      = current;
            current       = &current->getNeighborClosestTo(*B, avoid);
            path_distance = current->getDistance(*B);
        }

        // Adds new root node between current and previous
        Node *root = new Node();
        tree.addNode(*root);
        current->addBetween(*previous, *root, max_distance / 2.0 - path_distance);

        return tree.asRootedTree(*root);
    }



    RootedTree &
    UnrootedTree::asRootedTree() const {
        return midpointRoot();
    }



    RootedTree &
    UnrootedTree::asRootedTree(const size_t root) const {
        return asRootedTree(getNode(root));
    }



    RootedTree &
    UnrootedTree::asRootedTree(Node &root, const Node *parent) const {
        RootedTree *root_node = new RootedTree(root.getLabel());

        // Sets distance from parent, if parent is not null
        if (nullptr != parent) {
            root_node->setDistance(root.getDistance(*parent));
        }

        // Recursively adds children
        for (auto &entry : root.neighbors) {
            Node *neighbor = entry.second.first;
            if (neighbor != parent) {
                root_node->addChild(asRootedTree(*neighbor, &root));
            }
        }

        return *root_node;
    }



    UnrootedTree &
    UnrootedTree::asUnrootedTree() const {
        return *(const_cast<UnrootedTree *>(this));
    }



    UnrootedTree &
    UnrootedTree::parseNewick(const string &input) {
        RootedTree tree;

        tree.parseNewick(input);

        *this = tree.asUnrootedTree();

        return *this;
    }



    string
    UnrootedTree::asNewick() const {
        return asRootedTree().asNewick();
    }



    UnrootedTree &
    UnrootedTree::accept(Visitor &visitor) {
        visitor.visit(*this);
        return *this;
    }
}  // namespace Phylo
}  // namespace Victor