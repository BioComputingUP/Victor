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

#include <NewickParser.h>

#define SYNTAX_ERROR -1

using std::stringstream;

namespace Victor {
namespace Phylo {
    static const double NA = -1.0;


    NewickParser::~NewickParser() {
    }


    
    bool
    NewickParser::isEmpty() const {
        return input.empty();
    }



    bool
    NewickParser::isNext(string character_set) const {
        bool is_next = false;

        if (isEmpty()) {
            return false;
        }

        for (auto c : character_set) {
            is_next |= c == input[0];
        }

        return is_next;
    }



    char
    NewickParser::getNext() const {
        return (!isEmpty()) ? input[0] : '\0';
    }



    char
    NewickParser::shift() {
        char c = getNext();

        if (!isEmpty()){
            input = input.substr(1);
        }

        return c;
    }



    void
    NewickParser::error(const string &message) const {
        std::cerr << "[Victor::Phylo::NewickParser Syntax error]: "
                  << message
                  << " near: \"" << input << "\""
                  << std::endl;
        exit(SYNTAX_ERROR);
    }



    RootedTree &
    NewickParser::parse(const string &input) {
        this->input = input;
        return parseNewick();
    }



    RootedTree &
    NewickParser::parseNewick() {
        RootedTree *newick;

        // Empty Newick file
        if (isNext(";")) {
            newick = new RootedTree();
        }

        // Recursvely parses file
        else {
            newick = parseTree();
        }

        // Consumes ";" and checks for errors.
        shift();
        if (!isEmpty()) {
            error("Stray character after last character");
        }

        return *newick;
    }



    RootedTree *
    NewickParser::parseTree() {
         RootedTree *node;

        // Subtree
        if (isNext("(")) {
            vector<RootedTree *> children;

            shift();
            children.push_back(parseTree());
            parseSiblings(&children);
            shift();
            node = parseNode();
            
            for (auto child : children) {
                node->addChild(*child);
            }
        }

        // Length only
        else if (isNext(":")) {
            double length = parseLength();
            node = new RootedTree(length);
        }

        // Label and length
        else {
            string label = parseString();
            double length = parseLength();

            node = (length != NA)
                 ? new RootedTree(length, label)
                 : new RootedTree(label);
        }

        return node;
    }



    void
    NewickParser::parseSiblings(vector<RootedTree *> *siblings) {
        // No more siblings
        if (isNext(")")) {
            return;
        }

        // One more sibling
        else if (isNext(",")) {
            shift();
            siblings->push_back(parseTree());
            parseSiblings(siblings);
            return;
        }

        // Syntax error
        else {
            error("Error while parsing siblings");
        }
    }



    RootedTree *
    NewickParser::parseNode() {
        string label = parseLabel();
        double length = parseLength();
        RootedTree *node = (length != NA)
                         ? new RootedTree(length, label)
                         : new RootedTree(label);
        return node;
    }



    string
    NewickParser::parseLabel() {
        return (!isNext(";):,")) ? parseString() : "";
    }



    double
    NewickParser::parseLength() {
        // A length is specified
        if (isNext(":")) {
            shift();
            return parseNumber();
        }

        // No length is specified
        else if (isNext(";),")) {
            return NA;
        }

        // Syntax error
        else {
            error("Error while parsing number");
            return NA;
        }
    }



    string
    NewickParser::parseString() {
        string buffer = "";

        // String is enclosed by double quotes
        if (isNext("\"")) {
            shift();
            while (!isEmpty() && !isNext("\"")) {
                buffer += shift();
            }
            shift();
        }

        // String is not enclosed by double quotes
        else {
            while (!isEmpty() && !isNext("):, ;")) {
                buffer += shift();
            }
        }

        return buffer;
    }



    double
    NewickParser::parseNumber() {
        char *ptr;
        double length = strtod(input.c_str(), &ptr);

        input = string(ptr);

        return length;
    }

}  // namespace Phylo
}  // namespace Victor