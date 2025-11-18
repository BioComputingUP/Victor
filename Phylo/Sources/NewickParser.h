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

#ifndef _VICTOR_PHYLO_NEWICKPARSER_H_
#define _VICTOR_PHYLO_NEWICKPARSER_H_

#include <RootedTree.h>

namespace Victor {
namespace Phylo {

/** @brief Implements a parser for the Newick format.
 * 
 * This parser recognizes the following LL(1) grammar:
 * \verbatim
 Newick   ::= Tree;
            | ;
 Tree     ::= (Tree Siblings) Node
            | string Length
            | :number
 Siblings ::= , Tree Siblings
            | 
 Node     ::= Label Length
 Label    ::= string
            | 
 Length   ::= :number
            |
 \endverbatim
 * Where string represents any alphanumeric sequence or any sequence
 * of characters enclosed by double quotes, and number represents a
 * floating point number.
 * 
 * @note There is no official newick grammar. The one proposed here was
 * built with reverse engeneering techniques.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class NewickParser {
 public:
    /** @brief Destructor. */
    virtual ~NewickParser();


    
    /** @brief Parses a Newick string.
     * @param[in]       input     Input string
     * @return Rooted tree
     */
    RootedTree &parse(const string &input);



 private:
    string input;       ///< Input string


    /** @brief Tells whether input string is empty.
     * @retval true     Input string is empty
     * @retval false    Input string is not empty
     */
    bool isEmpty() const;

    /** @brief tells whether next input character belong to given set.
     * @param[in]       character_set Set of characters (as string)
     * @retval True if next character belongs to given set
     * @retval False if next character does not belong to given set
     */
    bool isNext(string character_set = "") const;

    /** @brief Returns next character in the input string.
     * @retval Next character from input string
     */
    char getNext() const;

    /** @brief Consumes first input character.
     *
     * Removes first character from input and returns it.
     * @retval Next input character
     */
    char shift();

    /** @brief Reports a parsing error and exits. */
    void error(const string &message = "Unknown error") const;

    /** @brief Parses a Newick file.
     * @return Rooted tree encoded by the Newick file
     */
    RootedTree &parseNewick();

    /** @brief Parses a tree.
     * @return Rooted tree
     */
    RootedTree *parseTree();

    /** @brief Parses siblings (sub)trees.
     * @param[in, out]  siblings  Siblings parsed so far
     */
    void parseSiblings(vector<RootedTree *> *siblings);

    /** @brief Parses a node.
     * @return Node of a rooted tree
     */
    RootedTree *parseNode();

    /** @brief Parses a label.
     * @return Label of the node
     */
    string parseLabel();

    /** @brief Parses an arc length.
     * @return Arc length
     */
    double parseLength();

    /** @brief Parses a string.
     *
     * Strings may be either spaces-free sequences of character or any
     * sequence of character enclosed by double quotes.
     * @return String
     */
    string parseString();

    /** @brief Parses a number.
     * @return Number (floating point)
     */
    double parseNumber();
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_NEWICKPARSER_H_