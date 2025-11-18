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

#ifndef _VICTOR_PHYLO_SEQUENCE_H_
#define _VICTOR_PHYLO_SEQUENCE_H_

#define _GLIBCXX_IOMANIP
#include <string>
#include <map>
#include <AlignmentBase.h>

using std::string;
using std::map;
using Victor::Align2::AlignmentBase;

namespace Victor {
namespace Phylo {

/** @brief A sequence.
 *
 * Contains information about a sequence.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class Sequence {
 public:
    /** @brief Builds a sequence with given identifier and residues
     * @param[in]       identifier Idenfifier of this sequence
     * @param[in]       residues Residues in this sequence
     */
    explicit Sequence(const string &identifier, const string &residues);

    /** @brief Destructor. */
    ~Sequence();



    /** Converts an alignment into a list of sequences.
     * @param[in]       alignment Alignment
     * @return Vector of Sequence
     */
    static vector<Sequence> &toSequences(const AlignmentBase &alignment);



    /** Converts an alignment into a map of sequences.
     * @param[in]       alignment Alignment
     * @return Map of sequences: identifier => residues
     */
    static map<string, string> &toMap(const AlignmentBase &alignment);



    /** @brief Tells whether this sequence is empty.
     * @retval true     This sequence is empty
     * @retval false    This sequence is not empty
     */
    bool isEmpty() const;



    /** @brief Returns identifier of this sequence.
     * @return Identifier of this sequence.
     */
    string getIdentifier() const;

    /** @brief Returns identifier of this sequence.
     * @return Identifier of this sequence.
     */
    string getID() const;

    /** @brief Returns residue at given position.
     * @param[in]       position  Position of the residue to return
     * @return Residue
     */
    char getResidue(const size_t position) const;

    /** @brief Returns residues in this sequence.
     * @return Residues in this sequence
     */
    string getResidues() const;


    
 private:
    string identifier;  ///< Identifier of this sequence
    string residues;   ///< Residues in this sequence



    /** @brief Returns a compact version of an identifier.
     *
     * Strips spaces and keeps only the first word of the identifier.
     * @param[in]       identifier Identifier
     * @return Compact identifier
     */
    static string formatIdentifier(const string &identifier);
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_SEQUENCE_H_
