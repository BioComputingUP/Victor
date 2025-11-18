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

#ifndef _VICTOR_PHYLO_MULTIPLEALIGNMENT_H_
#define _VICTOR_PHYLO_MULTIPLEALIGNMENT_H_

#include <string>
#include <vector>
#include <map>
#include <Sequence.h>

using std::string;
using std::pair;
using std::vector;
using std::map;

namespace Victor {
namespace Phylo {

/** @brief A multiple sequence alignment (MSA).
 *
 * MSAs holds information about alignments of a set of sequences.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class MultipleAlignment {
public:
    /** @brief Builds an MSA from a set of sequences.
     * @param[in]       sequences Sequences in the MSA
     */
    explicit MultipleAlignment(const vector<Sequence> &sequences);

    /** @brief Builds an MSA from a single sequence.
     * @param[in]       sequence  Sequence
     */
    explicit MultipleAlignment(const Sequence &sequence);

    /** @brief Default constructor.
     *
     * Produces an empty MSA.
     */
    MultipleAlignment();

    /** @brief Destructor. */
    ~MultipleAlignment();



    /** @brief Assignment operator.
     * @param[in]       other     MultipleAlignment to copy
     */
    void operator=(const MultipleAlignment &other);



    /** @brief Tells whether the MSA is empty.
     *
     * An MSA is empty when it contains no sequences.
     * @retval true     This MSA is empty
     * @retval false    This MSA is not empty
     */
    bool isEmpty() const;

    /** @brief Tells whether this MSA contains a sequence with given identifier.
     * @param[in]       identifier An identifier
     * @retval true     This MSA contains data about given sequence
     * @retval false    This MSA does not contain data about given sequence
     */
    bool hasSequence(const string &identifier) const;

    /** @brief Tells whether this MSA has a hydrophilic stretch.
     * @param[in]       position  Position to check
     * @retval true     This MSA contains a hydrophilic stretch
     * @retval false     This MSA does not contain a hydrophilic stretch
     */
    bool hasHydrophilicStretch(const size_t position) const;



    /** @brief Returns number of sequences in this MSA.
     * @return Number of sequences in this MSA.
     */
    size_t getSize() const;

    /** @brief Returns length of the MSA.
     *
     * The length of an MSA is the length of its sequences (including gaps).
     * Every sequence has the same length, by construction.
     * @return Length of this multiple alignment
     */
    size_t getLength() const;

    /** @brief Returns sequence at given index.
     * 
     * If index is invalid, returns an empty sequence.
     * @param[in]       index     Index of the sequence
     * @return Sequence at given index
     */
    Sequence getSequence(const size_t index) const;

    /** @brief Returns index of given sequence.
     *
     * Searches a sequence with given identifier and returns its index in this
     * MSA. If sequence is not found, returns an invalid index (greater than
     * or equal to the size of this MSA).
     * @param[in]       identifier Identifier of the sequence
     * @return Index of the found sequence, or an invalid position
     */
    size_t getSequence(const string &identifier) const;

    /** @brief Returns residues at given position.
     *
     * Returns a string containing residues (or gaps) at given column of this
     * MSA.
     *
     * If position is invalid, returns an empty string.
     * @param[in]       position  Index of the columng.
     * @return Residues at given position
     */
    string getResidues(const size_t position) const;

    /** @brief Returns number of residues (or gaps) at given position.
     *
     * If position is invalid, returns 0.
     * @param[in]       position  Index of the column
     * @return Number of residues at given column of this MSA
     */
    size_t getResiduesNumber(const size_t position) const;

    /** @brief Returns number of non-gaps at given position.
     *
     * If position is invalid, returns 0.
     * @param[in]       position  Index of the column
     * @return Number of non-gaps at given column of this MSA
     */
    size_t getNonGapNumber(const size_t position) const;

    /** @brief Returns number of gaps at given position.
     *
     * If position is invalid, returns 0.
     * @param[in]       position  Index of the column
     * @return Number of gaps at given column of this MSA
     */
    size_t getGapNumber(const size_t position) const;

    /** @brief Returns number of occurencies of given residues at given position.
     *
     * If residue or position are invalid, returns 0.
     * @param[in]       residue  Residue to find
     * @param[in]       position  Index of the column
     * @return Number of occurencies of given residue at given position
     */
    size_t getCount(const char residue, const size_t position) const;

    /** @brief Returns frequency of given residue at given position.
     *
     * If residue or position are invalid, returns 0.
     * @param[in]       residue  Residue to find
     * @param[in]       position  Index of the column
     * @return Frequency of given residue at given position
     */
    double getFrequency(const char residue, const size_t position) const;

    /** @brief Returns the map of frequencies for given position.
     *
     * Keys of the map are residues (or gap), values are their frequencies.
     * 
     * If position is invalid, frequencies in the map are 0.
     * @param[in]       position  Index of the column
     * @return Map of frequencies
     */
    map<char, double> &getFrequencies(const size_t position) const;

    /** @brief Returns percent identity of this multiple alignment.
     * @return Percent identity
     */
    double getPercentIdentity() const;

    /** @brief Returns consensus residue.
     * @param[in]       position  Position to use
     * @return Consensus
     */
    char getConsensus(const size_t position) const;

    /** @brief Returns consensus sequence.
     * @return Consensus
     */
    string &getConsensus() const;



    /** @brief Returns a string version of this multiple alignment.
     * @return This multiple sequence alignment as a string
     */
    string &asString() const;

    /** @brief Exports this multiple sequence alignment.
     * @param[out]      output    Stream to write output to
     */
    void saveClustalW(ostream &output) const;



private:
    static const string residues;     ///< Every possible residue
    vector<Sequence> sequences;  ///< Sequences in the MSA
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_MULTIPLEALIGNMENT_H_