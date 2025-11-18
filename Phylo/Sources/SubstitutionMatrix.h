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

#ifndef _VICTOR_PHYLO_SUBSTITUTIONMATRIX_H_
#define _VICTOR_PHYLO_SUBSTITUTIONMATRIX_H_

#include <map>
#include <SubMatrix.h>

using std::map;
using Victor::Align2::SubMatrix;

namespace Victor {
namespace Phylo {

/** @brief Implements a substitution matrix.
 *
 * This class is similar to SubMatrix in the Align2 package. This class
 * can interoperate with SubMatrix, since it offers methods to convert
 * a SubMatrix to a SubstitutionMatrix and vice-versa.
 *
 * This class uses Method Cascading (thought Method Chaining) for its
 * setter member functions.
 *
 * Most common substitution matrices are hardcoded, so that it is no
 * longer necessary to read from a stream.
 *
 * Values for the precalculated matrices are taken from the official
 * ClustalW implementation.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class SubstitutionMatrix {
 public:
    /** Identifiers of hardcoded substitution matrices. */
    typedef enum {
        BLOSUM30,       ///< BLOSUM 30
        BLOSUM35,       ///< BLOSUM 35
        BLOSUM40,       ///< BLOSUM 40
        BLOSUM45,       ///< BLOSUM 45
        BLOSUM50,       ///< BLOSUM 50
        BLOSUM55,       ///< BLOSUM 55
        BLOSUM62,       ///< BLOSUM 62
        BLOSUM65,       ///< BLOSUM 65
        BLOSUM70,       ///< BLOSUM 70
        BLOSUM75,       ///< BLOSUM 75
        BLOSUM80,       ///< BLOSUM 80
        BLOSUM90,       ///< BLOSUM 90

        PAM20,          ///< PAM 20
        PAM60,          ///< PAM 60
        PAM120,         ///< PAM 120
        PAM160,         ///< PAM 160
        PAM250,         ///< PAM 250
        PAM350,         ///< PAM 350

        MD40,           ///< MD 40
        MD120,          ///< MD 120
        MD250,          ///< MD 250
        MD350,          ///< MD 350

        IDENTITY,       ///< Identity matrix

        GON40,          ///< GON 40
        GON80,          ///< GON 80
        GON120,         ///< GON 120
        GON160,         ///< GON 160
        GON250,         ///< GON 250
        GON300,         ///< GON 300
        GON350          ///< GON 350
    } Identifier;



    /** @brief Builds a subsitution matrix from a map.
     * @param[in]       scores    Scores among residues
     */
    explicit SubstitutionMatrix(const map<char, map<char, int>> &scores);

    /** @brief Builds a pre-calculated substitution matrix.
     * @param[in]       identifer Identifier of the matrix
     */
    explicit SubstitutionMatrix(const Identifier identifer);

    /** @brief Builds a substitution matrix from a SubMatrix.
     * 
     * Converts a SubMatrix to a SubstitutionMatrix.
     * @param[in]       matrix    SubMatrix object
     */
    explicit SubstitutionMatrix(const SubMatrix &matrix);

    /** @brief Copy constructor.
     * @param[in]       matrix    Substitution matrix
     */
    explicit SubstitutionMatrix(const SubstitutionMatrix &matrix);

    /** @brief Destructor. */
    virtual ~SubstitutionMatrix();



    /** @brief Tells whether this matrix is equal to given one.
     * @param[in]       other     Other substitution matrix
     */
    bool operator==(const SubstitutionMatrix &other) const;

    /** @brief Returns score between given residues.
     * @param[in]       a         First amino acid
     * @param[in]       b         Second amino acid
     * @return Score between a and b
     */
    int operator()(const char a, const char b) const;

    /** @brief Computes sum of scores and a given constant.
     * @param[in]       value     Value to be added
     * @return A new substitution matrix
     */
    SubstitutionMatrix &operator+(const int value) const;

    /** @brief Computes difference between scores and a given constant.
     * @param[in]       value     Value to be subtracted
     * @return A new substitution matrix
     */
    SubstitutionMatrix &operator-(const int value) const;

    /** @brief Changes signess of scores.
     * @return A new substitution matrix
     */
    SubstitutionMatrix &operator-() const;

    /** @brief Computes product of scores and a given constant.
     * @param[in]       value     Value to be multiplied
     * @return A new substitution matrix
     */
    SubstitutionMatrix &operator*(const int value) const;



    /** @brief Tells whether this matrix is equal to given one.
     * @param[in]       other     Other substitution matrix
     */
    bool isEqual(const SubstitutionMatrix &other) const;



    /** @brief Returns size of this subsitution matrix.
     * @return Size of this matrix
     */
    size_t getSize() const;

    /** @brief Returns number of elements in this matrix.
     * @return Number of scores in this matrix.
     */
    size_t getScoresNumber() const;

    /** @brief Returns residues in this substitution matrix.
     * @return Residues
     */
    string getResidues() const;

    /** @brief Returns score between given residues.
     * @param[in]       a         First amino acid
     * @param[in]       b         Second amino acid
     * @return Score between a and b
     */
    int getScore(const char a, const char b) const;

    /** @brief Returns minimum score in this matrix.
     * @return Minimum score
     */
    int getMinScore() const;

    /** @brief Returns maximum score in this matrix.
     * @return Maximum score
     */
    int getMaxScore() const;

    /** @brief Returns average mismatch score.
     * @return Average mismatch score
     */
    double getAvgMismatchScore() const;

    /** @brief Returns scores as a map.
     * @return Scores as a map
     */
    map<char, map<char, int>> &getScoresAsMap() const;

    /** @brief Returns a SubMatrix reference.
     *
     * Converts this substitution matrix to a SubMatrix object.
     * @return This substitution matrix as a SubMatrix
     */
    SubMatrix &asSubMatrix() const;



    /** @brief Computes sum of scores and a given constant.
     * @param[in]       value     Value to be added
     * @return A new substitution matrix
     */
    SubstitutionMatrix &add(const int value) const;

    /** @brief Computes product of scores and a given constant.
     * @param[in]       value     Value to be multiplied
     * @return A new substitution matrix
     */
    SubstitutionMatrix &multiply(const int value) const;

    /** @brief Performs a fused multiply-add (fma) operation.
     *
     * Multiplies this matrix by given factor, then adds given addend to
     * the result at once. This is more efficient than a + b * M.
     * @param[in]       addend    Value to be added
     * @param[in]       factor    Value to be multiplied
     * @return A new substitution matrix
     */
    SubstitutionMatrix &fma(const int addend, const int factor) const;



 private:
    static const string residues;      ///< Residues

    int scores[26 * 26];         ///< Scores
    double avg_mismatch_score;   ///< Average mismatch score
    int min_score;               ///< Minimum score
    int max_score;               ///< Maximum score


    /** @brief Converts a string into an identifier.
     * 
     * If conversion fails, returns the identity identifier.
     * @param[in]       identifier Identifier as a string
     * @return Identifier
     */
    static Identifier stringToIdentifier(const string &identifier);

    /** @brief Returns an hased version of an amino acid.
     * 
     * Hash function maps amino acids into the interval [0, 25], minimizing
     * entropy.
     * @param[in]       amino     Amino acid
     * @return Hash value
     */
    static size_t hash(const char amino);

    /** @brief Returns a SubMatrix as a map.
     *
     * Converts a SubMatrix object into a map of scores.
     * @param[in]       submatrix SubMatrix object
     * @return Score map
     */
    static map<char, map<char, int>> &SubMatrixToMap(const SubMatrix &submatrix);
};



}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_SUBSTITUTIONMATRIX_H_