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

#ifndef _VICTOR_PHYLO_FENGDoolittleDISTANCE_H_
#define _VICTOR_PHYLO_FENGDoolittleDISTANCE_H_

#define _GLIBCXX_IOMANIP
#include <DistanceMatrixBuilder.h>
#include <SubMatrix.h>
#include <GapFunction.h>

namespace Victor {
namespace Phylo {

/** @brief Builds a distance matrix using the Feng-Doolittle algorithm.
 *
 * Compute pairwise alignment for every pair of sequences, then creates
 * the distance matrix using the formula proposed by Feng and Doolittle:
 * D(x, y) = -log( ((S(x, y) - Srand) / (S(x, x)/2 + S(y, y)/2 - Srand) ).
 *
 * This class follows the Strategy Design Pattern.
 * 
 * @author Marco Zanella <marco.zanella.9@studenti.unipd.it>
 */
class FengDoolittleDistance: public DistanceMatrixBuilder {
 public:
 	/** @brief Builds a Feng-Doolittle distance metrics with given paramters.
 	 *
 	 * Stores information about wich subsitution matrix to use (example:
     * PAM260, BLOSUM50...) and which gap function to use.
 	 * @param[in]       matrix    Substitution matrix to use
 	 * @param[in]       gap       Gap function to use
 	 */
 	explicit FengDoolittleDistance(SubMatrix &matrix, GapFunction &gap);

    /** @brief Copy constructor.
     * @param[in]       origin    Original Feng-Doolittle object
     */
    explicit FengDoolittleDistance(const FengDoolittleDistance &origin);

    /** @brief Destructor. */
    ~FengDoolittleDistance();



 private:
 	SubMatrix   &substitution_matrix;   ///< Substitution matrix
 	GapFunction &gap_function;          ///< Gap function


    /** @brief Returns distance between two sequences.
     * @param[in]       seq1      First sequence
     * @param[in]       seq2      Second sequence
     * @return Distance between seq1 and seq2
     * @todo S_r is not computed due to complete lack of documentation
     */
    virtual double computeDistance(const string &seq1, const string &seq2) const;



    /** @brief Performs alignment and returns score.
     * @param[in]       seq1      First sequence
     * @param[in]       seq2      Second sequence
     * @return Score of the alignment
     */
    double computeScore(const string &seq1, const string &seq2) const;
};

}  // namespace Phylo
}  // namespace Victor

#endif  // _VICTOR_PHYLO_FENGDoolittleDISTANCE_H_